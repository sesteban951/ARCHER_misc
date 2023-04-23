#include "../inc/QP.h"
#include <cassert>

// Use (void) to silence unused warnings.                                                            
#define assertm(exp, msg) assert(((void)msg, exp))

// MAIN QP Program
int QP::solve(Hopper hopper, vector_t &sol, vector_t u_ff) {
  /*
   * Evaluate the gradients of V and h
   * Eval dynamics through pinocchio through Hoppper class
   * compute Lie derivatives
   * eval V and h at the current state
   * fill matrix A and vector b_upper
   * solve the QP ("easy" just do .solve or something)
   * store the solution in some useful data type
  */

  // decision variables
  vector_3t u_opt;
  scalar_t  d_opt;
  u_opt.setZero();
  d_opt.setZero();

  // solution variable
  vector_4t sol_opt(0);
  sol_opt.setZero();

  // current state
  vector_t x(21); // x in Manifold
  x << hopper.q, hopper.v;

  vector_t s(20); // xi in Lie Algebra
  s = Log(x);

  vector_t xi(3), omega(3), eta(6); // somehow need to feed x_d
  xi << s(3), s(4), s(5);
  omega << x(14), x(15), x(16);
  eta << xi(0), xi(1), xi(2), omega(0), omega(1), omega(2); // not good enough, need \ddot qd

  // get f(x) and g(x) matrices
  vector_t f_x(20);
  matrix_t g_x(20,4);
  
  vector_t a(10);
  a.setZero();
  
  f_x = Hopper::f(hopper.q, hopper.v, a, hopper.dom); // full f(x) matrix
  g_x = Hopper::g(hopper.q); // full g(x) matrix (20x4)
  
  // to help compute Lie derivatives
  vector_t f_hat(6);
  matrix_t g_hat(6,3);           // i'm ommiting the input to the spring 
  matrix_t g_x_block(3,3);       // g_x submatrix
  f_hat << omega, f_x.segment(13,15); 
  g_x_block = g_x.block<3,3>(13,0);          // 3x3 block starting at (13,0) index
  g_hat << matrix_t::Zero(3,3), g_x_block    // stack zeros matrix with g_x submatrix
  
  // Lyapunov function
  matrix_6t P = lyap.P_lyap; // CTLE matrices
  matrix_6t Q = lyap.Q_lyap;
  lambda = lyap.lambd;
  V = eta.transpose() * P * eta; // V(eta)

  LfV = f_hat.transpose() * P * eta + eta.transpose() * P * f_hat; // Lie derivatives
  LgV = 2*eta.transpose() * P * g_hat;

  // Barrier function
  vector_t max_fly_vel(3), min_fly_vel(3);
  max_fly_vel << params.flywheel_max_vel, 
                 params.flywheel_max_vel,
                 params.flywheel_max_vel;
  min_fly_vel << params.flywheel_min_vel, 
                 params.flywheel_min_vel,
                 params.flywheel_min_vel;

  h1 << max_fly_vel(1) - x(17),
        max_fly_vel(2) - x(18),
        max_fly_vel(3) - x(19);

  h2 <<  x(17) - min_fly_vel(1),
         x(18) - min_fly_vel(2),
         x(19) - min_fly_vel(3);
  
  alpha1 = barr.alph1;
  alpha2 = barr.alph2;

  Lfh2 = f_x.segment(17,19);
  Lgh2 = g_x_block; 

  Lfh1 = -Lfh2;
  Lgh1 = -Lgh2;

  // Gradeints and hessians (fixed values)
  // OSQP: Cost = 1/2 x^T H x + F^T x -- https://robotology.github.io/osqp-eigen/md_pages_mpc.html
  // Hessian = H and Gradient = F; I dont know why F is called the gradeint
  hess = H;
  grad = F;

  // iterate in SQP fashion;
  for (int iter=0; iter<params.QP_SQP_iter; iter++) {

    updateConstraints(scalar_t lambda, scalar_t V, scalar_t LfV, matrix_t LgV,  
                      scalar_t alpha1, vector_t h1, vector_t Lfh1, matrix_t Lgh1,
                      scalar_t alpha2, vector_t h2, vector_t Lfh2, matrix_t Lgh2);
    solver.updateHessian(hess); // Dont need this w/ static H
    solver.updateGradient(grad);
    solver.updateLinearConstraintsMatrix(A);
    solver.updateBounds(b_lower, b_upper);
    
    // solve the QP problem
    solver.solve();
    sol = solver.getSolution();

    if (iter < params.QP_SQP_iter-1) {
      u_opt << sol.segment(0,2);
      d_opt << sol.segment(3);
      sol_opt << u_opt, d_opt; // insert optimal solution into vector
    }
  }
  return  0;
}

// reset all vars to zero 
void QP::reset() {
   
    LfV.setZero() // stability
    V.setZero();
    lambda.setZero();
    LgV.setZero();
    
    Lfh1.setZero(); // safety 1
    h1.setZero();
    alpha1.setZero();
    Lgh1.setZero();
    
    Lfh2.setZero(); // safety 2
    h2.setZero();
    alpha2.setZero();
    Lgh2.setZero();

    A.setZero();   // constraints as polytope 
    b_lower.setZero();
    b_upper.setZero();

    u.setZero();   // decision variables
    delta.setZero();
    u_bar.setZero();

    H.setZero();  // OSQP cost function matrices
    F.setZero();

    hess.setZero(); // hessian and gradient
    grad.setZero();
}

// build the cost function, onyl done once!
void QP::buildCost(vector_3t u_ff) { // need to pass in MPC u_ff

    int input_size = QP_inputScaling.rows(); 
    int delta_size = QP_deltaScaling.rows(); 

    // build W
    matrix_t W(input_size, input_size);
    for (int i=0; i<input_size; i++) {
      W(i,i) = QP_inputScaling(i);
    }

    // build -2*W*u_ff
    vector_t temp_vec(input_size);
    temp_vec = -2*W*u_ff;

    // populate H, daigonal matrix; H = [2*W, 0; 0, 2*rho]
    for (int i=0; i <  input_size-1; i++) {
      H.insert(i,i) = 2*params.QP_inputScaling(i); // input scaling
    }

    H.insert(input_size-1, input_size-1) = 2*params.QP_deltaScaling;

    // populate F, vector; F = [-2*W*uff; 0]
    for (int i=0; i < input_size-1; i++) {
      F(i) = temp_vec(i); // input scaling
    }

    F(input_size-1) = 0;

}

// build constraint list --  init
void QP::buildConstraints() { // use A.insert() -- follow Noel code
    
    // A matrix is not Sparse at all but need Sparse type to use w/ OSQP-EIGEN
    // fill A with random crap for the init
    int rows = A.rows();
    int cols = A.cols();

    for (int i=0; i<rows; i++) {
      for (int j=0; j<cols; j++) {
         A.insert(i,j) = 420; // arbitrarily fille with zeros?
      }
    }
    
    // fill b upper bounds with random crap
    int len_b = b_upper.rows();
    for (int i=0; i<len_b; i++) {
      b_upper(i) = 420;
    }

    // fill b lower bounds with infinity
    for (int i=0; i<len_b; i++) {
      b_lower(i) = -OsqpEigen::INFTY // or enter big ass number
    }

}

// update the constraint list
void QP::updateConstraints(scalar_t lambda, scalar_t V, scalar_t LfV, matrix_t LgV,  
                           scalar_t alph1, vector_t h1, vector_t Lfh1, matrix_t Lgh1,
                           scalar_t alph2, vector_t h2, vector_t Lfh2, matrix_t Lgh2) { 

  // use A.coeffRef()

  // indexing offset
  int r_offset = 0;
  int c_offset = 0;

  // LgV
  int n_rows = LgV.rows();
  int n_cols = LgV.cols();
  for (int i=0; i<n_rows; i++) {
    for (int j=0; j<n_cols; j++) {
      A.coeffRef(i+r_offset, j+c_offset) = LgV(i,j);
    }
  }

  // -Lgh1
  r_offset += LgV.rows();
  n_rows = Lgh1.rows();
  n_cols = Lgh1.cols();
  for (int i=0; i<n_rows; i++) {
    for (int j=0; j<n_cols; j++) {
      A.coeffRef(i+r_offset, j+c_offset) = -Lgh1(i,j);
    }
  }

  // -Lgh2
  r_offset += Lgh1.rows();
  n_rows = Lgh2.rows();
  n_cols = Lgh2.cols();
  for (int i=0; i<n_rows; i++) {
    for (int j=0; j<n_cols; j++) {
      A.coeffRef(i+r_offset, j+c_offset) = -Lgh2(i,j);
    }

  // scalar 1 for the delta
  r_offset = 0;
  c_offset += LgV.cols();
  A.coeffRef(r_offset,c_offset) = -1;

  // 0's for absence of delta in safety constraints
  r_offset += 1;
  n_rows = Lgh1.rows() + Lgh2.rows();
  for (int i=0; i<n_rows; i++) {
      A.coeffRef(i+r_offset, c_offset) = 0;
    }

  // update b_upper
  r_offset = 0;

  b_upper(r_offset) = -lambda * V - LfV;
  
  r_offset += LfV.rows();
  n_rows = h1.rows();
  for (i=0; i<n_rows; i++) {
    b_upper(i+r_offset) = alpha1 * h1(i) + Lfh1(i);  
  }

  r_offset += h1.rows();
  n_rows = h2.rows();
  for (i=0; i<n_rows; i++) {
    b_upper(i+r_offset) = alpha2 * h2(i) + Lfh2(i);  
  }
}

// Log operator on manifold element
vector_t QP::Log(vector_t x) { 
  
  vector_t g_frak(20); // Lie algebra variable

  quat_t quat(x(6), x(3), x(4), x(5));  // load in quaternion state
  
  auto quat_ = manif::SO3<scalar_t>(quat); 
  manif::SO3Tangent<scalar_t> xi = quat_.log(); // Manifold to Lie algebra
  
  g_frak << x.segment(0,3),xi.coeffs(),x.segment(7,4),x.segment(11,10); // full state as Lie algebra
  
  return g_frak; 
}

// Exp operator on Lie algebra element
vector_t QP::Exp(vector_t xi) {

  vector_t g(21); // manifold variable
  
  manif::SO3Tangent<scalar_t> xi_; // Lie algebra type
  xi_ << xi(3),xi(4),xi(5);        // populate Lie algebra var
  quat_t quat = xi_.exp().quat();  // Lie algerba to manifold
  
  g << xi.segment(0,3), quat.coeffs(), xi.segment(6,14); // full state on Manifold
  return g;
}




