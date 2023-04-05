#include "../inc/QP.h"
#include <cassert>

// Use (void) to silence unused warnings.                                                            
#define assertm(exp, msg) assert(((void)msg, exp))

// MAIN QP Program
int QP::solve(Hopper hopper, vector_t &sol) {

  // decision variables
  vector_3t = u_star;
  scalar_t = delta_star;

  // current state
  vector_t x_(21);
  x_ << hopper.q, hopper.v;

  // I need to access dynamic info somehow.
  // Need some lines from f(x) and g(x) to construct my Lie derivatives



  return  0;
}

// reset all vars to zero 
void QP::reset() {
   
    LfV.setZero()
    V.setZero();
    lambda.setZero();
    delta.setZero();
    
    Lfh1.setZero();
    h1.setZero();
    alph1.setZero();
    Lgh1.setZero();
    
    Lfh2.setZero();
    h2.setZero();
    alph2.setZero();
    Lgh2.setZero();

    A.setZero();
    b_lower.setZero();
    b_upper.setZero();

    H.setZero();
    F.setZero();

}

// build the cost function, onyl done once!
void QP::buildCost() {

    int quad_size = QP_quad_scaling.rows() - 1; // last index in vector
    int lin_size = QP_lin_scaling.rows() - 1;   // last index in vector

    // populate H, daigonal matrix
    for (int i=0; i <  quad_size; i++) {
      H.insert(i,i) = params.QP_quad_scaling(i); // input scaling
    }

    H.insert(quad_size, quad_size) = params.QP_quad_scaling(quad_size);

    // populate F, vector
    for (int i=0; i < lin_size; i++) {
      F.(i) = params.QP_input_scaling(i); // input scaling
    }

    F(lin_size) = params.QP_input_scaling(lin_size);

}

// build constraint list
void QP::buildConstraints() {



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




