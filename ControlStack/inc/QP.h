#include <Eigen/Dense>                                                                               
#include "OsqpEigen/OsqpEigen.h"
#include <iostream>
#include "Types.h"
#include "yaml-cpp/yaml.h"
#include "Hopper.h"
#include <manif/manif.h>

using namespace Hopper_t;
using namespace Eigen;

//  QP solver class for CLF-CBF (ONLY on the flywheels)
class QP {

    public:

        int nu; // dim of input (3 vector)
        int nd; // dim of delta (scalar)

        scalar_t LfV, V, lambda; // Stability criteria
        vector_t LgV; // 3x1

        scalar_t alph1; // Safety criteria 1
        vector_t h1;
        vector_t Lfh1; // 3x1
        vector_t Lgh1; // 3x3

        scalar_t alph2; // Safety criteria 2
        vector_t h2;
        vector_t Lfh2; // 3x1
        vector_t Lgh2; // 3x3 

        vector_t u_bar; // minimize [u; delta]
        vector_t u;
        scalar_t delta;

        // Safety plus stability criteria as polytope
        Eigen::SparseMatrix<double,ColMajor> A;    // to write cons as polytope
        vector_t b_lower, b_upper;                // to write cons as polytope

        // for cost function
        Eigen::SparseMatrix<double,ColMajor> SparseIdentity;
        Eigen::SparseMatrix<double,ColMajor> H;       // Hessian, quadratic term of cost
        vector_t F;                                   // gradient, linear term in cost

        // OSQP solver object
        OsqpEigen::Solver solver;

        // Parameters for the QP program
        struct QP_Params {
            
            int QP_SQP_iter;    // SQP iterations

            // argmin 0.5 u^T H u + F^T u s.t. stability1, safety1, safety2
            vector_t QP_quad_scaling; // input scaling, diagonal entries of H
            vector_t QP_lin_scaling;  // delta scaling, diagonal entry of F

            scalar_t u_max;  // upper input limit
            scalar_t u_min;  // lower input limit

            scalar_t flywheel_max; // max flywheel speed
            scalar_t flywheel_min; // min flyhweel speed

        } params;

        // Initialize the QP object
        QP (int nu, int nd, QP_Params &loaded_params) {
            
            // var dimension
            this -> nu = nu;
            this -> nd = nd;

            // load in the QP parameters (sent by YAML)
            params = loaded_params;

            // QP solver params
            int nvars = nu + nd; // u in R^3, delta in R
            bool warm_start = true;
            bool verbose = false;
            double tol = 1e-9;

            // QP messages
            std::cout << "QP Settings:" << std::endl;
            std::cout << "SQP Iterations: " << params.QP_SQP_iter << std::endl;

            // resize all QP variables. Or prespecify them above
            LgV.resize(1,nu); 
            
            h1.resize(3,1);
            Lfh1.resize(3,1);  // hardcoded for barrier 
            Lgh1.resize(3,nu);
            
            h2.resize(3,1);
            Lfh2.resize(3,1);
            Lgh2.resize(3,nu);

            u_bar.resize(num_vars,1);
            u.resize(nu,1);

            int A_rows = LgV.rows() + Lgh1.rows() + Lgh2.rows();
            int A_cols = nvars;

            A.resize(A_rows, A_cols);
            b_lower.resize(A_rows);
            b_upper.resize(A_rows);
            SparseIdentity.resize(A_rows,A_cols);
            H.resize(nvars,nvars);
            F.resize(nvars,1);

            // solver settings
            solver.settings() -> setWarmStart(warm_start);
            solver.settings() -> setVerbosity(verbose);
            //solver.settings()->setAbsoluteTolerance(tol); // hardcode the nuber of its
            solver.data() -> setNumberOfVariables(nvars);     // [u; delta];
            solver.data() -> setNumberOfConstraints(b_lower.rows());   // stability1; safety2l; safety2;
            
            reset();
            buildCost();
            buildConstraints(); // no updateConstraints?

            solver.data() -> setHessianMatrix(H); 
            solver.data() -> setGradient(F);
            solver.data() -> setLinearConstraintMatrix(A);
            solver.data() -> setLowerBound(b_lower);
            solver.data() -> setUpperBound(b_upper);

            // instantiate the solver
            solver.initSolver();
        };

        static vector_t Log(vector_t x):        // take manifold element to Lie algebra
        static vector_t Exp(vector_t xi);       // take Lie element to manifold
        
        void reset();    // reset all internal variables

        // QP solver
        int solve(Hopper hopper, vector_t &sol) //vector_3t &command, vector_2t &command_interp);

        void buildCost();         // build cost function
        void buildConstraints();  // build constraint list
        void updateConstraints(); // update constraints
};

