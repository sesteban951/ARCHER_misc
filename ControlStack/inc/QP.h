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
        int ncbf; //dim of CBF image

        scalar_t LfV, V, lambda; // Stability criteria
        matrix_t LgV; // 1x3

        scalar_t alpha1; // Safety criteria 1
        vector_t h1;   // 3x1
        vector_t Lfh1; // 3x1
        matrix_t Lgh1; // 3x3

        scalar_t alpha2; // Safety criteria 2
        vector_t h2;   // 3x1
        vector_t Lfh2; // 3x1
        matrix_t Lgh2; // 3x3 

        // Safety plus stability criteria as polytope
        Eigen::SparseMatrix<double,RowMajor> A;    // to write cons as polytope
        vector_t b_lower, b_upper;                // to write cons as polytope

        // OSQP cist function matrices
        Eigen::SparseMatrix<double,RowMajor> H;       // quadratic weights
        vector_t F;                                   // linear weights

        // gradient and Hessian
        Eigen::SparseMatrix<double,RowMajor> hess;       // Hessian, quadratic term of cost
        vector_t grad;                                   // gradient, linear term in cost

        // OSQP solver object
        OsqpEigen::Solver solver;

        // Parameters for the QP program
        struct QP_Params {
            
            int QP_SQP_iter;    // SQP iterations

            // argmin 0.5 u^T H u + F^T u s.t. stability1, safety1, safety2
            vector_3t QP_inputScaling; // input scaling, diagonal entries of W
            scalar_t  QP_deltaScaling;  // delta scaling, rho

            scalar_t u_max;  // upper input limit
            scalar_t u_min;  // lower input limit

            scalar_t flywheel_max_vel; // max flywheel speed
            scalar_t flywheel_min_vel; // min flyhweel speed

        } params;

        // Initialize the QP object
        QP (int nu, int nd, int ncbf, QP_Params &loaded_params) {
            
            // var dimension
            this -> nu = nu;
            this -> nd = nd;
            this -> ncbf = ncbf;   

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
            
            h1.resize(ncbf,1);
            Lfh1.resize(ncbf,1);  // hardcoded for barrier 
            Lgh1.resize(ncbf,nu);
            
            h2.resize(ncbf,1);
            Lfh2.resize(ncbf,1);
            Lgh2.resize(ncbf,nu);

            int A_rows = LgV.rows() + Lgh1.rows() + Lgh2.rows();
            int A_cols = nvars;

            A.resize(A_rows, A_cols);
            b_lower.resize(A_rows);
            b_upper.resize(A_rows);
            H.resize(nvars,nvars);
            F.resize(nvars,1);
            hess.resize(nvars, nvars);
            grad.resize(nvars,1);

            // solver settings
            solver.settings() -> setWarmStart(warm_start);
            solver.settings() -> setVerbosity(verbose);
            //solver.settings()->setAbsoluteTolerance(tol); // hardcode the nuber of its
            solver.data() -> setNumberOfVariables(nvars);     // [u; delta];
            solver.data() -> setNumberOfConstraints(b_lower.rows());   // stability1; safety2l; safety2;
            
            reset();
            // At some point, need to pass in u_ff into this -----------------------------------------
            vector_3t u_ff;
            u_ff << 0,0,0;
            buildCost(u_ff);
            buildConstraints();

            solver.data() -> setHessianMatrix(hess); 
            solver.data() -> setGradient(grad);
            solver.data() -> setLinearConstraintsMatrix(A);
            solver.data() -> setLowerBound(b_lower);
            solver.data() -> setUpperBound(b_upper);

            // instantiate the solver
            solver.initSolver();
        };

        void reset();    // reset all internal variables
        void buildCost(vector_t u_ff);         // build cost function
        void buildConstraints();  // build constraint list
        void updateConstraints(); // update constraints
        
        static vector_t Log(vector_t x); // Log(.) : Manif -> Lie
        static vector_t Exp(vector_t xi); // Exp(.) : LIe -> Manif
        
        // QP solver
        int solve(Hopper hopper, vector_t &sol, vector_t u_ff); //vector_3t &command, vector_2t &command_interp);
};

