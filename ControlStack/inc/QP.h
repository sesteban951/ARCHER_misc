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

        scalar_t LfV, V, lambda, delta; // Stability criteria
        row_vector_3t LgV;

        scalar_t Lfh1, h1, alph1; // Safety criteria 1
        row_vector_3t Lgh1;

        scalar_t Lfh2, h2, alph2; // Safety criteria 2
        row_vector_3t Lgh2;

        vector_4t u_bar; // minimize [u; delta]
        vector_3t u;

        // Safety plus stability criteria as polytope
        Eigen::SparseMatrix<double,ColMajor> A;                 // to write cons as polytope
        Eigen::SparseMatrix<double,ColMajor> b_lower, b_upper;  // to write cons as polytope

        // for ocst functions
        Eigen::SparseMatrix<double,ColMajor> SparseIdentity;
        Eigen::SparseMatrix<double,ColMajor> H;       // Hessian, quadratic term of cost
        vector_t F;                                   // gradient, linear term in cost

        // OSQP solver object
        Osqp:w
        Eigen::Solver solver;

        // Parameters for the QP program
        struct QP_Params {
            
            int QP_SQP_iter;    // SQP iterations

            // argmin 0.5 u^T H u + F^T u s.t. stability1, safety1, safety2
            vector_t QP_quad_scaling; // input scaling, diagonal entries of H
            vector_t QP_lin_scaling;  // delta scaling, diagonal entry of H

            scalar_t u_max;  // upper input limit
            scalar_t u_min;  // lower input limit

        } params;

        // Initialize the QP object
        QP (QP_Params &loaded_params) {
            
            // load in the QP parameters (sent by YAML)
            params = loaded_params;

            // QP solver params
            int num_vars = 4;
            int num_cons = 3;    
            bool verbose = false;
            double tol = 1e-9;

            // QP messages
            std::cout << "QP Settings:" << std::endl;
            std::cout << "SQP Iterations: " << params.SQP_iter << std::endl;

            // need all QP parameters

            // solver settings
            solver.settings()->setVerbosity(verbose);
            //solver.settings()->setAbsoluteTolerance(tol); // hardcode the nuber of its
            solver.data() -> setNumberOfVariables(num_vars);     // [u; delta];
            solver.data() -> setNumberOfConstraints(num_cons);   // stability1; safety2l; safety2;
            
            reset();
            buildCost();
            buildConstraints();

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

        void buildCost();        // build cost function
        void buildConstraints(); // build constraint list

};

