#include "../inc/Hopper.h"
#include "../inc/MPC.h"
//#include "../inc/QP.h"
#include <stdexcept>

#include <manif/manif.h>

#include "pinocchio/algorithm/cholesky.hpp"
#include "pinocchio/algorithm/joint-configuration.hpp"
#include "pinocchio/algorithm/aba.hpp"
#include "pinocchio/algorithm/aba-derivatives.hpp"
#include "pinocchio/algorithm/constrained-dynamics-derivatives.hpp"
#include "pinocchio/algorithm/impulse-dynamics.hpp"
#include "pinocchio/algorithm/impulse-dynamics-derivatives.hpp"
#include "pinocchio/algorithm/constrained-dynamics.hpp"


// added
#include <iostream>   
using namespace std;

using namespace pinocchio;

matrix_3t Hopper::cross(vector_3t q) {
    matrix_3t c;
    c << 0, -q(2), q(1),
            q(2), 0, -q(0),
            -q(1), q(0), 0;
    return c;
}

matrix_3t Hopper::quat2Rot(quat_t q) {
    scalar_t qw,qx,qy,qz;
    matrix_3t Rq;
    qw = q.w();
    qx = q.x();
    qy = q.y();
    qz = q.z();
    Rq << pow(qw,2)+pow(qx,2)-pow(qy,2)-pow(qz,2), 2*(qx*qy-qw*qz), 2*(qx*qz+qw*qy),
       2*(qx*qy+qw*qz), pow(qw,2)-pow(qx,2)+pow(qy,2)-pow(qz,2),2*(qy*qz-qw*qx),
       2*(qx*qz-qw*qy), 2*(qy*qz+qw*qx), pow(qw,2)-pow(qx,2)-pow(qy,2)+pow(qz,2);
    return Rq;
}

Hopper::Hopper() {
    // Construct Pinocchio model
    const std::string urdf_path_c = "../rsc/hopper.urdf";
    model = pinocchio::Model();
    pinocchio::urdf::buildModel(urdf_path_c, pinocchio::JointModelFreeFlyer(), model);
    data = Data(model);

    contact_model_ground.emplace_back(RigidConstraintModelTpl<scalar_t, 0>(CONTACT_3D, 2, SE3::Identity(), LOCAL_WORLD_ALIGNED));
    contact_data_ground.emplace_back(RigidConstraintDataTpl<scalar_t, 0>(contact_model_ground.at(0)));
    initConstraintDynamics(model, data, contact_model_ground);

    // Contact for when foot hits hard stops going to flight phase
    contact_model_flight.emplace_back(RigidConstraintModelTpl<scalar_t, 0>(CONTACT_3D, 2,LOCAL));
    contact_data_flight.emplace_back(RigidConstraintDataTpl<scalar_t, 0>(contact_model_flight.at(0)));
    initConstraintDynamics(model, data, contact_model_flight);

    // Read gain yaml
    YAML::Node config = YAML::LoadFile("../config/gains.yaml");
    std::vector<scalar_t> orientation_kp = config["LowLevel"]["Orientation"]["Kp"].as<std::vector<scalar_t>>();
    std::vector<scalar_t> orientation_kd = config["LowLevel"]["Orientation"]["Kd"].as<std::vector<scalar_t>>();
    gains.leg_kp = config["LowLevel"]["Leg"]["Kp"].as<scalar_t>();
    gains.leg_kd = config["LowLevel"]["Leg"]["Kd"].as<scalar_t>();

    gains.orientation_kp << orientation_kp[0], orientation_kp[1], orientation_kp[2];
    gains.orientation_kd << orientation_kd[0], orientation_kd[1], orientation_kd[2];

    springStiffness = config["SpringStiffness"].as<scalar_t>();

    quat_actuator = quat_t(0.8806, 0.3646, -0.2795, 0.1160);
    multiplier_on_deltaf = config["MPC"]["multiplier_on_deltaf"].as<scalar_t>();

    // size of state vector
    q.resize(model.nq);
    v.resize(model.nv);
    
    // get CLF tuning knobs
    int square_dim = 6;
    int vec_idx = 0;
    std::vector<scalar_t> temp;
    
    temp = config["CLF_CBF"]["CLF"]["P"].as<std::vector<scalar_t>>();
    for (int i=0; i<square_dim; i++) {
      for (int j=0; j<square_dim; j++) {
        lyap.P_lyap(i,j) = temp[vec_idx];   // P matrix
        vec_idx++;
      }
    }
    
    vec_idx = 0;
    temp = config["CLF_CBF"]["CLF"]["Q"].as<std::vector<scalar_t>>();
    for (int i=0; i<square_dim; i++) {
      for (int j=0; j<square_dim; j++) {
        lyap.Q_lyap(i,j) = temp[vec_idx]; // Q matrix
        vec_idx++;
      }
    }

    lyap.lambd = config["CLF_CBF"]["CLF"]["lambd"].as<scalar_t>(); // lambda

    // get CBF tuning knobs
    barr.alph1 = config["CLF_CBF"]["CBF"]["alph1"].as<scalar_t>(); // linear extended class K inf
    barr.alph2 = config["CLF_CBF"]["CBF"]["alph2"].as<scalar_t>(); // linear extended class K inf
}

void Hopper::updateState(vector_t state) {
    int ind = 0;

    t = state[ind];
    ind++;
    pos << state[ind], state[ind + 1], state[ind + 2];
    ind += 3;
    quat = Eigen::Quaternion<scalar_t>(state[ind], state[ind + 1], state[ind + 2], state[ind + 3]); // Loaded as w,x,y,z
    ind += 4;
    vel << state[ind], state[ind + 1], state[ind + 2];
    ind += 3;
    omega << state[ind], state[ind + 1], state[ind + 2];
    ind += 3;
    if (contact == 0 && state[ind] >= .1) {
	    last_impact_time = t;
    }
    if (contact == 1 && state[ind] <= .1) {
	    last_flight_time = t;
   }
    if (state[ind] <= .1) {
            contact = 0;
    } else {
            contact = 1;
    }
    ind++;
    leg_pos = state[ind];
    ind++;
    leg_vel = state[ind];
    ind++;
    wheel_vel << state[ind], state[ind + 1], state[ind + 2];
    ind++;
    q << pos, quat.coeffs(), leg_pos, 0,0,0;
    v << vel, omega, leg_vel, wheel_vel;
};

//////////////////////////////////////////////////////////////////////////////////////////////

 // Noel's PD controller on quaternion error
 void Hopper::computeTorque(quat_t quat_d_, vector_3t omega_d, scalar_t length_des, vector_t u_des) {
       
    vector_3t omega_a;

    quat_t quat_a_ = quat; // "quat" is variable in Hopper class

    vector_4t quat_d;  // convert to vector type to make it easier to work with
    vector_4t quat_a;
    quat_d << quat_d_.w(), quat_d_.x(), quat_d_.y(), quat_d_.z(); // "quat_d_" passed in
    quat_a << quat_a_.w(), quat_a_.x(), quat_a_.y(), quat_a_.z();

    vector_3t delta_quat;   //quaternion error, Noel does q_a^-1 * q_d
    delta_quat << quat_a[0] * quat_d.segment(1, 3) - quat_d[0] * quat_a.segment(1, 3) -
                  cross(quat_a.segment(1, 3)) * quat_d.segment(1, 3);
 
    matrix_3t Kp, Kd; 
    Kp.setZero();
    Kd.setZero();
    Kp.diagonal() << gains.orientation_kp; // "gains" is struct from Hopper class
    Kd.diagonal() << gains.orientation_kd;

    vector_3t tau;
    omega_a = -quat_actuator.inverse()._transformVector(omega); // omega from Hopper class. 
    tau = quat_actuator.inverse()._transformVector(Kp * delta_quat) - Kd * (omega_a - omega_d);

    scalar_t spring_f = (1 - contact) * (-gains.leg_kp * (leg_pos - length_des) - gains.leg_kd * leg_vel); // contact from Hopper Class
    torque << spring_f, tau;
    torque += u_des;
};

// Lyapunov based controller on xi and omega output functions
void Hopper::computeTorqueQP(quat_t quat_d_, vector_3t omega_d, scalar_t length_des, vector_t u_des, QP::QP_Params qp_params) {
    
    // create quaternion variables
    quat_t    quat_a_ = quat;   // quat is Hopper class variable
    vector_4t quat_d;           // quaternion desired
    vector_4t quat_a;           // quaternion actual
    vector_4t quat_e;           // quaternion error

    // create omega variable
    vector_3t omega_a;

    // quaternion error
    quat_d << quat_d_.w(), quat_d_.x(), quat_d_.y(), quat_d_.z(); // as vector
    quat_a << quat_a_.w(), quat_a_.x(), quat_a_.y(), quat_a_.z(); // as vector
    
    scalar_t  sd, sa;
    vector_3t vd, va;
    sd = quat_d(0);
    vd = quat_d.segment(1,3);
    sa = quat_a(0);
    va = quat_a.segment(1,3);
    quat_e << sa*sd + va.transpose() * vd,  // q_e = q_a^-1 * q_d
              sa*vd - sd*va - cross(va)*vd;
    
    // no idea what's going on here, but OK
    omega_a = -quat_actuator.inverse()._transformVector(omega); // omega from Hopper class. 

    // create QP object -- WARNING: recreating the QP onbject could be slow
    int nu, nd, ncbf;
    nu = 3;
    nd = 1;
    ncbf = 3;    
    QP QP_obj = QP(nu, nd, ncbf, qp_params, u_des.segment(1,3));

    // solve QP problem
    vector_t u_opt(nu);
    scalar_t d_opt;
    vector_t sol(nu+nd);
    QP_obj.solve(this, sol, u_des.segment(1,3), quat_e, omega_d);
    u_opt = sol.segment(0,3);
    d_opt = sol(3);

    vector_t tau(3);
    tau = u_opt;

    // update torque
    scalar_t spring_f = (1 - contact) * (-gains.leg_kp * (leg_pos - length_des) - gains.leg_kd * leg_vel); // contact from Hopper Class
    torque << spring_f, tau;
    torque += u_des;
};

//////////////////////////////////////////////////////////////////////////////////////////////

vector_t Hopper::f(const vector_t& q, const vector_t& v, const vector_t& a, const domain& d) {
    vector_t x_dot(2*model.nv);
    switch(d)
    {
        case flight: {
		//Not constrained dynamics for flight, because that would fix the foot position
		aba(model, data, q, v, a);
    		x_dot << v.segment(0,6),0,v.segment(7,3), data.ddq;
		break;
	}
        case ground: {
                initConstraintDynamics(model, data, contact_model_ground);
		constraintDynamics(model, data, q, v, a, contact_model_ground, contact_data_ground);
		vector_t springForce(10);
		springForce << 0,0,springStiffness*q(7),0,0,0,-springStiffness*q(7),0,0,0;
    		x_dot << v, data.ddq + springForce;
		break;
	}
	otherwise: {
	     throw std::invalid_argument("Invalid domain in f");
	     break;
	}
    }
    quat_t quat(q(6), q(3), q(4), q(5));
    x_dot.segment(0,3) += Hopper::cross(q.segment(0,3))*quat.inverse()._transformVector(v.segment(3,3));
    return x_dot;
};

void Hopper::Df(const vector_t q, const vector_t v, const vector_t a, const domain d,
                matrix_t &A, matrix_t &B, matrix_t &C, const vector_t q0) {
	matrix_t springJacobian(10,10);
	// Call dynamics first to populate q, v, a, for ComputeConstraintDynamicsDerivatives
	vector_t f = Hopper::f(q, v, a, d);
    switch(d)
    {
        case flight: {
		//Not constrained dynamics for flight, because that would fix the foot position
		computeABADerivatives(model, data, q,v,a);
		springJacobian.setZero();
		break;
	}
        case ground: {
                initConstraintDynamics(model, data, contact_model_ground);
		computeConstraintDynamicsDerivatives(model, data, contact_model_ground, contact_data_ground);
		springJacobian << matrix_t::Zero(2,10),
		       0,0,0,0,0,0,springStiffness,0,0,0,
		       matrix_t::Zero(3,10),
		       0,0,0,0,0,0,-springStiffness,0,0,0,
		       matrix_t::Zero(3,10);
		break;
        }
	otherwise: {
	     throw std::invalid_argument("Invalid domain in Df");
	     break;
	}
    }
    A << matrix_t::Zero(10,10),matrix_t::Identity(10,10), data.ddq_dq+springJacobian, data.ddq_dv;

    matrix_t B_mat(10,4);
    B_mat << 0,0,0,0,
        0,0,0,0,
        0,0,0,0,
        0,0,0,0,
        0,0,0,0,
        0,0,0,0,
        1,0,0,0,
        0,1,0,0,
        0,0,1,0,
        0,0,0,1;
    B << matrix_t::Zero(10,4), data.Minv*B_mat;

    vector_t s(20);
    vector_t x(21);
    x << q,v;
    s = MPC::qk_to_xik(x,q0);
	
    C << f - A * s - B * a.tail(4);
};

//////////////////////////////////////////////////////////////////////////////////////////////

// function to return the actuation matrix
matrix_t Hopper::g(const vector_t q) {
  
  // right now this might be fine -- but make sure you call the f_x function

  // do I need to query something to get an updated data.Minv?

  // acutation selection matrix for \ddot q
  matrix_t B_mat(10,4);
  B_mat << 0,0,0,0,
        0,0,0,0,
        0,0,0,0,
        0,0,0,0,
        0,0,0,0,
        0,0,0,0,
        1,0,0,0, // spring
        0,1,0,0, // flywheel
        0,0,1,0, // flywheel
I        0,0,0,1; // flywheel ?

  // build g(x)
  matrix_t g_act(20,4);
  g_act = matrix_t::Zero(10,4), data.Minv*B_mat; // do I need to update data to get latest Minv?

  return g_act;

};

//////////////////////////////////////////////////////////////////////////////////////////////

void Hopper::css2dss(const matrix_t &Ac, const matrix_t &Bc, const matrix_t &Cc, const float dt,
                                           matrix_t &Ad, matrix_t &Bd, matrix_t &Cd) {
    // Exact discretization
    int dim = Ac.cols() + Bc.cols() + Cc.cols();
    matrix_t M = matrix_t::Zero(dim, dim);
    M << Ac, Bc, Cc, Eigen::Matrix<scalar_t, Eigen::Dynamic,Eigen::Dynamic>::Zero( Bc.cols() + Cc.cols(), dim);
    matrix_t Me = matrix_t::Zero(dim, dim);
    M *= dt;
    Me = M.exp();
    Ad << Me.block(0, 0, Ac.rows(), Ac.cols());
    Bd << Me.block(0, Ac.cols(), Bc.rows(), Bc.cols());
    Cd << Me.block(0, Ac.cols() + Bc.cols(), Cc.rows(), Cc.cols());

    // One step euler prediction
    //Ad = (matrix_t::Identity(Ac.rows(), Ac.cols()) + Ac*dt);
    //Bd = Bc*dt;
    //Cd = Cc*dt;
};

vector_t Hopper::delta_f(const vector_t q, const vector_t v, const domain d){
    const double mu0 = 0.;
    const double r_coeff = 0; // restitution coeff -- assumes perfectly plastic
    switch(d)
    {
        case flight_ground: {
          initConstraintDynamics(model, data, contact_model_ground);
          impulseDynamics(model, data, q, v, contact_model_ground, contact_data_ground, r_coeff, mu0);
	  break;
        }
	case ground_flight: {
          initConstraintDynamics(model, data, contact_model_flight);
          impulseDynamics(model, data, q, v, contact_model_flight, contact_data_flight, r_coeff, mu0);
	  break;
	}
    }
    vector_t x_plus(21); 
    x_plus << q, (1-multiplier_on_deltaf)*v + multiplier_on_deltaf*data.dq_after;
    return x_plus;
};

void Hopper::Ddelta_f(const vector_t q, const vector_t v, const domain d,
                      matrix_t &A, matrix_t &B, matrix_t &C, const vector_t q0){
    vector_t df(21); 
    switch(d)
    {
        case flight_ground: {
          df = Hopper::delta_f(q, v, d);
          initConstraintDynamics(model, data, contact_model_ground);
          computeImpulseDynamicsDerivatives(model, data, contact_model_ground, contact_data_ground);
	  break;
	}
	case ground_flight: {
          df = Hopper::delta_f(q, v, d);
          initConstraintDynamics(model, data, contact_model_flight);
          computeImpulseDynamicsDerivatives(model, data, contact_model_flight, contact_data_flight);
	  break;
	}
    }
    A << matrix_t::Identity(10,10),matrix_t::Zero(10,10), data.ddq_dq, data.ddq_dv;

    B.setZero();

    vector_t s(20);
    vector_t x(21);
    x << q,v;
    s = MPC::qk_to_xik(x,q0);

    vector_t s_df(20);
    s_df << s.segment(0,10),df.segment(11,10);

    C << s_df - A * s;
};

void Hopper::DiscreteDynamics(const vector_t &x, const vector_t &u, const domain &d, const float dt,
                              matrix_t &Ac, matrix_t &Bc, matrix_t &Cc,
                              matrix_t &Ad, matrix_t &Bd, matrix_t &Cd,
			      const vector_t q0) {
    switch(d)
    {
        case flight: {
            vector_t a(10);
            a << 0, 0, 0, 0, 0, 0, u;
            Df(x.segment(0, 11), x.segment(11, 10), a, d, Ac, Bc, Cc,q0);
            css2dss(Ac, Bc, Cc, dt, Ad, Bd, Cd);
            break;
        }
        case ground: {
            vector_t a(10);
            a << 0, 0, 0, 0, 0, 0, u;
            Df(x.segment(0, 11), x.segment(11, 10), a, d, Ac, Bc, Cc,q0);
            css2dss(Ac, Bc, Cc, dt, Ad, Bd, Cd);
            break;
        }
        case flight_ground: {
            Ddelta_f(x.segment(0, 11), x.segment(11, 10), d, Ad, Bd, Cd,q0);
            break;
        }
        case ground_flight: {
            Ddelta_f(x.segment(0, 11), x.segment(11, 10), d, Ad, Bd, Cd,q0);
            break;
        }
    }
}
