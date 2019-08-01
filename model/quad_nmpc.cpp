#include <acado_optimal_control.hpp>
#include <acado_toolkit.hpp>
#include <acado_gnuplot.hpp>
#include <matrix_vector/matrix.hpp>
#include <iostream>

using namespace std;
USING_NAMESPACE_ACADO

int main(){

	//variables
	const double t_start = 0.0; //[s]
  	const double t_end = 50;//50.0;
  	const double dt = 0.25;//0.25;
  	const int 	 N = round(t_end/dt);
  	const double f_min = 0;
  	const double f_max = 4;
  	const double omega_min = -20; //dont really care
  	const double omega_max = omega_min*-1;
  	const double theta = 3.14/3;
  	const double phi_init = 0;
  	const double theta_init = 3.14/1.001;
  	const double psi_init = 0;
  	const double c = 0.03793; //prob wrong
  	const double mass = 0.468;
  	const double g = 9.81;


	//Differential states
	DifferentialState 	b1("",3,1);
	DifferentialState 	b2("",3,1);
	DifferentialState 	b3("",3,1);
	DifferentialState 	omega("",3,1);

	//Dummy state to introduce slack var
	DifferentialState dummy;

	//Online data, "references"
	OnlineData  psi_d;
	OnlineData 	omega_d("",3,1);
	OnlineData 	T_d("",3,1);
	OnlineData 	t;
	const int NR_ONLINE_DATA = 8;


	//Control signals
	Control u("",4,1);

	//Slack variable
	Control sv;

	const double I_xx 	= 0.0049; //Roll inertia
	const double I_yy 	= 0.0049; //Pitch inertia
	const double I_zz 	= 0.0088; //Yaw inertia
	const double l 		= 0.225;

	//Intermediate states
	IntermediateState T_phi 	= 	l*(u(3)-u(1));
	IntermediateState T_theta 	= 	l*(u(0)-u(2));
	IntermediateState T_psi 	= 	c*(-u(0)+u(1)-u(2)+u(3));
	
	IntermediateState T_d_norm = sqrt(T_d(0)*T_d(0)+T_d(1)*T_d(1)+T_d(2)*T_d(2));
	
	IntermediateState R_d("",3,3);

	//Caluclations of b1_d
	IntermediateState b1_d("",3,1);
	b1_d(0) = cos(psi_d);
	b1_d(1) = sin(psi_d);
	b1_d(2) = 0;

	//Calculations of b3_d
	IntermediateState b3_d("",3,1);
	
	b3_d(0) = T_d(0)/T_d_norm;
	b3_d(1) = T_d(1)/T_d_norm;
	b3_d(2) = T_d(2)/T_d_norm;
	
	//Calculations of b2_d
	IntermediateState b3xb1("",3,1);
	b3xb1(0) = b3_d(1)*b1_d(2) - b1_d(1)*b3_d(2);
	b3xb1(1) = b1_d(0)*b3_d(2) - b3_d(0)*b1_d(2);
	b3xb1(2) = b3_d(0)*b1_d(1) - b1_d(0)*b3_d(1);
	IntermediateState crossProduct = sqrt(b3xb1(0)*b3xb1(0)+b3xb1(1)*b3xb1(1)+b3xb1(2)*b3xb1(2));
	IntermediateState b2_d("",3,1);
	b2_d(0) = b3xb1(0)/crossProduct;
	b2_d(1) = b3xb1(1)/crossProduct;
	b2_d(2) = b3xb1(2)/crossProduct;

	//Update desired b1_d
	R_d(0,0) = b2_d(1)*b3_d(2) - b3_d(1)*b2_d(2);
	R_d(1,0) = b3_d(0)*b2_d(2) - b2_d(0)*b3_d(2);
	R_d(2,0) = b2_d(0)*b3_d(1) - b3_d(0)*b2_d(1);
	//Update desired b2_d
	R_d(0,1) = b2_d(0);
	R_d(1,1) = b2_d(1);
	R_d(2,1) = b2_d(2);
	//Update desired b3_d
	R_d(0,2) = b3_d(0);
	R_d(1,2) = b3_d(1);
	R_d(2,2) = b3_d(2);

	//Differential equation
	DifferentialEquation f;

	//Differential equation deffinitions
	f	<<	dot(b1(0)) == b2(0)*omega(2)-b3(0)*omega(1);
	f	<<	dot(b1(1)) == b2(1)*omega(2)-b3(1)*omega(1);
	f	<<	dot(b1(2)) == b2(2)*omega(2)-b3(2)*omega(1);

	f	<<	dot(b2(0)) == -b1(0)*omega(2)+b3(0)*omega(0);
	f	<<	dot(b2(1)) == -b1(1)*omega(2)+b3(1)*omega(0);
	f	<<	dot(b2(2)) == -b1(2)*omega(2)+b3(2)*omega(0);

	f	<<	dot(b3(0)) == b1(0)*omega(1)-b2(0)*omega(0);
	f	<<	dot(b3(1)) == b1(1)*omega(1)-b2(1)*omega(0);
	f	<<	dot(b3(2)) == b1(2)*omega(1)-b2(2)*omega(0);

	f	<<	dot(omega(0)) == -1*omega(1)*omega(2)*(I_zz-I_yy)/(I_xx) + l*(u(3)-u(1))/I_xx;
	f	<<	dot(omega(1)) == -1*omega(0)*omega(2)*(I_xx-I_zz)/(I_yy) + l*(u(0)-u(2))/I_yy;
	f	<<	dot(omega(2)) == -1*omega(0)*omega(1)*(I_yy-I_xx)/(I_zz) + c*(-u(0)+u(1)-u(2)+u(3))/I_zz;

	//mock dynamics of slack variable
	f << dot(dummy) == sv;

	//Defining least squares function
	IntermediateState R("",3,3);
	R(0,0) = b1(0);
	R(1,0) = b1(1);
	R(2,0) = b1(2);

	R(0,1) = b2(0);
	R(1,1) = b2(1);
	R(2,1) = b2(2);

	R(0,2) = b3(0);
	R(1,2) = b3(1);
	R(2,2) = b3(2);


	IntermediateState skew_error = 0.5*(R_d.transpose()*R - R.transpose()*R_d);
	IntermediateState e_R("",3,1);
	e_R(0) = skew_error(2,1);
	e_R(1) = skew_error(0,2);
	e_R(2) = skew_error(1,0);
	IntermediateState e_W = omega - R.transpose()*R_d*omega_d;
	IntermediateState e_T = (u(0)+u(1)+u(2)+u(3)) - T_d_norm;
	IntermediateState e_U("",4,1);
	e_U(0) = u(0) - (T_d_norm/4);
	e_U(1) = u(1) - (T_d_norm/4);
	e_U(2) = u(2) - (T_d_norm/4);
	e_U(3) = u(3) - (T_d_norm/4);
	
	//Defining cost function h
	Function h;
	h << e_R;
	h << e_W;
	h << e_T;
	h << e_U;
	h << sv;

	//Defining cost matrix Q
	DMatrix Q(h.getDim(), h.getDim());
  	Q.setIdentity();
  	//Q(6) = 0.001;
	
	//Defining end cost function hN
	//hN cannot depend on input u/sv
  	Function hN;
  	hN << e_R;
	hN << e_W;
	//hN << e_T;
	//hN << e_U;
	//hN << sv;

	//Defining cost matrix QN
	DMatrix QN(hN.getDim(), hN.getDim());
  	QN.setIdentity();

  	//Defining optimal control problem
  	OCP ocp(t_start, t_end, N);
  	ocp.minimizeLSQ(Q, h);
    ocp.minimizeLSQEndTerm(QN, hN);	

  	//Add system constraints
  	ocp.subjectTo(f); //System dynamics

  	//	Following type of constraints implemented:
  	ocp.subjectTo(f_min <= u(0) <= f_max);
  	ocp.subjectTo(f_min <= u(1) <= f_max);
  	ocp.subjectTo(f_min <= u(2) <= f_max);
	ocp.subjectTo(f_min <= u(3) <= f_max);

	/*ocp.subjectTo(omega_min <= omega(0) + sv(0));
	ocp.subjectTo(omega(0)-sv(0) <= omega_max);

	ocp.subjectTo(omega_min <= omega(1) + sv(1));
	ocp.subjectTo(omega(1)-sv(1) <= omega_max);

	ocp.subjectTo(omega_min <= omega(2) + sv(2));
	ocp.subjectTo(omega(2)-sv(2) <= omega_max);

	ocp.subjectTo(sv >= 0);
	*/

	ocp.subjectTo(omega_min <= omega(0) <= omega_max);
	ocp.subjectTo(omega_min <= omega(1) <= omega_max);
	ocp.subjectTo(omega_min <= omega(2) <= omega_max);
	//ocp.subjectTo(AT_START, sv == 1);

	ocp.subjectTo(cos(theta) <= t_end*b3(2)-sv);

	//Number of online data
	ocp.setNOD(NR_ONLINE_DATA);

	OCPexport mpc(ocp);

    mpc.set(INTEGRATOR_TYPE,        INT_IRK_GL4);
    mpc.set(HESSIAN_APPROXIMATION,  GAUSS_NEWTON);
    mpc.set(DISCRETIZATION_TYPE,    MULTIPLE_SHOOTING);
    mpc.set(SPARSE_QP_SOLUTION,     FULL_CONDENSING_N2);
    mpc.set(NUM_INTEGRATOR_STEPS,   N);
    mpc.set(QP_SOLVER,              QP_QPOASES);
    mpc.set(HOTSTART_QP,            YES);
    mpc.set( USE_SINGLE_PRECISION,        YES);

    mpc.set(GENERATE_TEST_FILE, 	YES);
    mpc.set(GENERATE_MAKE_FILE, 	YES);
    mpc.set( GENERATE_MATLAB_INTERFACE,   NO);
    mpc.set( GENERATE_SIMULINK_INTERFACE, NO);
    mpc.set(INTEGRATOR_TOLERANCE, 1e-6 );
    mpc.set( KKT_TOLERANCE, 1e-4 );



    if(mpc.exportCode("nonlinear_mpc_control_export") != SUCCESSFUL_RETURN){
    	exit( EXIT_FAILURE );
    	mpc.printDimensionsQP( );
  	}
	return EXIT_SUCCESS;
}