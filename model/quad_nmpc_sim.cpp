#include <acado_optimal_control.hpp>
#include <acado_toolkit.hpp>
#include <acado_gnuplot.hpp>
#include <matrix_vector/matrix.hpp>
#include <iostream>

using namespace std;
USING_NAMESPACE_ACADO

int main(){

	//variables
	const double t_start = 0.0;
  	const double t_end = 60;//50.0;
  	const double dt = 0.25;//0.25;
  	const int 	 N = round(t_end/dt);
  	const double f_min = 0;
  	const double f_max = 4;
  	const double omega_min = -20; //dont really care
  	const double omega_max = omega_min*-1;
  	const double theta = 3.14/9;
  	const double phi_init = 3.14/3;
  	const double theta_init = 3.839;
  	const double psi_init = 3.14/3;
  	const double psi_d = 3.14/3;
  	const double c = 0.03793; //prob wrong
  	const double mass = 0.468;
  	const double g = 9.81;


	//NMPC Differential states
	DifferentialState 	b1("",3,1);
	DifferentialState 	b2("",3,1);
	DifferentialState 	b3("",3,1);
	DifferentialState 	omega("",3,1);

	//Dummy state to introduce slack var
	DifferentialState dummy;

	//Online data, "references"
	/*OnlineData 	R_d("",3,3);
	OnlineData 	omega_d("",3,1);
	OnlineData 	T_d("",3,1);
	OnlineData 	t;
	const int NR_ONLINE_DATA = 4;
	*/
	//Control signals
	Control u("",4,1);

	//Slack variable
	Control sv;

	const double I_xx 	= 0.0049; //Roll inertia
	const double I_yy 	= 0.0049; //Pitch inertia
	const double I_zz 	= 0.0088; //Yaw inertia
	const double l 		= 0.225;

	IntermediateState T_phi 	= 	l*(-u(1)+u(3));
	IntermediateState T_theta 	= 	l*(-u(0)+u(2));
	IntermediateState T_psi 	= 	l*(u(0)-u(1)+u(2)-u(3));
	
	IntermediateState T_d("",3,1);
	T_d(0) = 0;
	T_d(1) = 0;
	T_d(2) = g;
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

	IntermediateState omega_d("",3,1);
	omega_d(0) = 0;
	omega_d(1) = 0;
	omega_d(2) = 0;

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

	f	<<	dot(omega(0)) == -1*omega(1)*omega(2)*(I_zz-I_yy)/(I_xx) + T_phi/I_xx;
	f	<<	dot(omega(1)) == -1*omega(0)*omega(2)*(I_xx-I_zz)/(I_yy) + T_theta/I_yy;
	f	<<	dot(omega(2)) == -1*omega(0)*omega(1)*(I_yy-I_xx)/(I_zz) + T_psi/I_zz;

	//mock dynamics of slack variable
	f << dot(dummy) == sv;

	//Setup simulated process
	OutputFcn identity;
	DynamicSystem dynamicSystem(f, identity);
	Process process(dynamicSystem, INT_RK78);

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
	IntermediateState e_omega = omega - R.transpose()*R_d*omega_d;
	IntermediateState e_T = (u(0)+u(1)+u(2)+u(3)) - T_d_norm;
	IntermediateState e_U("",4,1);
	e_U(0) = u(0) - (T_d_norm/4);
	e_U(1) = u(1) - (T_d_norm/4);
	e_U(2) = u(2) - (T_d_norm/4);
	e_U(3) = u(3) - (T_d_norm/4);
	
	//Defining cost (measurement) function h
	Function h;
	h << e_R;
	h << e_omega;
	h << e_T;
	h << e_U;
	h << sv;

	//Defining cost matrix Q
	DMatrix Q(h.getDim(), h.getDim());
  	Q.setIdentity();
  	Q(13) = 0.01;
	
  	// Set a reference for the analysis (if CODE_GEN is false).
  	DVector r(h.getDim());
  	r.setAll(0);
	
	//Defining end cost function hN
  	Function hN;
  	hN << e_R;
	hN << e_omega;
	//hN << e_T;
	//hN << e_U;
	//hN << sv;

	//Defining cost matrix QN
	DMatrix QN(hN.getDim(), hN.getDim());
  	QN.setIdentity();

  	// Set a reference for the analysis (if CODE_GEN is false).
  	DVector rN(h.getDim());
  	rN.setAll(0);

  	//Defining optimal control problem
  	OCP ocp(t_start, t_end, N);

  	ocp.minimizeLSQ(Q, h, r);
    ocp.minimizeLSQEndTerm(QN, hN, rN);	

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
	*/
	ocp.subjectTo(sv >= 0);
	

	ocp.subjectTo(omega_min <= omega(0) <= omega_max);
	ocp.subjectTo(omega_min <= omega(1) <= omega_max);
	ocp.subjectTo(omega_min <= omega(2) <= omega_max);

	ocp.subjectTo(cos(theta) <= b3(2)+sv);

	//Number of online data
	//ocp.setNOD(NR_ONLINE_DATA);

    RealTimeAlgorithm algorithm(ocp, dt);
    algorithm.set( INTEGRATOR_TYPE, INT_RK78 );
    algorithm.set( INTEGRATOR_TOLERANCE, 1e-6 );
    algorithm.set( KKT_TOLERANCE, 1e-4 );

    StaticReferenceTrajectory zeroReference;

    Controller controller(algorithm, zeroReference);

    //sim params
    double simStartTime = 0.0;
    double simEndTime = 8;

    SimulationEnvironment sim(simStartTime,simEndTime,process,controller);

    DVector x0(13);
    x0.setZero( );
    x0(0)=cos(psi_init)*cos(theta_init);
    x0(1)=sin(psi_init)*cos(theta_init);
    x0(2)=-sin(theta_init);
    x0(3)=cos(psi_init)*sin(theta_init)*sin(phi_init)-sin(psi_init)*cos(phi_init);
    x0(4)=sin(psi_init)*sin(theta_init)*sin(phi_init)+cos(psi_init)*cos(phi_init);
    x0(5)=cos(theta_init)*sin(phi_init);
    x0(6)=cos(psi_init)*sin(theta_init)*cos(phi_init)+sin(psi_init)*sin(phi_init);
    x0(7)=sin(psi_init)*sin(theta_init)*cos(phi_init)-cos(psi_init)*sin(phi_init);
    x0(8)=cos(theta_init)*cos(phi_init);
    x0(9)=0;
    x0(10)=0;
    x0(11)=0;
    x0(12)=1;

    if (sim.init( x0 ) != SUCCESSFUL_RETURN)
    	exit( EXIT_FAILURE );
    if (sim.run( ) != SUCCESSFUL_RETURN)
    	exit( EXIT_FAILURE );

    VariablesGrid diffStates;
    sim.getProcessDifferentialStates(diffStates);

    VariablesGrid feedbackControl;
    sim.getFeedbackControl(feedbackControl);

    GnuplotWindow window;
    window.addSubplot(diffStates(0),"b1(0)" );
    window.addSubplot(diffStates(3),"b2(0)" );
    window.addSubplot(diffStates(6),"b3(0)" );

    window.addSubplot(diffStates(1),"b1(1)" );
    window.addSubplot(diffStates(4),"b2(1)" );
	window.addSubplot(diffStates(7),"b3(1)" );

    window.addSubplot(diffStates(2),"b1(2)" );
	window.addSubplot(diffStates(5),"b2(2)" );
	window.addSubplot(diffStates(8),"b3(2)" );
		
    GnuplotWindow window2;
	window2.addSubplot(feedbackControl(0),"u(0)" );
	window2.addSubplot(feedbackControl(1),"u(1)" );
	window2.addSubplot(feedbackControl(2),"u(2)" );
	window2.addSubplot(feedbackControl(3),"u(3)"  );

	GnuplotWindow window3;
	window3.addSubplot(diffStates(9),"omega(0)" );
	window3.addSubplot(diffStates(10),"omega(1)" );
	window3.addSubplot(diffStates(11),"omega(2)"  );

    window.plot();
    window2.plot();
    window3.plot();

	return EXIT_SUCCESS;
}