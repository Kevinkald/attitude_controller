#include "ros/ros.h"
//#include "geometry_msgs/Point.h"
//#include <mavros_msgs/ActuatorControl.h>

#include "acado_common.h"
#include "acado_auxiliary_functions.h"

#include <stdio.h>

#include <iostream>

/* Some convenient definitions. */
#define NX          ACADO_NX  /* Number of differential state variables.  */
#define NXA         ACADO_NXA /* Number of algebraic variables. */
#define NU          ACADO_NU  /* Number of control inputs. */
#define NOD         ACADO_NOD  /* Number of online data values. */

#define NY          ACADO_NY  /* Number of measurements/references on nodes 0..N - 1. */
#define NYN         ACADO_NYN /* Number of measurements/references on node N. */

#define N           ACADO_N   /* Number of intervals in the horizon. */

#define NUM_STEPS   10        /* Number of real-time iterations. */
#define VERBOSE     1         /* Show iterations: 1, silent: 0.  */

/* Global variables used by the solver. */
ACADOvariables acadoVariables;
ACADOworkspace acadoWorkspace;

/* A template for testing of the solver. */
int main(int argc, char **argv )
{	
	// Reset all solver memory
	memset(&acadoWorkspace, 0, sizeof( acadoWorkspace ));
	memset(&acadoVariables, 0, sizeof( acadoVariables ));

	ros::init(argc, argv, "attitude_controller");

	ros::NodeHandle nh;

	ros::Rate loop_rate(50);
	/* Some temporary variables. */
	int    i, iter;
	acado_timer t;
	real_t t1=0, t2=0; 

	/* Initialize the solver. */
	acado_initializeSolver();

	/* Initialize the states and controls. */
	for (i = 0; i < N+1; ++i) {
		acadoVariables.x[0 + i*NX] = 1.0; 
		acadoVariables.x[1 + i*NX] = 0.0;
		acadoVariables.x[2 + i*NX] = 0.0;
		acadoVariables.x[3 + i*NX] = 0.0;
		acadoVariables.x[4 + i*NX] = 1.0;
		acadoVariables.x[5 + i*NX] = 0.0;
		acadoVariables.x[6 + i*NX] = 0.0;
		acadoVariables.x[7 + i*NX] = 0.0;
		acadoVariables.x[8 + i*NX] = 1.0;
		acadoVariables.x[9 + i*NX] = 0;
		acadoVariables.x[10 + i*NX] = 0;
		acadoVariables.x[11 + i*NX] = 0;
		acadoVariables.x[12 + i*NX] = 1.0;
	}

	
	for (i = 0; i < NU * N; ++i)  acadoVariables.u[ i ] = 0.0;




	/* Initialize the measurements/reference. */
	for (i = 0; i < NY ; ++i){
		acadoVariables.y[0 + i*NY] = 1.0;

	}
	for (i = 0; i < NYN; ++i)  acadoVariables.yN[ i ] = 0.0;





	/* MPC: initialize the current state feedback. */

	for (i = 0; i < NX; ++i) acadoVariables.x0[ i ] = acadoVariables.x[ i ];

	
	for (i = 0; i < NOD* (N + 1); ++i)  acadoVariables.od[ i ] = 0.0;

	if( VERBOSE ) acado_printHeader();

	acado_preparationStep();
	

	//acado_printControlVariables();
	//acado_printDifferentialVariables();

	//while(ros::ok()) {
	for(iter = 0; iter < NUM_STEPS; iter++ ) {
		acado_tic( &t );
		int status = acado_feedbackStep( );
		t2 = acado_toc( &t );
		if (status) {
			std::cout << "Iteration:" << iter << ", QP problem! QP status: " << status << std::endl;
			break;
		}

		acado_printDifferentialVariables();
		acado_printControlVariables();


		//
		// Prepare for the next iteration
		//

		// In this simple example, we feed the NMPC with an ideal feedback signal
		// i.e. what NMPC really expects in the next sampling interval
		for (int i = 0; i < NX; ++i)
			acadoVariables.x0[ i ] = acadoVariables.x[NX + i];

		// Shift states and control and prepare for the next iteration
		acado_shiftStates(2, 0, 0);
		acado_shiftControls( 0 );

		acado_tic( &t );
		acado_preparationStep();
		t1 = acado_toc( &t );
//
		//for (i = 0; i < NY * N; ++i)  acadoVariables.y[ i ] = 0.0;
		//for (i = 0; i < NYN; ++i)  acadoVariables.yN[ i ] = 0.0;
		//acadoVariables.x0[0] = 1.0; 
		//acadoVariables.x0[1] = 0.0;
		//acadoVariables.x0[2] = 0.0;
		//acadoVariables.x0[3] = 0.0;
		//acadoVariables.x0[4] = 1.0;
		//acadoVariables.x0[5] = 0.0;
		//acadoVariables.x0[6] = 0.0;
		//acadoVariables.x0[7] = 0.0;
		//acadoVariables.x0[8] = 1.0;
		//acadoVariables.x0[9] = 0.1;
		//acadoVariables.x0[10] = 0.1;
		//acadoVariables.x0[11] = 0.1;
//		//acadoVariables.x0[12] = 1.0;
//
//		//for (i = 0; i < NX; ++i)
//		//	acadoVariables.x[ i ] = acadoVariables.x0[ i ];
//
//		//std::cout << acadoVariables.x0[0] << std::endl;
//		//std::cout << acadoVariables.x[0] << std::endl;
//		///* Initialize the online data. */
//		///*
//		//	T_d is received from LQR, omega_d, psi_d and t is set to zero
//		//*/
//		//for (i = 0; i < N ; ++i){
//		//	acadoVariables.od[i*(NOD)+0] = 0.1;//T_d.x;
//		//	acadoVariables.od[i*(NOD)+1] = 0.1;//T_d.y;
//		//	acadoVariables.od[i*(NOD)+2] = 0.1;//T_d.z;
//		//	acadoVariables.od[i*(NOD)+3] = 0.1;//T_d.x;
//		//	acadoVariables.od[i*(NOD)+4] = 0.1;//T_d.x;
//		//	acadoVariables.od[i*(NOD)+5] = 0.1;//T_d.y;
//		//	acadoVariables.od[i*(NOD)+6] = 0.1;//T_d.z;
//		//	acadoVariables.od[i*(NOD)+7] = 0.1;//T_d.z;
//		//}
//		//acado_preparationStep();
//		///* Get the time before start of the loop. */
//		//acado_tic( &t );
//
//		///* The "real-time iterations" loop. */
//		////for(iter = 0; iter < NUM_STEPS; ++iter)
//		////{
//        ///* Perform the feedback step. */
//		// acado_feedbackStep( );
//
//		// real_t te = acado_toc( &t );
		

		//int status = acado_preparationStep();
		//std::cout << "Iteration:" << iter << ", QP problem! QP status: " << status << std::endl;

		/* Apply the new control immediately to the process, first NU components. */

		if( VERBOSE ) printf("\tReal-Time Iteration %d:  KKT Tolerance = %.3e\n\n", iter, acado_getKKT() );

		/* Optional: shift the initialization (look at acado_common.h). */
        /*acado_shiftStates(2, 0, 0); 
		acado_shiftControls( 0 ); 
		*/
		/* Prepare for the next step. */

		//}
		/* Read the elapsed time. */
		

		//if( VERBOSE ) printf("\n\nEnd of the RTI loop. \n\n\n");

		/* Eye-candy. */

		//if( !VERBOSE )
		//printf("\n\n Average time of one real-time iteration:   %.3g microseconds\n\n", 1e6 * te / NUM_STEPS);

		//ros::spinOnce();
		//loop_rate.sleep();
	}
    return 0;
}
