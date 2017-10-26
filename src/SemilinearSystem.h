/**
*	@file	SemilinearSystem
*
*	@class 	SemilienarSystem
*
*	@brief	The class allows for simualting transport of natural gas through a single pipe utilizing a semilienar system and for a user-provided friction coefficient and friction function.
*
*	The class 'SemilinearSystem' simulates the transport of natural gas through a single pipe utilizing a semilinear model that is given as <br>
*	<center>
*				p_t + q_x = 0 <br>
*				q_t + c_{SOS}p_x	=	Lamdba a(p,q) <br>
*	</center>
*	on the domain Omega = (0,T)x(x_L, x_R) with initial conditions <br>
*	<center>
*		p(0,x) = p_0(x) <br>
*		q(0,x) = q_0(x)	<br>
*	</center>
*	and boundary conditions <br>
*	<center>
*		q(t,x_L)	=	q_L(t)	<br>
*		q(t,x_R)	=	q_R(t)	<br>
*	</center>
*	The function a(p,q) is refferred to as friction function while Lambda denotes the spatially disctributed friction coefficient.
*	
*	@author	Dr. N. Strogies 	strogies@wias-berlin.de
*
*	@date 07.09.2017
*
*	@version .5
*/


#ifndef SEMILINEARSYSTEM_H
#define SEMILINEARSYSTEM_H

#include <iostream>
#include <fstream>
#include <cmath>
#include <string.h>


//	Including the functions specifying boundary data for Q at left and right boundary
#include "Boundary_Data_For_Q_L.h"
#include "Boundary_Data_For_Q_R.h"

//	Including the functions specifying intial data for P and Q
#include "Initial_Data_For_P.h"
#include "Initial_Data_For_Q.h"

//	Including the functions specifying the functions utilized in the expansion of the friction coefficient
#include "Value_Lambda_Base.h"

//	Including the function specifying the frictino function

#include "SourceTerm.h"


// #include <conio.h>


using namespace std;

class SemilinearSystem
{
 public:
  //	Declaration of functions
  SemilinearSystem(double D_SpeedOfSound, double D_TerminalTime,
		   double D_Boundary_Position_Left, double D_Boundary_Position_Right,
		   double D_Delta_x, int I_Lambda_Expansion_Length, double eps); // Constructor
		
  // giving info about pipe in the terminal
  void info();

  //	PUBLIC FUNCTIONS
  //	RETURN-TYPE void

  void EvalTest(  ); //	Testing function that is for development only and will be erased later

  /*
    Function that can be callled from main code that provides the
    coefficients of LAMBDA-Expansion and solves the main problem for
    these values
   */
  void Run( double DA_P_Lambda_Coefficients[], bool write2file_bool = true ); 

  /*
    Function that can be callled from main code that provides the
    coefficients of LAMBDA-Expansion and solves the main problem for these
    values
   */
  void Run_Text_Output( double DA_P_Lambda_Coefficients[] ); 
  
  void Write2File(char* filename, bool append = true);
  
  // 	RETURN-TYPE double
  //	Observation Operator
  double get_P_Left_Boundary( double D_EvaluationTime_Arg); //	Returns value of P on the left boundary for given time
  double get_P_Right_Boundary( double D_EvaluationTime_Arg); //	Returns value of P on the right boundary for given time
  double get_P_Difference( double D_EvaluationTime_Arg); //	Returns difference of P between boundaries.

  int CurrentTimeIndex();
  double* BoundaryValueP_Left();
  double* BoundaryValueP_Right();
  double* TimeSlices();
  double* LambdaAverage(double DA_P_Lambda_Coefficients_GIVEN[]);
  double* SendLambdaAverage();
  int NumberofCells();
 private:
  //	Declaration of VARIABLES
  //	TYPE: double - single variable
  ///	Width of the spatial Discretization

  /**	
	Double that stores the width of the spatial discretization, i.e., the
	width of the computational cells. It is connected to D_Delta_t and
	D_SpeedOfSound by D_Delta_t D_SpeedOfSound = D_Delta_x 
  */
  double	D_Delta_x;						//	size of spatial discretization
  ///	Speed of Soud in the underlying physical system.
  /**	Double that stores the speed of sound in the underlying physical system. This value determines the propagation-speed of information through the system.	*/
  double	D_SpeedOfSound;					//	speed of sound

  ///	Time Step width
  /**	Double that stores the time step width utilized in the
	discretization. It is connected to D_Delta_x and
	D_SpeedOfSound by D_Delta_t D_SpeedOfSound = D_Delta_x 
  */
  double 	D_Delta_t;						//	time step size
  ///	Time horizon of the dynamic problem
  /**	Double that stores the runtime T of the underlying physical
	system. At best, it is some multiple of D_Delta_t. 
  */
  double	D_TerminalTime;					//	Terminal Time
  ///	Coordinate of the left boundary of the pipe
  /**	Double that stores the coordinate of the left boundary of the
	pipe. Usually, it is set to zero. 
  */
  double	D_Boundary_Position_Left;		//	position of left boudnary
  ///	Coordinate of the right boundary of the pipe
  /**	Double that stores the right boundary of the pipe. Usually, it
	is set to the length of the pipe.  
  */
  double	D_Boundary_Position_Right;		//	position of right boundary
  
  //	TYPE: int
  ///	Number of summands in the Lambda-expansion
  /**	Integer that stores the number of summands in the expansion of
	Lambda after which this representation is truncated.  
  */
  int		I_Lambda_Expansion_Length;		//	Number of coefficients for the expansion that are considered (trunaction of the sum after this number)
  ///	Number of time steps
  /**	Integer that stores the number if time steps in the numerical
	scheme.  
  */
  int 	I_NumberOfTimeSteps;			//	Number of Time Steps for the forwad computation
  ///	Number of interior cells in the pipe
  /**	Integer that stores the number of interior cells in the
	numerical scheme. 
  */
  int		I_NumberOfCells;				//	Number of Cells in the discretization
  
  //	TYPE: double - arrays with dynamic memory allocation
  ///	Values of P in a single time step
  /**	Pointer to Double Array that contains the integral averages of
   *	P and values of P on the ghost cells in the active time
   *	slice. Since the boundary values of P are stored in a
   *	different array and no information on the actual evolution of
   *	P have to be stored, the values of this array are overwritten
   *	in every time step. 
   */
  double	* DA_P_Values_P;				//	Pointer on Internal p-Array
  ///	Values of Q in a single time step
  /**	Pointer to Double Array that contains the integral averages of
   *	Q and values of Q on the ghost cells in the active time
   *	slice. Since no information on the actual evolution of Q have
   *	to be stored, the values of this array are overwritten in
   *	every time step.  
   */
  double	* DA_P_Values_Q;				//	Pointer on Internal q-Array
  ///	Values of P in the left ghost cells for all time steps.
  /**	Pointer to Double Array that stores the values of P on the
	ghost cells left of the pipe, thus approximating the trace
	aveluation of P on this boundary.
  */
  double	* DA_P_Left_P;					//	Pointer on boundary values of p at left boundary
  ///	Values of P in the right ghost cells for all
  ///	time steps.
  /**	Pointer to Double Array that stores the values of P on the
	ghost cells right of the pipe, thus approximating the trace
	aveluation of P on this boundary.  */
  double 	* DA_P_Right_P;					//	Pointer on boundary values of p at right boundary
  ///	Integral averages of Lambda in all of the cells including the
  ///	ghost cells on the boundaries.
  /**	Pointer to Double Array that stores the integral averages of
   *	Lambda in every interior cell and the linear extrapolation
   *	depending on the two left- and right-most interior values on
   *	the ghost cells. The elements on this array are re-defined in
   *	every .Run call of the object since the runs only difer in the
   *	choice of Lambda.  */
  double 	* DA_P_Lambda_AV;				//	Pointer on Lambda-Averages on every cell of the discretization
  ///	Coefficients of the Lambda expansion defining the Lambda in a
  ///	certain .Run call.
  /**	Pointer to Double Array that stores the expectational value of
	Lambda on the first component and the coefficients of the
	Lambda-expansion on the following components.  */
  double * DA_P_Lambda_Coefficients;		//	Coefficients of the expansion of the random fiel (will be defined as a double-array with user-provided number of coefficients)
  
  double *DA_Centroid;
  
  /* current time */
  double D_Current_T;
  /* current time index */
  int I_Current_T_Index;
  /* epsilon for integration at the boundary */
  double epsilon_boundary;
  /* number of cells that fall into epsilon neighborhood */
  double n_epsilon;
  /* time slices */
  double * time_slices;
  /* length of the domain */
  double domain_length;
  
  //	Declaration of FUNCTIONS
		//	RETURN-VALUE void
  void Set_IC_Q( double* DA_P_Values_Q, int I_NumberOfCells, double D_Boundary_Position_Left, double D_Delta_x );	//	Stores the integral averages of the initial data for q on the array 'DA_P_Values_Q'
  void Set_IC_P( double* DA_P_Values_P, int I_NumberOfCells, double D_Boundary_Position_Left, double D_Delta_x );	//	Stores the integral averages of the initial data for p on the array 'DA_P_Values_P'
  void Set_Lambda_Averages( double*, double*, int, int, double, double  );
  
  //	RETURN-VALUE double		
  double 	Read_BoundaryValues_Q_Left(double D_EvaluationTime);		//	HANDSHAKE: evaluation of the boundary value of q at the 'left' boundary at time t
  double	Read_BoundaryValues_Q_Right(double D_EvaluationTime);		// 	HANDSHAKE: evaluation of the boundary value of q at the 'right' boundary at time t

};


#endif



