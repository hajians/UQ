#include "SemilinearSystem.h"

/**
*	@file	SemilinearSystem
*/

SemilinearSystem::SemilinearSystem(double D_SpeedOfSound_Arg,
				   double D_TerminalTime_Arg,
				   double D_Boundary_Position_Left_Arg,
				   double D_Boundary_Position_Right_Arg,
				   double D_Delta_x_Arg,
				   int I_Lambda_Expansion_Length_Arg,
				   double eps)
{
	/**
	*	@brief		Constructor of the class 'SemilinearSystem'
	*
	*	Any instance of 'SemilinearSystem' represents a
	*	certain setting of the gas-transportation model
	*	through a single pipe with repect to the underlying
	*	physics (such as 'speed of sound', 'length of the
	*	pipe', 'terminal time') as well as the discretization
	*	(such as 'mesh width', 'expansion size of Lambda')
	*	utilized for solving the problem. The only missing
	*	part, the coefficients in the expansion of Lambda, is
	*	user-provided in every run of the problem.\n The basic
	*	idea is, that by constructing an instance of this
	*	class, the physical and computational setting is fixed
	*	and only evaluated for several choices of Lambda.\n
	*
	*	@param[in] D_SpeedOfSound Speed of sound for the
	*	problem (defines the speed of propagation for any
	*	information in the semilinear setting)
	*
	*	@param[in] D_TerminalTime Time horizion, the dynamic
	*	problem has to be solved in.
	*
	*	@param[in] D_Boundary_Position_Left Coordinate of the
	*	left boundary of the pipe (provides length of the pipe
	*	with D_Boundary_Position_Right)
	*
	*	@param[in] D_Boundary_Position_Right Coordinate of the
	*	right boundary of the pipe (provides length of the
	*	pipe with D_Boundary_Position_Left)
	*
	*	@param[in] D_Delta_x Mesh width of the spatial
	*	discretiazetion @param[in] I_Lambda_Expansion_Length
	*	Number of coefficients in the expansion of Lambda
	*
	*	@version 1.0
	*
	*	@date 07.09.2017		
	*/
  
  //	Declaration of variables TYPE DOUBLE
  
  D_Delta_x = D_Delta_x_Arg;  //	Width of spatial discretization
  
  D_SpeedOfSound = D_SpeedOfSound_Arg; //	Speed of Sound
  D_Delta_t 					=	D_Delta_x/D_SpeedOfSound; //	time step size
  D_TerminalTime 				= 	D_TerminalTime_Arg; //	Terminal Time
  D_Boundary_Position_Left	=	D_Boundary_Position_Left_Arg; //	coordinate of the 'left' boundary
  D_Boundary_Position_Right 	=	D_Boundary_Position_Right_Arg; //	coordinate of the 'right' boundary
  
  
  //	Declaration of variables TYPE INT
  I_Lambda_Expansion_Length = I_Lambda_Expansion_Length_Arg; //	Number of coefficients that are included in the truncated representation of the friction coefficent Lambda
  I_NumberOfTimeSteps = D_TerminalTime/D_Delta_t;												//	#TimeSteps
  I_NumberOfCells = abs(D_Boundary_Position_Right - D_Boundary_Position_Left)/D_Delta_x;	//	#Cells

  //	DATA-Stores
  //	(the public .run - method will operate on DA_P_Values_P and DA_P_Values_Q)
  DA_P_Values_P = new double [I_NumberOfCells + 2];
  DA_P_Values_Q = new double [I_NumberOfCells + 2];
	
  //	(the public methods that represent the observation operator will act on the elements of these vectors)
  DA_P_Left_P  = new double [I_NumberOfTimeSteps + 1];
  DA_P_Right_P = new double [I_NumberOfTimeSteps + 1];
	
  //	LAMBDA-EXPANSION:
  //	The following vector will contain the integral average of LAMBDA (the truncated representation)
  DA_P_Lambda_AV = new double [I_NumberOfCells + 2];
  //	Initialization of the array storing the coefficients of the expansion of the radmon field
  DA_P_Lambda_Coefficients = new double [I_Lambda_Expansion_Length];

  DA_Centroid = new double[I_NumberOfCells];
  
  DA_Centroid[0] = D_Boundary_Position_Left_Arg + 0.5 * D_Delta_x;
  
  for (int i=1; i < I_NumberOfCells; i++)
    {
      DA_Centroid[i] = DA_Centroid[i-1] + D_Delta_x;
    }
  
  // current time
  D_Current_T = 0.0;

  // current time index
  I_Current_T_Index = 0;

  // epsilon for the boundary
  epsilon_boundary = eps;
  // number of cells that fall into epsilon neighborhood
  n_epsilon = (int) (eps / D_Delta_x);

  //
  time_slices = new double [I_NumberOfTimeSteps + 1];
  time_slices[0] = 0.0;
    
  for (int i=1; i<=I_NumberOfTimeSteps; i++){
    time_slices[i] = time_slices[i-1] + D_Delta_t;
  }

  /* domain length */
  domain_length = D_Boundary_Position_Right_Arg - D_Boundary_Position_Left_Arg;

}



//	Definition of functions that provide important data for the computation
//	NOTE: These functions are supposed to be defined by the user. We have impemented a 
//	1)	BOUNDARY DATA
//	The following functions provide boundary data for Q at the left and right boundary of the pipeline

double SemilinearSystem::Read_BoundaryValues_Q_Left( double D_EvaluationTime_Arg )				//	Evaluation Tested!
{
  /**
   *	@brief Reads the values of the user-provided function that
   *	specifies the boundary values of q at the left boundary.
   *
   *	The function calles an user-provided function that specifies
   *	the boundary data of Q at the left boundary and returns the
   *	result.
   *
   *	@param[in] D_EvaluationTime Time for which the right boundary
   *	value for Q should be evaluated
   *
   *	@param[out] D_ValueAtT Value of the boundary data of Q at the
   *	right boundary at time D_EvaluationTime
   *
   *	@date 07.09.2017
   */
	
	// Declaration of Variable
	double D_ValueAtT;
	
	// Evaluation of Boundary Data
	// At this point, the function for the boundary data at the left boundary is necessary
	D_ValueAtT 					=	Value_Q_L( D_EvaluationTime_Arg );
	
	//	Return the computed value
	return D_ValueAtT;
}

double SemilinearSystem::Read_BoundaryValues_Q_Right( double D_EvaluationTime_Arg )				//	Evaluation Tested
{
	/**
	*	@brief Reads the values of the user-provided function
	*	that specifies the boundary values of q at the right
	*	boundary.
	*
	*	The function calles an user-provided function that
	*	specifies the boundary data of Q at the right boundary
	*	and returns the result.
	*
	*	@param[in] D_EvaluationTime Time for which the right
	*	boundary value for Q should be evaluated

	*	@param[out] D_ValueAtT Value of the boundary data of Q
	*	at the right boundary at time D_EvaluationTime
	*
	*	@date 07.09.2017
	*/
	
	// Declaration of Variable
	double D_ValueAtT;
	
	// Evaluation of Boundary Data
	// At this point, the function for the boundary data at the left boundary is necessary
	D_ValueAtT 					=	Value_Q_R( D_EvaluationTime_Arg );
	
	//	Return the computed value
	return D_ValueAtT;
 } 

// 	2) INITIAL DATA
//	The following functions provide the initial data of P and Q, allready in integral-average
void SemilinearSystem::Set_IC_Q( double* DA_P_Values_Q_Arg, int I_NumberOfCells_Arg,
				 double D_Boundary_Position_Left_Arg,
				 double D_Delta_x_Arg )	//	Evaluation Tested
{
	/**
	*	@brief	Initializes the integral averages of Q based on the initial condition.
	*
	*	The function intializes the array containing the
	*	integral averages of Q_0 on the cosen
	*	discretization.\n For the internal cells, this is
	*	obtained by evaluating an user-provided function that
	*	specifies the initial data for Q.  The values for the
	*	ghost cells stem from the boundary data for Q.
	*
	*	@param[in] DA_P_Values_Q Memory location for the
	*	integral averages of Q for the internal and ghost
	*	cells frespectively. (size: Number of Cells + 2)
	*
	*	@param[in] I_NumberOfCells Number of cells in the
	*	discretization
	*
	*	@param[in] D_Boundary_Position_Left Coordinate of the
	*	left boundary of the pipe
	*
	*	@param[in] D_Delta_x Width of the spatial
	*	discretization
	*
	*	@date 07.09.2017
	*
	*	@todo	better integration scheme
	*/
	
	//	Fixing the values at the boundaries
	DA_P_Values_Q[0]					=	SemilinearSystem::Read_BoundaryValues_Q_Left(0.0);
	DA_P_Values_Q[I_NumberOfCells + 1]	=	SemilinearSystem::Read_BoundaryValues_Q_Right(0.0);
	
	for (int i = 1; i <= I_NumberOfCells; ++i)
	{
		DA_P_Values_Q[i] 				=	( Value_Q_0( D_Boundary_Position_Left + ( i-1 )*D_Delta_x ) + Value_Q_0( D_Boundary_Position_Left + ( i )*D_Delta_x ) )/2.0;
	}
}

void SemilinearSystem::Set_IC_P( double* DA_P_Values_P, int I_NumberOfCells, double D_Boundary_Position_Left, double D_Delta_x )
{
  /**
   *	@brief	Initializes the integral averages of P based on the initial condition.
   *
   *	The function intializes the array containing the integral
   *	averages of P_0 on the cosen discretization.\n For the
   *	internal cells, this is obtained by evaluating an
   *	user-provided function that specifies the initial data for P.
   *	The values for the ghost cells are obtained by a linear
   *	function based on the integral averages of P on the 2 most
   *	left and 2 most right internal cells.
   *
   *	@param[in] DA_P_Values_P Memory location for the
   *	integral averages of P for the internal and ghost
   *	cells frespectively. (size: Number of Cells + 2)
   *	@param[in] I_NumberOfCells Number of cells in the
   *	discretization @param[in] D_Boundary_Position_Left
   *	Coordinate of the left boundary of the pipe @param[in]
   *	D_Delta_x Width of the spatial discretization
   *
   *	@date 07.09.2017
   *
   *	@todo	better integration scheme
   */
  
  //	Fixing the values at the boundaries (Have to be given and to
  //	be consistent with the spatial averages)
  DA_P_Values_P[0]			= ( Value_P_0( D_Boundary_Position_Left - D_Delta_x ) + Value_P_0( D_Boundary_Position_Left ) )/2.0;
  DA_P_Values_P[I_NumberOfCells + 1]	= ( Value_P_0( D_Boundary_Position_Left + (I_NumberOfCells)*D_Delta_x ) 
					    + Value_P_0( D_Boundary_Position_Left + (I_NumberOfCells + 1)*D_Delta_x ) )/2.0;

  for (int i = 1; i <= I_NumberOfCells; ++i)
    {
      DA_P_Values_P[i] = ( Value_P_0( D_Boundary_Position_Left + ( i-1 )*D_Delta_x ) 
			   + Value_P_0( D_Boundary_Position_Left + ( i )*D_Delta_x ) )/2.0;
    }
}

//	Definition of the function that provides the integral average of the friction coefficient for the cells
void SemilinearSystem::Set_Lambda_Averages( double* DA_P_Lambda_AV, double* DA_P_Lambda_Coefficients,
					    int I_Lambda_Expansion_Length, int I_NumberOfCells, double D_Boundary_Position_Left, double D_Delta_x  )
{	
	/**
	*	@brief	calculates the integral averages of the friction coefficient in the computational cells
	*
	*	The function provides integral averages of the friction coeffficient in the compuational cells of the discretization.\n
	*	The arguments are:
	*
	*	@param[in]	 	DA_P_Lambda_AV				Memory position of the array, the cell-average of the friction coefficient is stored in
	*	@param[in]		DA_P_Lambda_Coefficients	Memary position of the array containing the coefficients of the expansion of Lambda
	*	@param[in]		I_Lambda_Expansion_Length	Integer at which the taylor expansion of the random field is truncated (length(DA_P_Lambda_Coefficients)-1)
	*	@param[in]		I_NumberOfCells				Number of cells in the discretization
	*	@param[in]		D_Boundary_Position_Left	Coordinate of the left boundary of the pipe
	*	@param[in]		D_Delta_x					Width of the spatial discretization
	*
	*	@date 07.09.2017
	*
	*	@version 1.0
	*		Tested for constant Lambda, i.e., the coefficient array contains one element only.
	*
	*	@todo	better integration routine instead of average
	*/
	//	Initialization of the increment double
	double	Increment;
	//	Left boundary of the computaing cell
	double 	xx_L					=	D_Boundary_Position_Left;
	//	right boundary of the computaing cell
	double	xx_R					=	D_Boundary_Position_Left + D_Delta_x;
	//	The method computes the integral average of the friction coefficient
	//	FIRST VERSION: Simple average of Lamdba at the boudnary of the corresponding integral.
	//	FURTERH PLANS: Increase the accuracy of the integration process. Since Lambda is expected to be in L^\infty, discontinuities will happem!

	int counter = 0;
	
	//	Calculating the integral average on the cells within the pipe.
	//	Loop over interior cells
	for (int i = 1; i <= I_NumberOfCells; ++i)
	{
		//	Initialization of the integral average in cell 'i'
		DA_P_Lambda_AV[i]	=	0.0;
		//	Loop over coefficients of the expansion

		//	j = 0
		Increment				=	( Value_Lambda_Base(0, xx_L) + Value_Lambda_Base(0, xx_R) )/2.0;
		//	Update of Lamdba-Average
		DA_P_Lambda_AV[i]		=	DA_P_Lambda_AV[i] + DA_P_Lambda_Coefficients[0]*Increment;

		for (int j = 1; j < I_Lambda_Expansion_Length; j+=2)
		{
		  counter += 1;
		//	Evaluating the average of basis function 'j' in interval 'i'
		  Increment			=	( Value_Lambda_Base_Sin(counter, xx_L, domain_length) + Value_Lambda_Base_Sin(counter, xx_R, domain_length) )/2.0;
		//	Update of Lamdba-Average
		  DA_P_Lambda_AV[i]		=	DA_P_Lambda_AV[i] + DA_P_Lambda_Coefficients[j]*Increment * ( 1.0 / pow(counter,2) );
		}
		counter = 0;
		
		/* for cos */
		for (int j = 2; j < I_Lambda_Expansion_Length + 1; j+=2)
		{
		  counter += 1;
		  //	Evaluating the average of basis function 'j' in interval 'i'
		  Increment			=	( Value_Lambda_Base_Cos(counter, xx_L, domain_length) + Value_Lambda_Base_Cos(counter, xx_R, domain_length) )/2.0;
		  //	Update of Lamdba-Average
		  DA_P_Lambda_AV[i]		=	DA_P_Lambda_AV[i] + DA_P_Lambda_Coefficients[j]*Increment * ( 1.0 / pow(counter,2) );
		}

		/**
		 * A bug was fixed using this script.
		  if (abs(DA_P_Lambda_AV[i]) > 1000.0){
		    cout << "error in Set_... : " << I_Lambda_Expansion_Length << endl;
		    for (int j=1; j < I_Lambda_Expansion_Length + 1; j+=2){ cout << "j: " << j << "," << DA_P_Lambda_Coefficients[j] << " "; }
		    for (int j=2; j < I_Lambda_Expansion_Length + 1; j+=2){ cout << "j: " << j << "," << DA_P_Lambda_Coefficients[j] << " "; }
		    cout << endl;
		  }
		*/

		counter = 0;

		//	Re-Defining the values of the interval boundaries
		xx_L						=	xx_R;
		xx_R						=	xx_R + D_Delta_x;
		
	}
	
	//	Fixing the friction coefficient on the ghost cells on left and right boundary of the pipeline
	//	HEURISTIK: Linear Interpolation of 2 adjacent interior cells and maximum with 0 to ensure non-negativity of the friction coefficient
	DA_P_Lambda_AV[0]					=	max( 0.0, 2.0*DA_P_Lambda_AV[1] - DA_P_Lambda_AV[2]);
	DA_P_Lambda_AV[I_NumberOfCells + 1]	=	max( 0.0, 2.0*DA_P_Lambda_AV[I_NumberOfCells] - DA_P_Lambda_AV[I_NumberOfCells - 1]);
}


void SemilinearSystem::Run(double DA_P_Lambda_Coefficients[], bool write2file_bool, bool progress_bool)
{
  /**
   *	@brief function that solves the PDE-System for a
   *	friction coefficient defined by the coefficients
   *	passed by caller and storing the boundary values of
   *	the unknown P at all time steps
   *
   *	This member function is the essential handshake for
   *	the user.\n Here, the coefficients of the
   *	Lambda-Expansion are provided, the semilienar system
   *	is solved and the boundary values of P are stored in
   *	the corresponding arrays.\n This method calls:\n
   *	Set_Lambda_Averages' with the provided coefficients\n
   *	
   *	@param[in] DA_P_Lambda_Coefficients Coefficients of
   *	the Lambda-Expansion with the expectation value of
   *	lambda on the first array element (length = Number of
   *	Coefficients + 1)
   *
   *	@date 07.09.2017
   *
   *	@version 1.0
   *
   * @todo check the scheme again! - done! (11.09.2017)\n test
   *					without source-term
   *					for traveling waves
   */
  
  //	Reading the initial conditions for P and Q
  Set_IC_Q(DA_P_Values_Q, I_NumberOfCells, D_Boundary_Position_Left, D_Delta_x);	
  Set_IC_P(DA_P_Values_P, I_NumberOfCells, D_Boundary_Position_Left, D_Delta_x);

  /* filling the initial pressure boundary data */
  DA_P_Left_P[0]                    = DA_P_Values_P[0];
  DA_P_Right_P[0] = DA_P_Values_P[I_NumberOfCells+1];
  
  //	Declaring arrays for P,Q that store intermediate results
  double DA_Values_P_Intermediate[I_NumberOfCells + 2];
  double DA_Values_Q_Intermediate[I_NumberOfCells + 2];
  
  //	Calling the method that sets the integral averages of the chosen friction coefficient
  Set_Lambda_Averages(DA_P_Lambda_AV, DA_P_Lambda_Coefficients,
		      I_Lambda_Expansion_Length, I_NumberOfCells,
		      D_Boundary_Position_Left, D_Delta_x);
  
  
  // setting time to zero
  D_Current_T = 0.0;
  I_Current_T_Index = 0;
  
  // writing the values
  char filename[100] = "output";
  if (write2file_bool){
    Write2File(filename, false);
  }

  //	Solving the forward problem with current friction coefficient
  //	Loop over timesteps of the discretization
  for (int I_TimeStepCount = 1; I_TimeStepCount <= I_NumberOfTimeSteps; ++I_TimeStepCount)//1; ++I_TimeStepCount)//
    {

      if ( I_TimeStepCount % 50 == 0 && progress_bool){
	cout << "computing the solution " << std::fixed 
	     << std::setprecision(2) 
	     << 100 * I_TimeStepCount / float(I_NumberOfTimeSteps) << "%" << " from cpp" << "\r";
      }

      //	Updating P and Q on the interior cells of the
      //	discretization.  Loop over interior cells of the
      //	discretization (all except '0' and 'I_NumberOfCells +
      //	1' reserved for the ghost cells)
      for (int I_CellCount = 1; I_CellCount <= I_NumberOfCells; ++I_CellCount)
	{
	  //	Calculating the value of P
	  DA_Values_P_Intermediate[I_CellCount] = (1.0/2.0)*( DA_P_Values_P[I_CellCount - 1] + DA_P_Values_P[I_CellCount + 1] )
	    - (1.0/(2.0*D_SpeedOfSound))*( DA_P_Values_Q[I_CellCount + 1] - DA_P_Values_Q[I_CellCount - 1] )
	    + (D_Delta_t/(2.0*D_SpeedOfSound))*( DA_P_Lambda_AV[I_CellCount-1] *
						 FrictionFunction( DA_P_Values_P[I_CellCount - 1], DA_P_Values_Q[I_CellCount - 1] )
						 - DA_P_Lambda_AV[I_CellCount + 1] *
						 FrictionFunction( DA_P_Values_P[I_CellCount + 1], DA_P_Values_Q[I_CellCount + 1] ) );
	  //	Calculating the value of Q
	  DA_Values_Q_Intermediate[I_CellCount] = (1.0/2.0)*( DA_P_Values_Q[I_CellCount - 1] + DA_P_Values_Q[I_CellCount + 1] )
	    - (D_SpeedOfSound/2.0)*( DA_P_Values_P[I_CellCount + 1] - DA_P_Values_P[I_CellCount - 1] )
	    + (D_Delta_t/2)*( DA_P_Lambda_AV[I_CellCount-1] *
			      FrictionFunction( DA_P_Values_P[I_CellCount - 1], DA_P_Values_Q[I_CellCount - 1] )
			      + DA_P_Lambda_AV[I_CellCount + 1] *
			      FrictionFunction( DA_P_Values_P[I_CellCount + 1], DA_P_Values_Q[I_CellCount + 1] ) );
	}
		
      //	Updating Q on the ghost cells by calling the functions
      //	specifying the boudnary values for Q
      DA_Values_Q_Intermediate[0]			=	Value_Q_L( I_TimeStepCount*D_Delta_t );
      DA_Values_Q_Intermediate[I_NumberOfCells + 1]	=	Value_Q_R( I_TimeStepCount*D_Delta_t );
		
      //	Updating P on the ghost cells according to the flow of P and Q
      DA_Values_P_Intermediate[0]	=	DA_P_Values_P[1] - DA_P_Values_Q[1]/D_SpeedOfSound 
	+ DA_Values_Q_Intermediate[0]/D_SpeedOfSound 
	- (D_Delta_t/D_SpeedOfSound)*DA_P_Lambda_AV[1]*FrictionFunction( DA_P_Values_P[1], DA_P_Values_Q[1] );

      DA_Values_P_Intermediate[I_NumberOfCells + 1]	=	DA_P_Values_P[I_NumberOfCells] + DA_P_Values_Q[I_NumberOfCells]/D_SpeedOfSound 
	- DA_Values_Q_Intermediate[I_NumberOfCells + 1]/D_SpeedOfSound 
	+ (D_Delta_t/D_SpeedOfSound)*DA_P_Lambda_AV[I_NumberOfCells]
	*FrictionFunction( DA_P_Values_P[I_NumberOfCells], DA_P_Values_Q[I_NumberOfCells] );
	        
      
      
      //	Overwriting 'DA_P_Values_Q' and 'DA_P_Values_P'
      for (int I_CellCount = 0; I_CellCount <= I_NumberOfCells + 1; ++I_CellCount)
	{
	  //	Calculating the value of P
	  DA_P_Values_P[I_CellCount]					=	DA_Values_P_Intermediate[I_CellCount];
	  //	Calculating the value of Q
	  DA_P_Values_Q[I_CellCount]					=	DA_Values_Q_Intermediate[I_CellCount];
	}

      // Storing the boundary values of P in the corresponding vectors
      DA_P_Left_P[I_TimeStepCount] = 0.0;
      DA_P_Right_P[I_TimeStepCount] = 0.0;
	
      for (int i=0; i<=n_epsilon; i++){
	DA_P_Left_P[I_TimeStepCount]  += pow( DA_P_Values_P[i], 2);
	DA_P_Right_P[I_TimeStepCount] += pow( DA_P_Values_P[I_NumberOfCells + 1 - i], 2);
      }
      
      DA_P_Left_P[I_TimeStepCount] = sqrt(DA_P_Left_P[I_TimeStepCount] * D_Delta_x);
      DA_P_Right_P[I_TimeStepCount] = sqrt(DA_P_Right_P[I_TimeStepCount] * D_Delta_x);

      /* pointwise storage of P */
      // DA_P_Left_P[I_TimeStepCount]		=	DA_Values_P_Intermediate[0];
      // DA_P_Right_P[I_TimeStepCount]		=	DA_Values_P_Intermediate[I_NumberOfCells + 1];

      // upodate the time
      D_Current_T = D_Current_T + D_Delta_t;
      I_Current_T_Index += 1;
      // writing the values
      if (write2file_bool){
	Write2File(filename, true);
      }
     
    }

  if (progress_bool){
    cout << endl;			// move the cursor to the next line
  }

}

void SemilinearSystem::Write2File(char* filename, bool append)
{
  /** 
   * @brief writing quantities into a file.
   * @author Soheil Hajian
   */

  //	Creating instances of the Data-File-Operation - Classes
  fstream	P_Store;
  fstream	Q_Store;

  char urlP[100];
  char urlQ[100];

  strcpy(urlP, filename);
  strcpy(urlQ, filename);  

  strcat(urlP, "_P.dat");
  strcat(urlQ, "_Q.dat");
  
  //	Open the objects
  if (append)
    {
      P_Store.open(urlP, fstream::app|fstream::out);
      Q_Store.open(urlQ, fstream::app|fstream::out);
    }
  else
    {
      P_Store.open(urlP, std::ios_base::out);
      Q_Store.open(urlQ, std::ios_base::out);
    }

  for (int i=0; i<I_NumberOfCells; i++)
    {
      P_Store << DA_Centroid[i] << "," << DA_P_Values_P[i+1] << "," << D_Current_T << endl;
      Q_Store << DA_Centroid[i] << "," << DA_P_Values_Q[i+1] << "," << D_Current_T << endl;
    }

  //	closing
  P_Store.close();
  Q_Store.close();
    
}


void SemilinearSystem::Run_Text_Output( double DA_P_Lambda_Coefficients[] )
{
  /**
   *	@brief function that solves the PDE-System for a friction
   *	coefficient defined by the coefficients passed by caller and
   *	storing the boundary values of the unknown P at all time
   *	steps.
   *
   *    Moreover, the values of P and Q are stored in a Matlab-readable format.
   *
   *	This member function is the essential handshake for the
   *	user.\n Here, the coefficients of the Lambda-Expansion are
   *	provided, the semilienar system is solved and the boundary
   *	values of P are stored in the corresponding arrays.\n All
   *	intermediate values of P and Q are stored in the
   *	Matlab-readable form Pressure.m and VolumeFlow.m. The call
   *	from Matlab has the form P = Pressure(); and Q =
   *	VolumeFlow();\n ensuring, that for example the Matrix P
   *	((Number Of Time Steps + 1)x(Number of Cells + 2)) contains
   *	all values of the computation including values in the ghost
   *	cells and initial state.\n This method calls:\n
   *	Set_Lambda_Averages' with the provided coefficients\n
   *	
   *	@param[in] DA_P_Lambda_Coefficients Coefficients of the
   *	Lambda-Expansion with the expectation value of lambda on the
   *	first array element (length = Number of Coefficients + 1)
   *
   *	@date 			21.09.2017
   *
   *	@version 		1.0
   *
   *	@todo check the scheme again! - done! (11.09.2017)\n test
   *					without source-term for
   *					traveling waves
   */
	
	
	
	//	Reading the initial conditions for P and Q
	Set_IC_Q(DA_P_Values_Q, I_NumberOfCells, D_Boundary_Position_Left, D_Delta_x);	
	Set_IC_P(DA_P_Values_P, I_NumberOfCells, D_Boundary_Position_Left, D_Delta_x);
	
	//	Creating instances of the Data-File-Operation - Classes
	fstream	P_Store;
	fstream	Q_Store;
	
	//	Open the objects
	P_Store.open("Pressure.m", ios::out);
	Q_Store.open("VolumeFlow.m", ios::out);
	
	
	//	Initial lines in the text files to make the call-able by matlab
	P_Store << "function [ P ] = Call_Pressure( input_args )" << endl;
	P_Store	<< " P = [ " << endl;
	
	Q_Store << "function [ Q ] = Call_VolumeFlow( input_args )" << endl;
	Q_Store	<< " Q = [ " << endl;
	
	
	//	Storing Initial Data
	for (int I_CellCount = 0; I_CellCount <= I_NumberOfCells + 1; ++I_CellCount)
	{
		//	Storing P
		P_Store	<<  DA_P_Values_P[I_CellCount] << "  ";
		//	Storing Q
		Q_Store <<  DA_P_Values_Q[I_CellCount] << "  ";
 	}
 	P_Store << ";" << endl;
 	Q_Store << ";" << endl;
 	
	
	//	Declaring arrays for P,Q that store intermediate results
	double DA_Values_P_Intermediate[I_NumberOfCells + 2];
	double DA_Values_Q_Intermediate[I_NumberOfCells + 2];
	
	//	Calling the method that sets the integral averages of the chosen friction coefficient
	Set_Lambda_Averages(DA_P_Lambda_AV, DA_P_Lambda_Coefficients, I_Lambda_Expansion_Length, I_NumberOfCells, D_Boundary_Position_Left, D_Delta_x);
	
	//	Solving the forward problem with current friction coefficient
	//	Loop over timesteps of the discretization
	for (int I_TimeStepCount = 1; I_TimeStepCount <= I_NumberOfTimeSteps; ++I_TimeStepCount)//1; ++I_TimeStepCount)//
	{
		//	Updating P and Q on the interior cells of the discretization
		//	Loop over interior cells of the discretization (all except '0' and 'I_NumberOfCells + 1' reserved for the ghost cells)
		for (int I_CellCount = 1; I_CellCount <= I_NumberOfCells; ++I_CellCount)
		{
		  //	Calculating the value of P
		  DA_Values_P_Intermediate[I_CellCount]		=	(1.0/2.0)*( DA_P_Values_P[I_CellCount - 1] + DA_P_Values_P[I_CellCount + 1] ) - (1.0/(2.0*D_SpeedOfSound))*( DA_P_Values_Q[I_CellCount + 1] - DA_P_Values_Q[I_CellCount - 1] ) + (D_Delta_t/(2.0*D_SpeedOfSound))*( DA_P_Lambda_AV[I_CellCount-1]*FrictionFunction( DA_P_Values_P[I_CellCount - 1], DA_P_Values_Q[I_CellCount - 1] ) - DA_P_Lambda_AV[I_CellCount + 1]*FrictionFunction( DA_P_Values_P[I_CellCount + 1], DA_P_Values_Q[I_CellCount + 1] ) );
		  //	Calculating the value of Q
		  DA_Values_Q_Intermediate[I_CellCount]		=	(1.0/2.0)*( DA_P_Values_Q[I_CellCount - 1] + DA_P_Values_Q[I_CellCount + 1] ) - (D_SpeedOfSound/2.0)*( DA_P_Values_P[I_CellCount + 1] - DA_P_Values_P[I_CellCount - 1] ) + (D_Delta_t/2)*( DA_P_Lambda_AV[I_CellCount-1]*FrictionFunction( DA_P_Values_P[I_CellCount - 1], DA_P_Values_Q[I_CellCount - 1] ) + DA_P_Lambda_AV[I_CellCount + 1]*FrictionFunction( DA_P_Values_P[I_CellCount + 1], DA_P_Values_Q[I_CellCount + 1] ) );
 		}
		
		//	Updating Q on the ghost cells by calling the functions specifying the boudnary values for Q
		DA_Values_Q_Intermediate[0]						=	Value_Q_L( I_TimeStepCount*D_Delta_t );
		DA_Values_Q_Intermediate[I_NumberOfCells + 1]	=	Value_Q_R( I_TimeStepCount*D_Delta_t );
		
		//	Updating P on the ghost cells according to the flow of P and Q
		DA_Values_P_Intermediate[0]						=	DA_P_Values_P[1] - DA_P_Values_Q[1]/D_SpeedOfSound + DA_Values_Q_Intermediate[0]/D_SpeedOfSound - (D_Delta_t/D_SpeedOfSound)*DA_P_Lambda_AV[1]*FrictionFunction( DA_P_Values_P[1], DA_P_Values_Q[1] );
		DA_Values_P_Intermediate[I_NumberOfCells + 1]	=	DA_P_Values_P[I_NumberOfCells] + DA_P_Values_Q[I_NumberOfCells]/D_SpeedOfSound - DA_Values_Q_Intermediate[I_NumberOfCells + 1]/D_SpeedOfSound + (D_Delta_t/D_SpeedOfSound)*DA_P_Lambda_AV[I_NumberOfCells]*FrictionFunction( DA_P_Values_P[I_NumberOfCells], DA_P_Values_Q[I_NumberOfCells] );
		
		//	Storing the boundary values of P in the corresponding vectors
		DA_P_Left_P[I_TimeStepCount]					=	DA_Values_P_Intermediate[0];
		DA_P_Right_P[I_TimeStepCount]					=	DA_Values_P_Intermediate[I_NumberOfCells + 1];
		
		//	Overwriting 'DA_P_Values_Q' and 'DA_P_Values_P'
		for (int I_CellCount = 0; I_CellCount <= I_NumberOfCells + 1; ++I_CellCount)
		{
			//	Calculating the value of P
			DA_P_Values_P[I_CellCount]					=	DA_Values_P_Intermediate[I_CellCount];
			
			DA_Values_P_Intermediate[I_CellCount] 		=	0.0;
			//	Calculating the value of Q
			DA_P_Values_Q[I_CellCount]					=	DA_Values_Q_Intermediate[I_CellCount];
			DA_Values_Q_Intermediate[I_CellCount] 		=	0.0;
			
			//	Storing P
			P_Store	<<  DA_P_Values_P[I_CellCount] << "  ";
			//	Storing Q
			Q_Store <<  DA_P_Values_Q[I_CellCount] << "  ";
 		}
	 	P_Store << ";" << endl;
	 	Q_Store << ";" << endl;
	 	
	 	
	}
	P_Store << "];" << endl;
	Q_Store << "];" << endl;
	
	//	closing
	P_Store.close();
	Q_Store.close();
	
}



double SemilinearSystem::get_P_Left_Boundary( double D_EvaluationTime_Arg )
{
	/**
	*	@brief returns the trace evaluation of P at the left
	*	boundary at a user-specified time
	*
	*	The member function provides the boundary value of P
	*	for a user-spcified poiint in time at the left
	*	boundary of the pipe.  This trace evaluation is based
	*	on the following basic assumptions.\n First, function
	*	value of P at the boundary (atctually within the
	*	entire pipe) is piecewise constant with respect to
	*	time in that for every time step n, we have P(t,0) =
	*	P(n*D_Delta_t, 0) for t \in ((n-1/2)D_Delta_t,
	*	(n+1/2)D_Delta_t].\n Second, the trace evaluation is
	*	approximated by the value of P on the left ghost
	*	cell.\n To find the integer value corresponding to the
	*	array-entry, D_EvaluationTime_Arg/D_Delta_t is rounded
	*	to the next integer value.
	*
	*	@param[in] D_EvaluationTime_Arg Time at which the
	*	boundary trace of P at the left boundary is evaluated
	*	@param[out] D_BoundaryValue Boundary value of P at the
	*	left boundary
	*
	*	@date		09.09.2017
	*
	*	@version	1.0
	*
	*	@todo		Schreiben der Auswertungsfunktion
	*/
	
	int I_Element = (int)( D_EvaluationTime_Arg/D_Delta_t + 0.5 );	//	Rounding to the next integer value
	
	double D_BoundaryValue = DA_P_Left_P[I_Element];				//	Reading the corresponing element
	
	return D_BoundaryValue;											//	Return the result
	
}


double SemilinearSystem::get_P_Right_Boundary( double D_EvaluationTime_Arg )
{
  /**
   *	@brief		returns the trace evaluation of P at the right boundary at a user-specified  time
   *
   *	The member function provides the boundary value of P
   *	for a user-spcified point in time at the right
   *	boundary of the pipe.  This trace evaluation is based
   *	on the following basic assumptions.\n First, function
   *	value of P at the boundary (atctually within the
   *	entire pipe) is piecewise constant with respect to
   *	time in that for every time step n, we have P(t,L) =
   *	P(n*D_Delta_t, L) for t \in ((n-1/2)D_Delta_t,
   *	(n+1/2)D_Delta_t].\n Second, the trace evaluation is
   *	approximated by the value of P on the right ghost
   *	cell.
   *
   *	@param[in] D_EvaluationTime_Arg Time at which the
   *	boundary trace of P at the right boundary is evaluated
   *	@param[out] D_BoundaryValue Boundary value of P at the
   *	right boundary
   *
   *	@date		09.09.2017
   *
   *	@version	1.0
   *
   *	@todo		Schreiben der Auswertungsfunktion
   */
	
	
	int I_Element = (int)( D_EvaluationTime_Arg/D_Delta_t + 0.5 );	//	Rounding to the next integer value
	
	double D_BoundaryValue = DA_P_Right_P[I_Element];				//	Reading the corresponing element
	
	return D_BoundaryValue;											//	Return the result
}


double SemilinearSystem::get_P_Difference( double D_EvaluationTime_Arg )
{
  /**
   *	@brief returns difference of the trace evaluation of P at
   *	right and left boundary at a user-specified time
   *
   *	The member function provides the difference of the boundary
   *	values of P for a user-spcified poiint in time between right
   *	and left boundary of the pipe.  This trace evaluation is based
   *	on the following basic assumptions.\n First, function value of
   *	P at the boundary (atctually within the entire pipe) is
   *	piecewise constant with respect to time in that for every time
   *	step n, we have P(t,0) = P(n*D_Delta_t, 0) for t \in
   *	((n-1/2)D_Delta_t, (n+1/2)D_Delta_t].\n Second, the trace
   *	evaluation is approximated by the value of P on the left ghost
   *	cell.\n To find the integer value corresponding to the
   *	array-entry, D_EvaluationTime_Arg/D_Delta_t is rounded to the
   *	next integer value.
   *
   *	@param[in] D_EvaluationTime_Arg Time at which the boundary
   *	trace of P at the left boundary is evaluated
   *	@param[out] D_BoundaryDifference Boundary value of P at the
   *	left boundary
   *
   *	@date		09.09.2017
   *
   *	@version	1.0
   *
   *	@todo		Schreiben der Auswertungsfunktion
   */
	
  int I_Element	= (int)( D_EvaluationTime_Arg/D_Delta_t + 0.5 );		//	Rounding to the next integer value
  
  double D_BoundaryDifference = DA_P_Right_P[I_Element] - DA_P_Left_P[I_Element]; //	Reading the corresponing element

  return D_BoundaryDifference;	//	Return the result
	
}


void SemilinearSystem::EvalTest(   )
{
	/**
	*	@brief	function for testing purposes.
	*
	*/

	
	Set_IC_Q(DA_P_Values_Q, I_NumberOfCells, D_Boundary_Position_Left, D_Delta_x);
	
	Set_IC_P(DA_P_Values_P, I_NumberOfCells, D_Boundary_Position_Left, D_Delta_x);
	

	
	for (int i = 0; i < I_Lambda_Expansion_Length + 1; ++i)
	{
		DA_P_Lambda_Coefficients[i] = 0.3;
	}

	
	Set_Lambda_Averages(DA_P_Lambda_AV, DA_P_Lambda_Coefficients, I_Lambda_Expansion_Length, I_NumberOfCells, D_Boundary_Position_Left, D_Delta_x);
	
	//	return the 0!!!!
	return; 
}

void SemilinearSystem::info()
{
  /**
   * @brief giving info about the pipe in the terminal
   * @author Soheil Hajian
   */
  cout << setfill ('=') << setw(50) << right << " info ==" << endl;
  cout << setfill ('_') << setw(35) << left << "Omega "
       << left << "(" << D_Boundary_Position_Left << "," << D_Boundary_Position_Right << ")" << endl;
  cout << setw(35) << left << "Length of the domain "
       << left << domain_length << endl;
  cout << setw(35) << left << "Dt "
       << left << D_Delta_t << endl;
  cout << setw(35) << left << "Dx "
       << left << D_Delta_x << endl;
  cout << setw(35) << left << "Number of cells "
       << left << I_NumberOfCells << endl;
  cout << setw(35) << left << "Current time "
       << left << D_Current_T << endl;
  cout << setw(35) << left << "Final time "
       << left << D_TerminalTime << endl;
  cout << setw(35) << left << "Epsilon for the boundary "
       << left << epsilon_boundary << endl;
  cout << setw(35) << left << "N_epsilon "
       << left << n_epsilon << endl;
  cout << setw(35) << left << "Size of friction coef. vec. "
       << left << I_Lambda_Expansion_Length << endl;
  cout << setfill ('=') << setw(50) << right << "" << endl;
}

int SemilinearSystem::CurrentTimeIndex()
{
  /**
   * @brief returns current time index.
   * @author Soheil Hajian
   */
  return I_Current_T_Index;
}

double* SemilinearSystem::BoundaryValueP_Left()
{
  /**
   * @brief returns the Boundary Value of P as an array.
   * @author Soheil Hajian
   */

  return DA_P_Left_P;
}

double* SemilinearSystem::BoundaryValueP_Right()
{
  /**
   * @brief returns the Boundary Value of P as an array.
   * @author Soheil Hajian
   */

  return DA_P_Right_P;
}

double* SemilinearSystem::TimeSlices()
{
  /**
   * @brief returns the time slices
   * @author Soheil Hajian
   */

  return time_slices;
}

double* SemilinearSystem::LambdaAverage(double DA_P_Lambda_Coefficients_GIVEN[])
{
  /**
   * return compute and return LambdaAverage
   */

  Set_Lambda_Averages( DA_P_Lambda_AV, DA_P_Lambda_Coefficients_GIVEN,
		       I_Lambda_Expansion_Length, I_NumberOfCells,
		       D_Boundary_Position_Left, D_Delta_x  );

  return DA_P_Lambda_AV;

}

double* SemilinearSystem::SendLambdaAverage()
{
  /**
   * send the current DA_P_Lambda_AV: It is updated after every SemilinearSystem::Run
   */
  return DA_P_Lambda_AV;
}

int SemilinearSystem::NumberofCells()
{
  /** 
   * send back the number of cells
   */
  return I_NumberOfCells;
  
}


