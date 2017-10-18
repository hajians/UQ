/*This Code evaluates the forward problem of gas transport through a
  single pipe given by the system of PDE's

  p_t + 	  q_x  = 0
  q_t + c^2 p_x  = \lambda f(p,q)

where c denotes the (constant) speed of sound. The semilinear system
not only depends on the spatially distributed friction parameter
\lambda but also on a friction function f: R^2->R^1 that simulates
nonlinear effects.*/

#include <iostream>

// #include <conio.h>

#include "SemilinearSystem.h"

using namespace std;

int main()
{
  //  SemilinearSystem(double SpeedOfSound, double T, double x_L,
  //  double x_R, double Dx, int ExNum)
  // SemilinearSystem SemilinearSystem(1.0, 10.0, 0.0, 1.0, .01, 0);

  SemilinearSystem pipe1 = SemilinearSystem(1.0, 10.0, 0.0, 1.0, .01, 0);
  SemilinearSystem pipe2 = SemilinearSystem(1.0, 10.0, 0.0, 1.0, .01, 0);

  pipe1.info();

  // SemilinearSystem.EvalTest( );
	
  //	Run over time steps
  // SemilinearSystem.Run()

  cout<< "> Starting <" << endl;
  
  double a[1] = {0.0};
  pipe1.Run(a);
	
  // SemilinearSystem.Run(a);

  // SemilinearSystem.Write2File("output", false);

	
  // getch ();
	
  // cout << SemilinearSystem.get_P_Left_Boundary( 0.5 ) << endl;
	
  // double b[1] = {.7};
	
  // SemilinearSystem.Run(b);
	
	
  // cout << SemilinearSystem.get_P_Left_Boundary( 0.5 ) << endl;
  
  return 0;
}
