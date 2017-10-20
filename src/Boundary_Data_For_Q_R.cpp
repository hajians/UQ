
#include "Boundary_Data_For_Q_R.h"

const double pi = atan2(0,-1);

double Value_Q_R( double Time )
{
  //return 0.5;
  return 10.0 + sin( -2.0*pi*Time );
  // // Testinf purposes
  // double ReturnValue;
  // if (Time < 1/300.0)
  //   {
  //     ReturnValue = 120.0;
  //   }
  // else
  //   {
  //     ReturnValue = 24*sin( 10*Time );
  //    }

  // ReturnValue = 24*sin( 2*Time );
  
  // 	//cout << ReturnValue << " Bounadry Q R " << Time - 1/300.0 << "  Time " << sin( 5.0*Time - 1/300.0 ) << " Value of sin "   << endl;
	
  // return ReturnValue; //return 120.0 + sin( Time ); //
  // //return 180 + pow(Time,2);
}
