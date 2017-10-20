
#include "Boundary_Data_For_Q_L.h"

const double pi = atan2(0,-1);

double Value_Q_L( double Time )
{
	
  // Testing Purposes
  // cout << 120.0 + 40.0*sin( 5.0*Time ) << " Bounadry Q L " << Time << "  Time " << sin( Time ) << " Value of sin "  <<  endl;
  return 10.0 + sin( -2*pi*Time );
  //return 1.0;
    // 120.0 + 40.0*sin( Time ); //;
    // return 180.0 + pow(Time,2);
}
