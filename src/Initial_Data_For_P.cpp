
#include "Initial_Data_For_P.h"

const double pi = atan2(0,-1);

double Value_P_0( double coordinate )
{
  // return sqrt( pow(50.0,2) ); //- 2*0.3*120.0*coordinate );
  return 10.0 + sin( 2*pi*coordinate );
}
