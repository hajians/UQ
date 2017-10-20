
#include "Initial_Data_For_Q.h"

const double pi = atan2(0,-1);

double Value_Q_0( double coordinate )
{
  //	return 120.0;
  return 10.0 + sin( 2*pi*coordinate );
}
