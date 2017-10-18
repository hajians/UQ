
#include "SourceTerm.h"


double FrictionFunction( double D_Value_P, double D_Value_Q )
{
  //	return  -abs(D_Value_Q)*D_Value_Q/D_Value_P;			//abs(D_Value_P*D_Value_Q);
  return  -abs(D_Value_Q)*D_Value_Q;
}
