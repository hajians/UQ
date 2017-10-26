
#include "Value_Lambda_Base.h"

const double pi = atan2(0,-1);

double Value_Lambda_Base( int Base_Function_Number, double Coordinate )
{
	return 1.0;
}

double Value_Lambda_Base_Sin(int Base_Function_Number, double Coordinate, double length)
{
  return sin(2.0 * pi * Base_Function_Number * Coordinate / length);
}

double Value_Lambda_Base_Cos(int Base_Function_Number, double Coordinate, double length)
{
  return cos(2.0 * pi * Base_Function_Number * Coordinate / length);
}
