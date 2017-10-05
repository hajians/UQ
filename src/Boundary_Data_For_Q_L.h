/*
	This function specifies the boundary data for Q
	It is a continuous function that, for a given
					double Time
	provides the value of Q at the left boudnary
	
	
	
	The default value is
					Q(t) = 180 + t^2
*/


#include <iostream>

#include <cmath>

using namespace std;




double Value_Q_L( double Time );
