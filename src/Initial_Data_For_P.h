/**
*	@brief 		Header-File of user-provided function specifying the initial data for P.
*
*				The function specifies the initial data for P in the semilinear model of gas transport through a single pipe.\n
*				The inital choice is\n
*				<center>
*					P(x)	=	sqrt( P(0)^2  )
*				</center>
*
*	@param[in]	Coordinate		Double that represents the coordiante in the spatial domain
*	@param[out]	Value_P_0		Double that represents the value of the initial datum for P at the provided coordinate.
*
*	@author		Dr. Nikolai Strogies
*
*	@date		08.09.2017
*
*	@version	1.0
*	
*/


#include <iostream>

#include <cmath>

using namespace std;




double Value_P_0( double Coordinate );
