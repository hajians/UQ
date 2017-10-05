/**
*	@brief 		Header-File of user-provided function specifying the source term of the semilinear system.
*
*				The function specifies the friction function in the semilinear model of gas transport through a single pipe.\n
*				The function has to be specified by the user and assigns the determininstic influence of Pressure and Volume flow
*				to the latter values. The interface is FrictionFunction, takes two double values and returns the corresponding
*				value of the friction function.\n
*				The inital choice is\n
*				<center>
*					a(P,Q)	=	abs(P,Q)
*				</center>
*
*	@param[in]	D_Value_P			Value of the Pressure P
*	@param[in]	D_Value_Q			Value of the Volume Flow Q
*	@param[out]	D_FrictionValue		Value of the user-provided friction function
*
*	@author		Dr. Nikolai Strogies
*
*	@date		08.09.2017
*
*	@version	1.0
*
*	@todo		different choices for the friction function
*	
*/


#include <iostream>

#include <cmath>

using namespace std;




double FrictionFunction( double D_Value_P, double D_Value_Q );
