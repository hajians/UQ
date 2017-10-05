
#include "Boundary_Data_For_Q_R.h"


double Value_Q_R( double Time )
{
	/*
	Testinf purposes
	double ReturnValue;
	if (Time < 1/300.0)
	{
		ReturnValue = 120.0;
	}
	else
	{
		ReturnValue = 120.0 + 40.0*sin( Time - 1/300.0 );
	}
	
	cout << ReturnValue << " Bounadry Q R " << Time - 1/300.0 << "  Time " << sin( 5.0*Time - 1/300.0 ) << " Value of sin "   << endl;
	
	return ReturnValue; //return 120.0 + sin( Time ); //
	*/
	return 180 + pow(Time,2);
}
