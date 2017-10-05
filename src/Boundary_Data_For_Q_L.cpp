
#include "Boundary_Data_For_Q_L.h"


double Value_Q_L( double Time )
{
	/*
	Testing Purposes
	cout << 120.0 + 40.0*sin( 5.0*Time ) << " Bounadry Q L " << Time << "  Time " << sin( Time ) << " Value of sin "  <<  endl;
	return 120.0 + 40.0*sin( Time ); //;
	*/
	return 180.0 + pow(Time,2);
}
