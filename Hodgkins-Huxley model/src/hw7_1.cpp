#include "hw7.h"
#include<cmath>
#include<functional>
#include<vector>
using namespace std;

#define EK -80
#define STEPSIZE 0.01

void K_Euler( vector<double> &current, vector<double> &vol, int time_length )
{
	// Potassium channel simulation using Forward Euler solver

	// initialize
	vol.push_back(0);
	double n = 0;

	// create the solver for ODEs
	function<double (double,double)> fn = f_n;
	Euler<double> euler_dn( STEPSIZE, fn );
	function<double (double,double)> fv = f_v_K;
	Euler<double> euler_dv( STEPSIZE, fv );

	// iterate through the time period and collect the current and voltage
	for( double h = 0, i = 0; h < time_length; h += STEPSIZE, ++i ){
		n = euler_dn.solve( vol[i], n );
		current.push_back( pow(n, 4)*(vol[i] - EK) );
		vol.push_back( euler_dv.solve( n, vol[i]) );
	}
}