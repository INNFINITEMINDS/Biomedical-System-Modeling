#include "hw7.h"
#include<cmath>
#include<functional>
#include<vector>
using namespace std;

#define ENa 60
#define STEPSIZE 0.01

void Na_Euler( vector<double> &current, vector<double> &vol, int time_length )
{
	// Sodium channel simulation using Forward Euler solver

	// initialize
	vol.push_back(0);
	double m = 0, h = 1;

	// create the solver for ODEs
	function<double (double,double)> fm = f_m;
	Euler<double> euler_dm( STEPSIZE, fm );
	function<double (double,double)> fh = f_h;
	Euler<double> euler_dh( STEPSIZE, fh );
	function<double (double,double,double)> fv = f_v_Na;
	Euler<double> euler_dv( STEPSIZE, fv );

	// iterate through the time period and collect the current and voltage
	for( double t = 0, i = 0; t < time_length; t += STEPSIZE, ++i ){
		m = euler_dm.solve( vol[i], m );
		h = euler_dh.solve( vol[i], h );
		current.push_back( pow(m, 3)*h*(vol[i] - ENa) );
		vol.push_back( euler_dv.solve( m, h, vol[i]) );
	}
}