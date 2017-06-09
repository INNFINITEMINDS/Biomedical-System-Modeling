#include "hw7.h"
#include<cmath>
#include<functional>
#include<vector>
using namespace std;

#define EK -80
#define ENa 60
#define STEPSIZE 0.01

void Both_RK4( vector<double> &current_K, vector<double> &current_Na,
				 vector<double> &vol, int time_length )
{
	// both channels simulation using 4th order Runge-Kutta solver

	// initialize
	vol.push_back(0);
	double n = 0, m = 0, h = 1;

	// create the solver for ODEs
	function<double (double,double)> fn = f_n;
	RK4<double> RK4_dn( STEPSIZE, fn );
	function<double (double,double)> fm = f_m;
	RK4<double> RK4_dm( STEPSIZE, fm );
	function<double (double,double)> fh = f_h;
	RK4<double> RK4_dh( STEPSIZE, fh );
	function<double (double,double,double,double)> fv = f_v_both;
	RK4<double> RK4_dv( STEPSIZE, fv );

	// iterate through the time period and collect the current and voltage
	for( double t = 0, i = 0; t < time_length; t += STEPSIZE, ++i ){
		n = RK4_dn.solve( vol[i], n );
		m = RK4_dm.solve( vol[i], m );
		h = RK4_dh.solve( vol[i], h );
		current_K.push_back( pow(n, 4)*(vol[i] - EK) );
		current_Na.push_back( pow(m, 3)*h*(vol[i] - ENa) );
		vol.push_back( RK4_dv.solve( m, h, n, vol[i]) );
	}
}