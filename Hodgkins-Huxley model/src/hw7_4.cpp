#include "hw7.h"
#include<cmath>
#include<functional>
#include<vector>
using namespace std;

#define EK -80
#define ENa 60
#define STEPSIZE 0.01

void Both_RK5( vector<double> &current_K, vector<double> &current_Na,
				 vector<double> &vol, int time_length )
{
	// both channels simulation using 5th order Runge-Kutta solver
	// initialize
	vol.push_back(0);
	double n = 0, m = 0, h = 1;

	// create the solver for ODEs
	function<double (double,double)> fn = f_n;
	RK5<double> RK5_dn( STEPSIZE, fn );
	function<double (double,double)> fm = f_m;
	RK5<double> RK5_dm( STEPSIZE, fm );
	function<double (double,double)> fh = f_h;
	RK5<double> RK5_dh( STEPSIZE, fh );
	function<double (double,double,double,double)> fv = f_v_both;
	RK5<double> RK5_dv( STEPSIZE, fv );

	// iterate through the time period and collect the current and voltage and time step
    vector<double> time_step;
	for( double t = RK5_dv.stepSize, i = 0; t < time_length; t += RK5_dv.stepSize, ++i ){
		n = RK5_dn.solve( vol[i], n );
		m = RK5_dm.solve( vol[i], m );
		h = RK5_dh.solve( vol[i], h );
		current_K.push_back( pow(n, 4)*(vol[i] - EK) );
		current_Na.push_back( pow(m, 3)*h*(vol[i] - ENa) );
		vol.push_back( RK5_dv.solve( m, h, n, vol[i]) );
        double minStep = min(min(min(RK5_dn.stepSize,RK5_dm.stepSize),RK5_dh.stepSize),RK5_dv.stepSize);
        RK5_dn.stepSize = minStep;
        RK5_dm.stepSize = minStep;
        RK5_dh.stepSize = minStep;
        RK5_dv.stepSize = minStep;
        time_step.push_back( minStep );
	}
}