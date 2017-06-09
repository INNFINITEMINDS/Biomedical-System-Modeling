#include "hw7.h"
#include<cmath>
#include<functional>
#include<vector>
using namespace std;

#define EK -80
#define ENa 60
#define g_Na 1.2
#define g_K 0.36
#define STEPSIZE 0.01

void HH_RK4( vector<double> &current_K, vector<double> &current_Na,
				 vector<double> &vol , int time_length )
{
	// complete Hodgkins-Huxley model
	// conductance
	// double g_Na = 1.2, g_K = 0.36, g_l = 0.003;
	vol.push_back(-60);
	double n, m, h;
	n = 0.01*(vol[0]+50)/(1-exp(-(vol[0]+50)/10)) / 
		(0.01*(vol[0]+50)/(1-exp(-(vol[0]+50)/10))+0.125*exp(-(vol[0]+60)/80));
	m = 0.1*(vol[0]+35)/(1-exp(-(vol[0]+35)/10)) /
	 (0.1*(vol[0]+35)/(1-exp(-(vol[0]+35)/10))+4*exp(-(vol[0]+60)/18));
	h = 0.07*exp(-(vol[0]+60)/20) /
		(0.07*exp(-(vol[0]+60)/20)+1/(1+exp(-(vol[0]+30)/10)));

	function<double (double,double)> fn = f_n;
	RK4<double> RK4_dn( STEPSIZE, fn );
	function<double (double,double)> fm = f_m;
	RK4<double> RK4_dm( STEPSIZE, fm );
	function<double (double,double)> fh = f_h;
	RK4<double> RK4_dh( STEPSIZE, fh );
	function<double (double,double,double,double)> fv = f_v_HH;
	RK4<double> RK4_dv( STEPSIZE, fv );
	for( double t = 0, i = 0; t < time_length; t += STEPSIZE, ++i ){
		n = RK4_dn.solve( vol[i], n );
		m = RK4_dm.solve( vol[i], m );
		h = RK4_dh.solve( vol[i], h );
		current_K.push_back( g_K*pow(n, 4)*(vol[i] - EK) );
		current_Na.push_back( g_Na*pow(m, 3)*h*(vol[i] - ENa) );
		vol.push_back( RK4_dv.solve( m, h, n, vol[i]) );
	}
}