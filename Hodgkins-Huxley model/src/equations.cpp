#include "hw7.h"
#include<cmath>

#define EK -80
#define ENa 60
#define El -50
#define g_Na 1.2
#define g_K 0.36
#define g_l 0.003

// version for double
double f_n( double v, double n )
{
	return 0.01*(v+50)*(1-n)/(1-exp(-(v+50)/10)) - 0.125*exp(-(v+60)/80)*n;
}

double f_m( double v, double m )
{
	return 0.1*(v+35)*(1-m)/(1-exp(-(v+35)/10)) - 4*exp(-(v+60)/18)*m;
}

double f_h( double v, double h )
{
	return 0.07*exp(-(v+60)/20)*(1-h) - h/(1+exp(-(v+30)/10));
}

double f_v_K( double n, double v )
{
	return -pow(n, 4)*(v - EK);
}

double f_v_Na( double m, double h, double v )
{
	return -pow(m, 3)*h*(v - ENa);
}

double f_v_both( double m, double h, double n, double v)
{
	return -pow(m, 3)*h*(v - ENa) - pow(n, 4)*(v - EK);
}

double f_v_HH( double m, double h, double n, double v)
{
	return (0.1 - g_Na*pow(m, 3)*h*(v - ENa) - g_K*pow(n, 4)*(v - EK) - g_l*(v - El))/0.01;
}

// version for float
float f_n( float v, float n )
{
	return 0.01*(v+50)*(1-n)/(1-exp(-(v+50)/10)) - 0.125*exp(-(v+60)/80)*n;
}

float f_m( float v, float m )
{
	return 0.1*(v+35)*(1-m)/(1-exp(-(v+35)/10)) - 4*exp(-(v+60)/18)*m;
}

float f_h( float v, float h )
{
	return 0.07*exp(-(v+60)/20)*(1-h) - h/(1+exp(-(v+30)/10));
}

float f_v_K( float n, float v )
{
	return -pow(n, 4)*(v - EK);
}

float f_v_Na( float m, float h, float v )
{
	return -pow(m, 3)*h*(v - ENa);
}

float f_v_both( float m, float h, float n, float v)
{
	return -pow(m, 3)*h*(v - ENa) + pow(n, 4)*(v - EK);
}

float f_v_HH( float m, float h, float n, float v)
{
	return (0.1 - g_Na*pow(m, 3)*h*(v - ENa) - g_K*pow(n, 4)*(v - EK) - g_l*(v - El))/0.01;
}