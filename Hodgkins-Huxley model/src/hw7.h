#ifndef HW7_H
#define HW7_H
#include<iostream>
#include<fstream>
#include<cmath>
#include<functional>
#include<vector>
using namespace std;

#define eps 1.0e-3
#define STEPSIZE 0.01

double f_n( double v, double n );
double f_m( double v, double m );
double f_h( double v, double h );
double f_v_K( double n, double v );
double f_v_Na( double m, double h, double v );
double f_v_both( double m, double h, double n, double v);
double f_v_HH( double m, double h, double n, double v);

void K_Euler( vector<double> &current, vector<double> &vol, int time_length );
void Na_Euler( vector<double> &current, vector<double> &vol, int time_length );
void Both_RK4( vector<double> &current_K, vector<double> &current_Na,
				 vector<double> &vol, int time_length );
void Both_RK5( vector<double> &current_K, vector<double> &current_Na,
				 vector<double> &vol, int time_length );
void RK5_interp3( vector<double> &current_K, vector<double> &current_Na,
				 vector<double> &vol, int time_length );
void HH_RK4( vector<double> &current_K, vector<double> &current_Na,
				 vector<double> &vol , int time_length );

// base class
template<class T> class Solver{
	public:
    	T stepSize;
		function<T (T,T)> f;
		function<T (T,T,T)> f_Na;
		function<T (T,T,T,T)> f_both;
		Solver( T stepSize, function<T (T,T)> f );
		Solver( T stepSize, function<T (T,T,T)> f_Na );
		Solver( T stepSize, function<T (T,T,T,T)> f_both );
		virtual T solve( T t, T y ) = 0;
		void setStepSize( T stepsize ){ stepSize = stepSize; }
};

// constructor
template<class T> Solver<T>::Solver( T stepSize, function<T (T,T)> f ): stepSize(stepSize), f(f){};
template<class T> Solver<T>::Solver( T stepSize, function<T (T,T,T)> f_Na ): stepSize(stepSize), f_Na(f_Na){};
template<class T> Solver<T>::Solver( T stepSize, function<T (T,T,T,T)> f_both ): stepSize(stepSize), f_both(f_both){};

// Forward Euler
template<class T> class Euler : public Solver<T>{
	public:
		Euler( T stepSize, function<T (T,T)> f ): Solver<T>( stepSize, f ){}
		Euler( T stepSize, function<T (T,T,T)> f_Na ): Solver<T>( stepSize, f_Na ){}
		Euler( T stepSize, function<T (T,T,T,T)> f_both ): Solver<T>( stepSize, f_both ){}
		T solve( T t, T y );
		T solve( T t, T h, T y );
		T solve( T t, T h, T n, T y );
};

template<class T> T Euler<T>::solve( T t, T y ){
	return y + this->stepSize*this->f(t,y);
}

template<class T> T Euler<T>::solve( T t, T h, T y ){
	return y + this->stepSize*this->f_Na(t,h,y);
}

template<class T> T Euler<T>::solve( T t, T h, T n, T y ){
	return y + this->stepSize*this->f_both(t,h,n,y);
}

// 4th Runge Kutta
template<class T> class RK4 : public Solver<T>{
	public:
		RK4( T stepSize, function<T (T,T)> f ): Solver<T>( stepSize, f ){}
		RK4( T stepSize, function<T (T,T,T)> f_Na ): Solver<T>( stepSize, f_Na ){}
		RK4( T stepSize, function<T (T,T,T,T)> f_both ): Solver<T>( stepSize, f_both){}
		T solve( T t, T y );
		T solve( T t, T h, T y );
		T solve( T t, T h, T n, T y );
};

template<class T> T RK4<T>::solve( T t, T y ){
    T stepSize = this->stepSize;
	T k1 = this->f(t,y);
	T k2 = this->f(t+stepSize/2,y+k1*stepSize/2);
	T k3 = this->f(t+stepSize/2,y+k2*stepSize/2);
	T k4 = this->f(t+stepSize,y+stepSize*k3);
	return y + stepSize*(k1+2*k2+2*k3+k4)/6;
}

template<class T> T RK4<T>::solve( T t, T h, T y ){
    T stepSize = this->stepSize;
	T k1 = this->f_Na(t,h,y);
	T k2 = this->f_Na(t+stepSize/2,h+stepSize/2,y+k1*stepSize/2);
	T k3 = this->f_Na(t+stepSize/2,h+stepSize/2,y+k2*stepSize/2);
	T k4 = this->f_Na(t+stepSize,h+stepSize,y+k3*stepSize);
	return y + stepSize*(k1+2*k2+2*k3+k4)/6;
}

template<class T> T RK4<T>::solve( T t, T h, T n, T y ){
    T stepSize = this->stepSize;
	T k1 = this->f_both(t,h,n,y);
	T k2 = this->f_both(t+stepSize/2,h+stepSize/2,n+stepSize/2,y+k1*stepSize/2);
	T k3 = this->f_both(t+stepSize/2,h+stepSize/2,n+stepSize/2,y+k2*stepSize/2);
	T k4 = this->f_both(t+stepSize,h+stepSize,n+stepSize/2,y+k3*stepSize);
	return y + stepSize*(k1+2*k2+2*k3+k4)/6;
}


// 5th Runge Kutta
template<class T> class RK5 : public Solver<T>{
	public:
		RK5( T stepSize, function<T (T,T)> f ): Solver<T>( stepSize, f ){}
		RK5( T stepSize, function<T (T,T,T)> f_Na ): Solver<T>( stepSize, f_Na ){}
		RK5( T stepSize, function<T (T,T,T,T)> f_both ): Solver<T>( stepSize, f_both){}
		T solve( T t, T y );
		T solve( T t, T h, T y );
		T solve( T t, T h, T n, T y );
};

template<class T> T RK5<T>::solve( T t, T y ){
	T stepSize = this->stepSize;
	T k1 = this->f(t,y);
	T k2 = this->f(t+stepSize/5,y+k1/5);
	T k3 = this->f(t+stepSize*3/10,y+k1*3/40+k2*9/40);
	T k4 = this->f(t+stepSize*3/5,y+k1*3/10-k2*9/10+k3*6/5);
	T k5 = this->f(t+stepSize,y-k1*11/54+k2*5/2-k3*70/27+k4*35/27);
	T k6 = this->f(t+stepSize*7/8,y+k1*1631/55296+k2*175/512+k3*575/13824+k4*44275/110592-k5*253/4096);
	double res = y + stepSize*(37*k1/378+250*k3/621+125*k4/594+512*k6/1771);

	// adaptive step size
	// compute the embedded 4th order formula for error calculation
	double res_s = y + stepSize*(2825*k1/27648+18575*k3/48384+13525*k4/55296-277*k5/14336+k6/4);
	double error = pow(abs(res-res_s), 0.2);

	double new_step = max(0.9 * stepSize * pow(eps,0.2) / error, 0.2*STEPSIZE);
    new_step = min(new_step, 5*STEPSIZE);

	this->stepSize = new_step;

	return res;
}

template<class T> T RK5<T>::solve( T t, T h, T y ){
	T stepSize = this->stepSize;
	T k1 = this->f_Na(t,h,y);
	T k2 = this->f_Na(t+stepSize/5,h+stepSize/5,y+k1/5);
	T k3 = this->f_Na(t+stepSize*3/10,h+stepSize*3/10,y+k1*3/40+k2*9/40);
	T k4 = this->f_Na(t+stepSize*3/5,h+stepSize*3/5,y+k1*3/10-k2*9/10+k3*6/5);
	T k5 = this->f_Na(t+stepSize,h+stepSize,y-k1*11/54+k2*5/2-k3*70/27+k4*35/27);
	T k6 = this->f_Na(t+stepSize*7/8,h+stepSize*7/8,y+k1*1631/55296+k2*175/512+k3*575/13824+k4*44275/110592-k5*253/4096);
	double res = y + stepSize*(37*k1/378+250*k3/621+125*k4/594+512*k6/1771);

	// adaptive step size
	// compute the embedded 4th order formula for error calculation
	double res_s = y + stepSize*(2825*k1/27648+18575*k3/48384+13525*k4/55296-277*k5/14336+k6/4);
    double error = pow(abs(res-res_s), 0.2);

	double new_step = max(0.9 * stepSize * pow(eps,0.2) / error, 0.2*STEPSIZE);
    new_step = min(new_step, 5*STEPSIZE);

	this->stepSize = new_step;

	return res;
}

template<class T> T RK5<T>::solve( T t, T h, T n, T y ){
	T stepSize = this->stepSize;
	T k1 = this->f_both(t,h,n,y);
	T k2 = this->f_both(t+stepSize/5,h+stepSize/5,n+stepSize/5,y+k1/5);
	T k3 = this->f_both(t+stepSize*3/10,h+stepSize*3/10,n+stepSize*3/10,y+k1*3/40+k2*9/40);
	T k4 = this->f_both(t+stepSize*3/5,h+stepSize*3/5,n+stepSize*3/5,y+k1*3/10-k2*9/10+k3*6/5);
	T k5 = this->f_both(t+stepSize,h+stepSize,n+stepSize,y-k1*11/54+k2*5/2-k3*70/27+k4*35/27);
	T k6 = this->f_both(t+stepSize*7/8,h+stepSize*7/8,n+stepSize*7/8,y+k1*1631/55296+k2*175/512+k3*575/13824+k4*44275/110592-k5*253/4096);
	double res = y + stepSize*(37*k1/378+250*k3/621+125*k4/594+512*k6/1771);

	// adaptive step size
	// compute the embedded 4th order formula for error calculation
	double res_s = y + stepSize*(2825*k1/27648+18575*k3/48384+13525*k4/55296-277*k5/14336+k6/4);
	double error = pow(abs(res-res_s), 0.2);

	double new_step = max(0.9 * stepSize * pow(eps,0.2) / error, 0.2*STEPSIZE);
	new_step = min(new_step, 5*STEPSIZE);

	this->stepSize = new_step;

	return res;
}
#endif