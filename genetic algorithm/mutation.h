#ifndef MUTATION_H
#define MUTATION_H

#include<iostream>
#include<random>
#include<vector>
#include<cmath>
#include <chrono>
using namespace std;

// base class for mutation
template<class T> class mutation{

public:
    T mrate;
	T msize = 0.025;
	int size_elites;
	vector<T> upper_bound;
	vector<T> lower_bound;
	vector<vector<T> > newpop;
	int type;
	mutation(int,vector<T>,vector<T>,vector<vector<T> >,int);
	vector<vector<T> > mutate();
};

// constructor
template<class T> mutation<T>::mutation(int size_elites,vector<T> upper_bound,
	vector<T> lower_bound,vector<vector<T> > newpop,int type):size_elites(size_elites),
	upper_bound(upper_bound),lower_bound(lower_bound),newpop(newpop),type(type){
	mrate = 1/(T)newpop.size();
}

// One-point mutation
template<class T> class onepoint : public mutation<T>{

public:
	onepoint(int size_elites,vector<T> upper_bound,vector<T> lower_bound, vector<vector<T> > newpop, int type):
			mutation<T>(size_elites,upper_bound,lower_bound,newpop,type){}
	vector<vector<T> > run();
};

template<class T> vector<vector<T> > onepoint<T>::run(){

	unsigned seed = static_cast<int> (chrono::system_clock::now().time_since_epoch().count());
	mt19937 mt(seed); uniform_real_distribution<T> u(0,1);
	for( size_t i = 0; i < this->newpop.size() - this->size_elites; ++i ){
		// for each parameter
		for( size_t j = 0; j < this->upper_bound.size(); ++j ){
			T asize = this->msize*(this->upper_bound[j]-this->lower_bound[j])*(u(mt)-0.5);
			if(u(mt) < this->mrate){
				this->newpop[i][j] += asize;
				// Make sure were staying in bounds
				if( this->newpop[i][j] < this->lower_bound[j] )
					this->newpop[i][j] = this->lower_bound[j];
				if( this->newpop[i][j] > this->upper_bound[j] )
					this->newpop[i][j] = this->upper_bound[j];
			}
		}
	}

	return this->newpop;
}

template<class T> vector<vector<T> > mutation<T>::mutate(){

	if( this->type == 1){
		onepoint<T> one(size_elites,upper_bound,lower_bound,newpop,type);
		return one.run();
	}
	else{
		cout << "Mutation error!" << endl;
	}
}

#endif
