#ifndef RECOMBINATION_H
#define RECOMBINATION_H

#include<iostream>
#include<random>
#include<vector>
#include<cmath>
#include <chrono>
using namespace std;

// base class for recombination
template<class T> class recombination{

public:
    vector<T> parents;
	vector<vector<T> > tribe;
	vector<T> upper_bound;
	vector<T> lower_bound;
	vector<T> R;
	int type;
	int size_elites;
	recombination(vector<T>,vector<vector<T> >,vector<T>,vector<T>,vector<T>,int,int);
	vector<vector<T> > recombine();
};

// constructor
template<class T> recombination<T>::recombination(vector<T> parents,
									vector<vector<T> > tribe,vector<T> upper_bound,
									vector<T> lower_bound,vector<T> R,int type,
									int size_elites):parents(parents),tribe(tribe),
									upper_bound(upper_bound),lower_bound(lower_bound),
									R(R),type(type),size_elites(size_elites){}

// Elitism
template<class T> class elitism : public recombination<T>{

public:
	elitism(vector<T> parents,vector<vector<T> > tribe,vector<T> upper_bound,
			vector<T> lower_bound, vector<T> R, int size_elites, int type) :
			recombination<T>(parents,tribe,upper_bound,lower_bound,R,size_elites,type){}
	vector<vector<T> > run();
};

template<class T> vector<vector<T> > elitism<T>::run(){

    unsigned seed = static_cast<int> (chrono::system_clock::now().time_since_epoch().count());
	mt19937 mt(seed); uniform_real_distribution<T> u(0,1);
	vector<vector<T> > newpop(this->tribe.size());

	for ( size_t i = 0; i < newpop.size()-this->size_elites; ++i ){
		newpop[i].resize(this->tribe[0].size());
		vector<T> parent1(newpop[i].size()), parent2(newpop[i].size());
		int modnum = 2*(i % this->parents.size()/2);
		for( size_t j = 0; j < parent1.size(); ++j ){
			parent1[j] = this->tribe[this->parents[modnum]][j];
			parent2[j] = this->tribe[this->parents[modnum+1]][j];
		}

		// Find new parameters
		T plb, pub;
		for( size_t j = 0; j < newpop[0].size(); ++j ){
			if( parent1[j] < parent2[j] ){
				plb = parent1[j] - 0.25*abs(parent1[j]);
				pub = parent2[j] + 0.25*abs(parent2[j]);
			}
			else{
				plb = parent2[j] - 0.25*abs(parent2[j]);
				pub = parent1[j] + 0.25*abs(parent1[j]);
			}
			// Stay in bounds
			if( plb < this->lower_bound[j] )
				plb = this->lower_bound[j];
			if( pub > this->upper_bound[j] )
				pub = this->upper_bound[j];

			// Pick new parameter
			newpop[i][j] = plb + (pub-plb)*u(mt);
		}
	}

	// Re-insert elite members
	for( size_t i = this->tribe.size() - this->size_elites; i < this->tribe.size(); ++i ){
		newpop[i].resize(this->tribe[0].size());
		for( size_t j = 0; j < this->tribe[0].size(); ++j )
			newpop[i][j] = this->tribe[this->R[i-(this->tribe.size()-this->size_elites)]][j];
	}
	return newpop;
}

template<class T> vector<vector<T> > recombination<T>::recombine(){
	if( this->type == 1 ){
		elitism<T> elite(parents,tribe,upper_bound,lower_bound,R,type,size_elites);
		return elite.run();
	}
	else{
		cout << "Recombination error!" << endl;
	}
}

#endif
