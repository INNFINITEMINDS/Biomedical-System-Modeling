#ifndef GA_H
#define GA_H

#include<iostream>
#include<random>
#include<vector>
#include<cmath>
#include <chrono>
#include "selection.h"
#include "recombination.h"
#include "mutation.h"
using namespace std;

float minfun(vector<float> x){

	return pow(100*(x[0]-0.11),2)
	       +pow(10*(x[1]-7.3),2)
	       +pow(x[2]-89.4,2);
}

double minfun(vector<double> x){

	return pow(100*(x[0]-0.11),2)
	       +pow(10*(x[1]-7.3),2)
	       +pow(x[2]-89.4,2);
}

template<class T> class GA{

public:
    vector<T> upper_bound;
	vector<T> lower_bound;
	int population;
	int size_parent;
	int generation;
	GA(vector<T>,vector<T>,int,int,int);
	vector< vector<T> > tribe;
	vector<T> run();
};

// constructor
template<class T> GA<T>::GA(vector<T> upper_bound, vector<T> lower_bound,
							int population, int size_parent, int generation):
							upper_bound(upper_bound),lower_bound(lower_bound),
							population(population),size_parent(size_parent),
							generation(generation){
	// initialize the population
	tribe.resize(population);
	unsigned seed = static_cast<int> (chrono::system_clock::now().time_since_epoch().count());
	mt19937 mt(seed); uniform_real_distribution<T> u(0,1);
	for( auto &p : tribe ){
		// Allocate space for the parameters
		p.resize(upper_bound.size());
		for( size_t i = 0; i < upper_bound.size(); ++i )
			p[i] = lower_bound[i] + (upper_bound[i]-lower_bound[i])*u(mt);
	}

	// Output initial population
	cout << "Initial population:" << endl;
	for( auto p : tribe ) {
		for ( auto q : p ) cout << q << " ";
		cout << endl;
	}
	cout << endl;
}

template<class T> vector<T> GA<T>::run(){

	mt19937 mt(0); uniform_real_distribution<T> u(0,1);
	int best;
	for( int age = 0; age < generation; ++age ){

		cout << "Generation: " << age + 1 << endl;

		// compute the objective values
		vector<T> vals(this->population);
		for( int i = 0; i < this->population; ++i )
			vals[i] = minfun(this->tribe[i]);

		vector <int> nums(vals.size());
		for( size_t i = 0; i < vals.size(); ++i )
			nums[i] = i;

		// Vector to hold ranks
		vector<T> R(vals.size());

		// Find max
		T max = 0;
		for( auto val : vals )
			if( val > max )
				max = val;

		// sort rank in ascending order of objective values
		for( size_t i = 0; i < vals.size()-1; ++i ){
			T min = max; int nmin;
			for(auto &n : nums){
				if( vals[n] < min ){
					min = vals[n];
					nmin = n;
				}
			}

			//  Delete nmin from nums
			for(size_t j = 0; j < nums.size(); ++j ){
				if( nums[j] == nmin ){
					nums.erase(nums.begin()+j);
					break;
				}
			}

			// Set R[i] to nmin
			R[i] = nmin;
		}

		// Set the last value
		R[R.size()-1] = nums[0];

		// Compute Fitness
		double S = 2;
		vector<T> F(vals.size());
		int iter = 0;
		for(auto r : R){
			F[r] = 2 - S + 2*(S-1)*((T)(vals.size()-iter-1)/(T)(vals.size()-1));
			++iter;
		}

		// Selection
		selection<T> s(F,size_parent,1);
		vector<T> parents(size_parent);
		parents = s.select();

		// Recombination
		recombination<T> rec(parents,tribe,upper_bound,lower_bound,R,1,5);
		vector<vector<T> > newpop(population);
		newpop = rec.recombine();

		// Mutation
		mutation<T> mut(5,upper_bound,lower_bound,newpop,1);
		newpop = mut.mutate();

		for( size_t i = 0; i < tribe.size(); ++i ){
			for( size_t j = 0; j < upper_bound.size(); ++j )
				tribe[i][j] = newpop[i][j];
		}
		best = R[0];
	}

	return tribe[best];
}

#endif