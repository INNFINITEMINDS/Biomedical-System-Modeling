#ifndef SELECTION_H
#define SELECTION_H

#include<iostream>
#include<random>
#include<vector>
#include<cmath>
#include <chrono>
using namespace std;

// base class for selection operator
template<class T> class selection{

public:
    vector<T> fitness;
	int size_parents;
	int type;
	selection(vector<T>,int,int);
	vector<T> select();
};

// constructor
template<class T> selection<T>::selection(vector<T> fitness,int size_parents,int type):
								fitness(fitness),size_parents(size_parents),type(type){}

// Russian Roulette
template<class T> class roulette : public selection<T>{
public:
	roulette(vector<T> fitness,int size_parents,int type) :
			selection<T>(fitness,size_parents,type){}
	vector<T> run();
};

template<class T> vector<T> roulette<T>::run(){

	unsigned seed = static_cast<int> (chrono::system_clock::now().time_since_epoch().count());
	mt19937 mt(seed); uniform_real_distribution<T> u(0,1);
	vector<T> parents(this->size_parents);

	// How long is our line
	T linelength = 0;
	for( auto f : this->fitness )
        linelength += f;

	for( size_t i = 0; i < parents.size(); ++i ){
		T pick = linelength*u(mt);
		T popsum = 0;
		for( size_t j = 0; j < this->fitness.size(); ++j ){
			popsum += this->fitness[j];
			if( pick < popsum ){
				parents[i] = j;
				break;
			}
		}
	}

	return parents;
}

// select
template<class T> vector<T> selection<T>::select(){
	if( this->type == 1 ){
		roulette<T> rou(fitness,size_parents,type);
		return rou.run();
	}
	else{
		cout << "Selection error!" << endl;
	}
}

#endif
