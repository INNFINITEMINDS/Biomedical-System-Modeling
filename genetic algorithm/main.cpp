#include<iostream>
#include<random>
#include<vector>
#include<cmath>
#include "GA.h"
#include "selection.h"
#include "recombination.h"
#include "mutation.h"
using namespace std;

int main(){
    // set parameters
    vector<float> ub={1,10,100}, lb={-1,-10,-100};
	int population = 50, generation = 1000, size_parent = 2;

	// run GA
    GA<float> optimizer(ub, lb, population, size_parent, generation);
    vector<float> result;
    result = optimizer.run();

    cout << "Result: " << endl;
    for( auto val : result )
        cout << val << " ";
    cout << endl;
    return 0;
}