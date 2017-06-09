#include "hw7.h"
#include<iostream>
#include<fstream>
#include<cmath>
#include<functional>
#include<vector>
using namespace std;

int main()
{
	cout << "Type the number to get corresponding outputs:" << endl;
	cout << "1 for Potassium channel simulation using Forward Euler solver;" << endl;
	cout << "2 for Sodium channel simulation using Forward Euler solver;" << endl;
	cout << "3 for both channels simulation using 4th order Runge-Kutta solver;" << endl;
	cout << "4 for both channels simulation using 5th order Runge-Kutta solver with adaptive time step;" << endl;
	cout << "5 for 5th order Runge-Kutta soler with 3rd order cubic interpolation;" << endl;
	cout << "6 for complete Hodgkins-Huxley model (using 4th order Runge-Kutta solver)." << endl;
	cout << "> ";
	int ans;
	cin >> ans;
	cout << "Type the time length (integer) of simulation: " << endl << "> ";
	int time_length;
	cin >> time_length;
	vector<double> current, vol, current_K, current_Na;
	switch(ans){
		case 1:
			K_Euler( current, vol, time_length );
			break;
		case 2:
			Na_Euler( current, vol, time_length );
			break;
		case 3:
			Both_RK4( current_K, current_Na, vol, time_length );
			break;
		case 4:
			Both_RK5( current_K, current_Na, vol, time_length );
			break;
		case 5:
			RK5_interp3( current_K, current_Na, vol, time_length );
			break;
		case 6:
			HH_RK4( current_K, current_Na, vol, time_length );
			break;
	}

	// output the data (voltage and current)
	if( ans < 3 ){
		ofstream file;
    	file.open("output.txt");
    	for( size_t i = 0; i < current.size(); ++i ){
        	file << vol[i] << "\t" << current[i] << endl;
    	}
    	file.close();
    	cout << "Data output completed in 'output.txt' in the same directory." << endl;
    	cout << "1st column: voltage; 2nd column: current." << endl;
	}
	else{
		ofstream file;
    	file.open("output.txt");
    	for( size_t i = 0; i < vol.size(); ++i ){
        	file << vol[i] << "\t" << current_K[i] << "\t" << current_Na[i] << endl;
    	}
    	file.close();
    	cout << "Data output completed in 'output.txt' in the same directory." << endl;
    	cout << "1st column: voltage; 2nd column: current of Potassium; 3rd column: current of Sodium." << endl;
	}
	return 0;
}