//============================================================================
//
// 	Bidirectional_work-distribution_builder.cpp
// 	Copyright 2014  © Mostafa Nategholeslam
//
//	This file is part of Bidirectional_work-distribution_builder.
//
//    Bidirectional_work-distribution_builder is free software: you can redistribute it and/or modify
//    it under the terms of the GNU General Public License as published by
//    the Free Software Foundation, either version 3 of the License, or
//    (at your option) any later version.
//
//    Bidirectional_work-distribution_builder is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU General Public License for more details.
//
//    You should have received a copy of the GNU General Public License
//    along with Foobar.  If not, see <http://www.gnu.org/licenses/>.
//
 //============================================================================


#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include <vector>
using namespace std;
#include <iomanip>
#include <locale>
#include <sstream>
#include <sys/stat.h>
#include <time.h>
//#include <sys/types.h>
#include <cstring>
#include "distributionsCalculation.h"


int main(int argc, char* argv[]) {

//  gtk_init(&argc, &argv);
// Win();
//     gtk_main();

	double xInit, xFin, dt=2.0000;
	int numberOfBins, numberOfWorkBins, numberOfVelocityBins;
	cout<<endl<<"Enter the initial value of the range of reaction coordinate:  ";
	cin>>xInit;
	cout<<"Enter the final value of the range of reaction coordinate:  ";
	cin>>xFin;
	cout<<"How many bins do you want along this range of reaction coordinate?  ";
	cin>>numberOfBins;
	cout <<"Enter the length of the time-step used for this simulation (in femtoseconds): ";
	cin >> dt;
	dt *= 1.00000e-6;
	cout <<"Enter the desired number of work bins for each position bin: ";
		cin >> numberOfWorkBins;
	cout <<"Enter the desired number of velocity bins for each position bin: ";
		cin >> numberOfVelocityBins;


	double binSize=(xFin-xInit)/(double)numberOfBins;   //ATTENTION: This can be negative.
	cout <<"\n Bin size = "<<binSize<<endl<<endl;

	vector <xBin> forward, reverse;
	xBinVectorCreator (numberOfBins, numberOfWorkBins, numberOfVelocityBins, xInit, binSize, forward, dt);
	xBinVectorCreator (numberOfBins, numberOfWorkBins, numberOfVelocityBins, xInit, binSize, reverse, dt);


	time_t seconds;
	seconds = time (NULL);
	cout.precision(18);

	stages (forward, reverse);

	seconds = time (NULL) - seconds;
	cout<<endl<<endl<<"Done!"<<endl<<"And it took "<<seconds<<" seconds to do it!"<<endl<<endl;


	return 0;
}
