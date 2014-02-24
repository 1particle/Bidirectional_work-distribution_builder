//============================================================================
//
// 	distributionsCalculation.cpp
// 	Copyright 2014  © Mostafa Nategholeslam
//
//    This file is part of Bidirectional_work-distribution_builder.
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
#include <time.h>
#include <iomanip>
#include <locale>
#include <sstream>
#include <sys/stat.h>
//#include <sys/types.h>
#include <cstring>
#include "distributionsCalculation.h"


double maximumAbs (double a, double b)	{
	if (fabs(a) > fabs(b))
			return fabs(a);
	else
		return fabs(b);
}

double minimumAbs (double a, double b)	{
	if (fabs(a) < fabs(b))
			return fabs(a);
	else
		return fabs(b);
}




long file_size ()  {
  long begin,end;
  ifstream input ("FR.txt");
  begin = input.tellg();
  input.seekg (0, ios::end);
  end = input.tellg();
  input.close();
  return end-begin;
}


void velocityBoundariesDeterminator(vector<xBin> & bin, double X1, double X2)
{
  // Positive value for variable 'binSize' means forward direction is positive (from smaller to larger values of x) and negative means the reverse.
      double dX = X2 - X1;
      int i= (int) floor ((X1 - bin[0].xi)/bin[0].binSize);

      if (i>=0 && i<(int)bin.size())
      if ( (bin[i].binSize > 0 && bin[i].xi < X1 && X1  < bin[i].xf &&   bin[i].xi < X2  &&  X2  < bin[i].xf )
            || (bin[i].binSize < 0 && bin[i].xi > X1 && X1  > bin[i].xf &&   bin[i].xi > X2  &&  X2  > bin[i].xf ))
            {
                    if(bin[i].vMax < dX/bin[i].dt)  bin[i].vMax= dX/bin[i].dt;
                    if(bin[i].vMin > dX/bin[i].dt)  bin[i].vMin= dX/bin[i].dt;
                          }
}


void workBoundariesDeterminator(vector<xBin> & bin, double X1, double X2, double work)
{
  // Positive value for variable 'binSize' means forward direction is positive (from smaller to larger values of x) and negative means the reverse.
      double dX = X2 - X1;
      int i= (int) floor ((X1 - bin[0].xi)/bin[0].binSize);

      if (i>=0 && i<(int)bin.size())
      if ( (bin[i].binSize > 0 && bin[i].xi < X1 && X1  < bin[i].xf &&   bin[i].xi < X2  &&  X2  < bin[i].xf )
            || (bin[i].binSize < 0 && bin[i].xi > X1 && X1  > bin[i].xf &&   bin[i].xi > X2  &&  X2  > bin[i].xf ))
            {
    	  	  	  	double scaledWork = (work*fabs(bin[i].binSize / dX ));
                    if(bin[i].wMax < scaledWork)  bin[i].wMax= scaledWork ;
                    if(bin[i].wMin > scaledWork)  bin[i].wMin= scaledWork ;
                          }
}


void xBinPrepare (vector<xBin> & b1, vector<xBin> & b2)
{
for (int i=0; i < (int)b1.size() ; i++)             {
	b1[i].dv = (b1[i].vMax-b1[i].vMin)/(double)b1[i].numberOfVelocityBins;
	b2[i].dv = (b2[i].vMax-b2[i].vMin)/(double)b2[i].numberOfVelocityBins;
	b1[i].velocityVectorGenerator();
	b2[i].velocityVectorGenerator();

	b1[i].dw = (b1[i].wMax- b1[i].wMin)/(double)b1[i].numberOfWorkBins;
	b2[i].dw = (b2[i].wMax- b2[i].wMin)/(double)b2[i].numberOfWorkBins;

	b1[i].workVectorGenerator();
	b2[i].workVectorGenerator();
  }
}


void readInput_BoundaryFinder ( vector<xBin>  & forward,  vector<xBin> & reverse)              {

  long fileSize = file_size ();
  long current=1;
  int percent = 0;
  cout <<endl<<"Stage 1 --> Reading from FR.txt to determine minimum and maximum \n work and velocity values in each bin ....\n";
  string line;
  char* st;
  ifstream FR;
  FR.open("FR.txt");
  FR.precision(15);
  while ( ! FR.eof() )
          {
                  getline (FR,line);
                  current += (int)line.size()+1;  //Need to add 1 for the newline character which is not included in the line character. (?????)
                  if ( ((current*100)/fileSize) > percent)       {
                      percent ++;
                      cout <<"\r"<<percent<<"% read..... "<< std::flush;
                      }
                  double x1 , x2, F,  xx1, xx2, FF;
                  x1= strtod(line.c_str(),&st);
                  x2= strtod(st,&st);
                  F= strtod(st,&st);
                  xx1= strtod(st,&st);
                  xx2= strtod(st,&st);
                  FF= strtod(st,&st);

                  double X1 = x1-xx1;
                  double X2= x2-xx2;

                  double dX= X2-X1;
                  if (dX*forward[0].binSize > 0)	{
                	  velocityBoundariesDeterminator(forward, X1, X2);
                	  workBoundariesDeterminator(forward, X1, X2, F*(x2-x1)+FF*(xx2-xx1));
                  }
                  else	{
                	  velocityBoundariesDeterminator(reverse, X1, X2);
                	  workBoundariesDeterminator(reverse, X1, X2, F*(x2-x1)+FF*(xx2-xx1));
                  }
          }

  FR.close();

}


void xBinVectorCreator (int numberOfBins, int numberOfWorkBins, int numberOfVelocityBins, double xInit, double binSize, vector<xBin> & bins, double dt)  {
	for (int a=0; a< numberOfBins; a++)		{
		xBin *b = new xBin(xInit+(double)a*binSize, xInit+(double)(a+1)*binSize, dt, numberOfVelocityBins, numberOfWorkBins);
		bins.push_back(*b);
	}

}


void velocityCounter(vector<xBin> & bin, double X1, double X2, double work)    {
	// Positive value for variable 'binSize' means forward direction is positive (from smaller to larger values of x)
	//	and negative means the reverse.

	double dX = X2 - X1;
	double v = dX/bin[0].dt;
    int i= (int) floor ((X1 - bin[0].xi)/bin[0].binSize);

      if (i>=0 && i<(int)bin.size())
      if ( (bin[i].binSize > 0 && bin[i].xi < X1 && X1  < bin[i].xf &&   bin[i].xi < X2  &&  X2  < bin[i].xf )
            || (bin[i].binSize < 0 && bin[i].xi > X1 && X1  > bin[i].xf &&   bin[i].xi > X2  &&  X2  > bin[i].xf ))
            {
    	  int j= (int) floor ((v - bin[i].vMin)/bin[i].dv);
    	  if (j>=0 && j<(int)bin[i].V.size())	{
        	  bin[i].V[j].count ++;
        	  bin[i].V[j].work += work;
        	  bin[i].V[j].weight += fabs( dX / bin[i].binSize );
    	  }
                          }
}


void workCounter(vector<xBin> & bin, double X1, double X2, double work)    {
	// Positive value for variable 'binSize' means forward direction is positive (from smaller to larger values of x)
	//	and negative means the reverse.


	double dX = X2 - X1;

      int i= (int) floor ((X1 - bin[0].xi)/bin[0].binSize);

      if (i>=0 && i<(int)bin.size())
      if ( (bin[i].binSize > 0 && bin[i].xi < X1 && X1  < bin[i].xf &&   bin[i].xi < X2  &&  X2  < bin[i].xf )
            || (bin[i].binSize < 0 && bin[i].xi > X1 && X1  > bin[i].xf &&   bin[i].xi > X2  &&  X2  > bin[i].xf ))
            {
    	  	double scaledWork = ((work) * fabs(bin[i].binSize / dX) );
    	  	int j= (int) floor ((scaledWork - bin[i].wMin)/bin[i].dw);
    	  	//double v = dX/bin[i].dt;
    	  	if (j>=0 && j<(int)bin[i].W.size() /*&& bin[i].vMin<= v && v <= bin[i].vMax*/)	{
        	  bin[i].W[j].count ++;
        	  bin[i].W[j].weight += fabs( dX / bin[i].binSize );
        	  bin[i].wMean += scaledWork;
        	  bin[i].tsCount ++;
    	  }
                          }
}


void readInput_sampler ( vector<xBin>  & forward,  vector<xBin> & reverse, int stage)    {

  long fileSize = file_size ();
  long current=1;
  int percent = 0;
  string line;
  char* st;
  ifstream FR;
  FR.open("FR.txt");
  FR.precision(15);
  while ( ! FR.eof() )
          {
	  getline (FR,line);
	  current += (int)line.size()+1;  //Need to add 1 for the newline character which is not included in the line character. (?????)
	  if ( ((current*100)/fileSize) > percent)       {
		  percent ++;
		  cout <<"\r"<<"Stage "<<stage<<" -->  "<<percent<<"% read..... "<< std::flush;
		  }
	  double x1 , x2, F,  xx1, xx2, FF;
	  x1= strtod(line.c_str(),&st);
	  x2= strtod(st,&st);
	  F= strtod(st,&st);
	  xx1= strtod(st,&st);
	  xx2= strtod(st,&st);
	  FF= strtod(st,&st);

	  double X1 = x1-xx1;
	  double X2= x2-xx2;

	  double dX= (X2-X1);
	  if (dX*forward[0].binSize >= 0)	{
		  velocityCounter(forward, X1, X2, F*(x2-x1)+ FF*(xx2-xx1));
		  workCounter(forward, X1, X2, F*(x2-x1)+ FF*(xx2-xx1));
	  }
	  else	{
		  velocityCounter(reverse, X1, X2, F*(x2-x1)+ FF*(xx2-xx1));
		  workCounter(reverse, X1, X2, F*(x2-x1)+ FF*(xx2-xx1));
	  }

          }

  FR.close();

}


void outputVelocityDistributions(vector<xBin> forward, vector <xBin> reverse)	{

	mkdir("velocity_distributions",S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);

	for(int i=0; i<(int)forward.size(); i++)	{
		ofstream results;
		  results.precision(18);
		  ostringstream x1, x2;
		  x1 << forward[i].xi;
		  x2 << forward[i].xf;
		  string filename = "velocity_distributions/Bin_"+ x1.str() + "_to_"+ x2.str() + "_Angstrom.txt";
		  results.open(filename.c_str());

		  double wForward =0, wReverse=0;
		  results<<"vForward   wForward   vForwardCount   vReverse   wReverse   vReverseCount \n";
		  for (int j=0; j< (int) forward[i].V.size(); j++)		{
			  if (forward[i].V[j].weight > 0)
				  wForward= forward[i].V[j].work/forward[i].V[j].weight;
			  else wForward = 0;
			  results <<(forward[i].V[j].vMin + forward[i].V[j].vMax)/2<<"  "<<wForward <<"   " <<forward[i].V[j].count<<"  ";

			  if (reverse[i].V[j].weight > 0)
				  wReverse= reverse[i].V[j].work/reverse[i].V[j].weight;
			  else  wReverse = 0;
			  results <<(reverse[i].V[j].vMin + reverse[i].V[j].vMax)/2<<"  "<<wReverse<<"   "<<reverse[i].V[j].count<<endl;
		  }
		  results.close();
	}
}


void outputWorkDistributions(vector<xBin> forward, vector <xBin> reverse)	{

	mkdir("work_distributions",S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
	for(int i=0; i<(int)forward.size(); i++)	{
		ofstream results;
		results.precision(18);
		ostringstream x1, x2;
		x1 << forward[i].xi;
		x2 << forward[i].xf;
		string filename = "work_distributions/Bin_"+ x1.str() + "_to_"+ x2.str() + "_Angstrom.txt";
		results.open(filename.c_str());

		results<<"wForward   wForwardCount  wReverse   wReverseCount\n";
		for (int j=0; j< (int)forward[i].W.size() ; j++)		{
		results <<(forward[i].W[j].wMin + forward[i].W[j].wMax)/2<<"  "<<forward[i].W[j].count <<"   ";
		results <<(reverse[i].W[j].wMin + reverse[i].W[j].wMax)/2<<"  "<<reverse[i].W[j].count <<endl;
		}
/////////////////////////////////////////////////////////////////////////////

//		for (int j=0; j< (int) reverse[i].W.size(); j++)	{
//			if (reverse[i].W[j].count >= 100) { minReverse = j; break;}		}
//
//		for (int j=(int)(reverse[i].W.size()) - 1 ; j >= 0; j--)	{
//			if (reverse[i].W[j].count >= 100) { maxReverse = j; break;}		}



		results.close();

}
}


void velocityPeakFinder (vector<xBin> & b)   {

for (int i=0; i < (int)b.size() ; i++)             {
	long max=0;
		for (int j=0; j< (int)b[i].V.size() ; j++ )	{
				if (max < b[i].V[j].count) {
				max = b[i].V[j].count;
				b[i].vPeak= (b[i].V[j].vMax + b[i].V[j].vMin)/2;
			}
		}
	}
}


void workAndVelocityBoundarySetter ( vector<xBin> & b , double fraction)   {

for (int i=0; i < (int)b.size() ; i++)             {
	for (int j=0; j< (int)b[i].W.size() ; j++ )	{
		if (b[i].W[j].count >=  fraction * b[i].wPeakCount)  {
			b[i].wMin = b[i].W[j].wMin;
			break;	}
	}
	for (int j= (int)b[i].W.size() - 1 ; j>=0 ; j-- )	{
		if (b[i].W[j].count >= fraction * b[i].wPeakCount)  {
			b[i].wMax = b[i].W[j].wMax;
			break;	}
	}
}

for (int i=0; i < (int)b.size() ; i++)             {
	for (int j=0; j< (int)b[i].V.size() ; j++ )	{
		if (b[i].V[j].count >= fraction * b[i].vPeakCount)  {
			b[i].vMin = b[i].V[j].vMin;
			break;	}
	}
	for (int j= (int)b[i].V.size() - 1 ; j>=0 ; j-- )	{
		if (b[i].V[j].count >= fraction * b[i].vPeakCount)  {
			b[i].vMax = b[i].V[j].vMax;
			break;	}
	}
}

}





void workAndVelocityPeakFinder (vector<xBin> & b)   {

for (int i=0; i < (int)b.size() ; i++)             {
		for (int j=0; j< (int)b[i].W.size() ; j++ )	{
			if (b[i].wPeakCount < b[i].W[j].count) {
				b[i].wPeakCount = b[i].W[j].count;
				b[i].wPeak= (b[i].W[j].wMax + b[i].W[j].wMin)/2;
			}
		}

		for (int j=0; j< (int)b[i].V.size() ; j++ )	{
				if (b[i].vPeakCount < b[i].V[j].count) {
					b[i].vPeakCount = b[i].V[j].count;
					b[i].vPeak= (b[i].V[j].vMax + b[i].V[j].vMin)/2;
			}
		}
	}

}


void workAndVelocityVectorResetter(vector<xBin> & b)   {

	for (int i=0; i < (int)b.size() ; i++)   {

			b[i].SDwMean = b[i].wPeak = b[i].wMean = b[i].wSkewness = b[i].wPeakCount = b[i].tsCount = b[i].vPeak = b[i].vPeakCount = 0;

			b[i].dw = (b[i].wMax-b[i].wMin) / (double)b[i].numberOfWorkBins;

			for (int j=0; j <  b[i].numberOfWorkBins ; j++)   {
				b[i].W[j].wMin= b[i].wMin + (double)j*b[i].dw;
				b[i].W[j].wMax= b[i].wMin + (double)(j+1)*b[i].dw;
				b[i].W[j].weight = 0;
				b[i].W[j].count = 0;
			}


			b[i].dv = (b[i].vMax-b[i].vMin) / (double)b[i].numberOfVelocityBins;

			for (int j=0; j <  b[i].numberOfVelocityBins ; j++)   {
				b[i].V[j].vMin =  b[i].vMin + (double)j*b[i].dv ;
				b[i].V[j].vMax =  b[i].vMin + (double)(j+1)*b[i].dv ;
				b[i].V[j].count = 0;
				b[i].V[j].weight = 0;
				b[i].V[j].work = 0;
				}

		}

}




void workMeansOutputter (vector<xBin> & f , vector<xBin> & r)  {

	ofstream results;
	results.precision(18);
	results.open("work_distributions/workPeaksAndAverages.txt");
	results<<"x_i  PMF_wMean  PMF_wMode   wForwardMean SDwForward wForwardSkewness  wReverseMean SDwReverse wReverseSkewness \n";
	double PMFmean = 0, PMFmode=0;
	for(int i=0; i<(int)f.size(); i++)	{
			results <<f[i].xi <<"  "<<PMFmean<<"   "<<PMFmode<<"   "<<f[i].wMean<<"   "<<f[i].SDwMean<<"   "<<f[i].wSkewness<<"   "
										                     <<r[i].wMean<<"   "<<r[i].SDwMean<<"   "<<r[i].wSkewness<<endl;
			PMFmean += (f[i].wMean - r[i].wMean)/2;
			PMFmode += (f[i].wPeak - r[i].wPeak)/2;
		}
			results.close();

}




void workAverageCalculator(vector<xBin> & b)   {

	for (int i=0; i < (int)b.size() ; i++)   {
		(b[i].tsCount > 0)  ?   b[i].wMean /= (double) b[i].tsCount  :  b[i].wMean=0;
		}

	for (int i=0; i < (int)b.size() ; i++)   {
		long count = 0;
		double mu3 = 0;
		for (int j=0; j <  b[i].numberOfWorkBins ; j++)   {
		b[i].SDwMean += pow((b[i].W[j].wMin + b[i].W[j].wMax)/2  -  b[i].wMean , 2) *(double) b[i].W[j].count  ;
		mu3 += pow((b[i].W[j].wMin + b[i].W[j].wMax)/2  -  b[i].wMean , 3) *(double) b[i].W[j].count  ;
		count += b[i].W[j].count;
		}
		b[i].SDwMean = sqrt(b[i].SDwMean)/(double) count;
		mu3 /= (double) count;
		b[i].wSkewness =  mu3 / ( pow(b[i].SDwMean * sqrt((double)count) , 3) );

		}


}


void stages ( vector<xBin>  & forward,  vector<xBin> & reverse)	{
	int numberOfStages = 3;
	do {
	cout<<endl<<"Enter the desired number of stages (minimum of 3): ";
	cin>>numberOfStages;
	}
	while (numberOfStages < 3);

	double * fraction = new double[numberOfStages+1];

	for (int i = 2; i<= numberOfStages-1 ; i++)	{
		cout<<"Enter the fraction (of peak value) for stage "<<i<<" : ";
		cin>> fraction[i];
	}

	for (int i=1; i<= numberOfStages ; i++)	{
		if (i == 1)	{
			readInput_BoundaryFinder ( forward,reverse);
			xBinPrepare (forward, reverse);
		}   else if (i < numberOfStages ) 	{
			readInput_sampler ( forward, reverse, i);
			workAndVelocityPeakFinder (forward);
			workAndVelocityPeakFinder (reverse);
			workAndVelocityBoundarySetter (forward, fraction[i]);
			workAndVelocityBoundarySetter (reverse, fraction[i]);
			workAndVelocityVectorResetter(forward);
			workAndVelocityVectorResetter(reverse);
			}	else		{
			readInput_sampler ( forward , reverse , i);
			workAndVelocityPeakFinder (forward);
			workAndVelocityPeakFinder (reverse);

			workAverageCalculator(forward);
			workAverageCalculator(reverse);

			outputVelocityDistributions(forward , reverse);
			outputWorkDistributions( forward , reverse);
			workMeansOutputter (forward , reverse);
		}

	}

}
