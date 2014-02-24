//============================================================================
//
// 	distributionsCalculation.h
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

#ifndef DISTRIBUTIONSCALCULATION_H_
#define DISTRIBUTIONSCALCULATION_H_



class velocityBin {
public:
	double vMin, vMax;
	double work, weight;
	long count;
	velocityBin(double v1, double v2): vMin(v1), vMax(v2), work(0), weight(0), count(0) {};
};

class workBin {
public:
	double wMin, wMax;
	double weight;
	long count;
	workBin(double w1, double w2): wMin(w1), wMax(w2), weight(0), count (0) {};
};




class xBin {

public:
  double xi, xf, vMax, vMin, dv, vPeak, dt, binSize;
  double wMax, wMin, dw, wPeak, wMean, SDwMean, wSkewness;
  long tsCount;
  int numberOfVelocityBins, numberOfWorkBins, wPeakCount, vPeakCount;
  vector <velocityBin> V;
  vector <workBin> W;

  xBin(double x1, double x2, double d_t, int nOfVelocityBins, int nOfWorkBins) {
	  xi=x1;  xf=x2;  tsCount=0;  vMin=3000000.000;  vMax= - 3000000.000;
	  dv=fabs(vMax-vMin)/numberOfVelocityBins;
	  vPeak = wMean = SDwMean = wSkewness = 0;
	  dt = d_t;
	  binSize= xf - xi;
	  wMax = -1000000,  wMin = 1000000, dw = wPeak = 0;
	  numberOfVelocityBins=nOfVelocityBins;
	  numberOfWorkBins=nOfWorkBins;
	  wPeakCount= vPeakCount= 0;
	  vector <velocityBin> V ;
  };

  void velocityVectorGenerator() {
	  dv=(vMax-vMin)/numberOfVelocityBins;
	  for (int i=0; i<numberOfVelocityBins; i++)     {
		  velocityBin *v = new velocityBin(vMin+(double)i*dv, vMin+(double)(i+1)*dv);
          V.push_back(*v);
	  }
  }

  void workVectorGenerator() {
	  dw=(wMax-wMin)/numberOfWorkBins;
	  for (int i=0; i<numberOfWorkBins; i++)     {
		  workBin *w = new workBin(wMin+(double)i*dw, wMin+(double)(i+1)*dw);
          W.push_back(*w);
	  }
  }

};


long file_size ();
void velocityBoundariesDeterminator(vector<xBin> & bin, double X1, double X2);
void xBinPrepare (vector<xBin> & b1, vector<xBin> & b2);
void readInput_BoundaryFinder ( vector<xBin>  & forward,  vector<xBin> & reverse);
void xBinVectorCreator (int numberOfBins, int numberOfWorkBins, int numberOfVelocityBins, double xInit, double binSize, vector<xBin> & bins, double dt);
void readInput_sampler ( vector<xBin>  & forward,  vector<xBin> & reverse, int stage);
void velocityCounter(vector<xBin> & bin, double X1, double X2, double work);
void workCounter(vector<xBin> & bin, double X1, double X2, double work);
void outputVelocityDistributions(vector<xBin> forward, vector <xBin> reverse);
void outputWorkDistributions(vector<xBin> forward, vector <xBin> reverse);
void workAndVelocityPeakFinder (vector<xBin> & b);
void workAndVelocityVectorResetter(vector<xBin> & b);
void workMeansOutputter (vector<xBin> & f , vector<xBin> & r);
void workAndVelocityBoundarySetter ( vector<xBin> & b , double fraction);
void workAverageCalculator(vector<xBin> & b);
void stages ( vector<xBin>  & forward,  vector<xBin> & reverse);



#endif /* DISTRIBUTIONSCALCULATION_H_ */
