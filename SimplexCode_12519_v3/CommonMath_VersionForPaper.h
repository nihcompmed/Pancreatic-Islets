#include <iostream>                  // for std::cout
#include <fstream>
#include <utility>                   // for std::pair
#include <algorithm>                 // for std::for_each
#include <vector>
#include <time.h>
#include <cstdlib>
#include <complex>
#include <cmath>
  
using namespace std;

/* void CalculateAutoCorrelation(vector<vector<double> > * , vector<double> *); */
/* double  CalculateKullback_Leibler(vector<double>*, vector<double>*, double); */
/* void CalculateMeanAndStdev(vector<double> * , double* , double*); */
/* double CalculatePearsonCorrelation(vector<double>*, vector<double>* ); */
/* void CalculatePeriodogram(vector<vector<double> > *, vector<vector<double> > *); */
/* void CalculateSpectrum(vector<vector<double> > *, vector<vector<double> > *,double f_interval=0.5); */
/* void CreateDistribution (vector<double>* , vector<double>*, double, double); */
/* double Distance(vector<double>& ,vector<double>& ); */
/* double DistanceFromPointToLineSegment(double , double , double , double , double , double ); */
 double Distancexy(double ,double , double , double ); 
/* double Distancexyz(double ,double , double , double, double, double ); */
/* void FitExponential(vector<vector<double> > * , double* , double* ); */
 double GetAngleOfLine(double , double , double , double ); 
 bool IntersectionPointBetweenLines(double, double ,double , double ,double , double ,double , double, double *, double *, bool); 
 double mod(double , double ); 
/* int imod(int , int ); */
/* void Shuffle2d(vector<vector<int> >* , vector<vector<int> > *); */
/* bool TwoLinesIntersect(double, double,double , double ,double , double ,double, double); */
