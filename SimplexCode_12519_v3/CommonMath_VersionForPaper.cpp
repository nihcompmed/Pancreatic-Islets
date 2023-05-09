#include "./CommonMath_VersionForPaper.h"
#define PI 3.14159

double Distancexy(double x1,double y1, double x2, double y2){
  return sqrt(pow(x2-x1,2.0)+pow(y2-y1,2));
}

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

double GetAngleOfLine(double p1x, double p1y, double p2x, double p2y){
  //uses p1 as origin and returns the degree a horizontal line has to turn to create the line from p1 to p2 
  double xDiff = p2x - p1x;
  double yDiff = p2y - p1y;
  //cout << "xdiff and ydiff are " << xDiff << " and " << yDiff << endl;
  double degree = atan2(yDiff,xDiff);
  if(degree < 0) degree = 2*PI + degree;
  return degree;
}

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
bool IntersectionPointBetweenLines(double L1x1, double L1y1,double L1x2, double L1y2,double L2x1, double L2y1,double L2x2, double L2y2, double *Ix, double *Iy, bool verbose){
  
  //This function should be used only when v1, v2, v3 or v4 is not the intersection point
  if((L1x1 == L2x1 && L1y1 ==L2y1) ||(L1x1 == L2x2 && L1y1 ==L2y2)||(L1x2 == L2x1 && L1y2 ==L2y1) ||(L1x2 == L2x2 && L1y2 ==L2y2))return 0;
  bool returnval = 0;
  //Find slopes of two lines
  double m1 = 1.0e14;
  if(fabs(L1x2-L1x1) > 1e-8){
    m1 = (L1y2-L1y1)/(L1x2-L1x1);
  }
  double m2 = 1.0e14;
  if(fabs(L2x2-L2x1) > 1e-8){
    m2 = (L2y2-L2y1)/(L2x2-L2x1);
  }
  if(verbose)cout << "The slopes are " << m1 << " and " << m2 << endl;
  if(fabs(m1-m2)<1e-4)return 0;
  else{
    if(m1 == 1e14){
      (*Ix) = L1x1;
      (*Iy) = m2*((*Ix)-L2x1)+L2y1;
    }
    else if(m2 == 1e14){
      (*Ix) = L2x1;
      (*Iy) = m1*((*Ix)-L1x1)+L1y1;
    }
    else{      
      (*Ix) = (L1y1-L2y1+m2*L2x1-m1*L1x1)/(m2-m1);
      (*Iy) = m2*((*Ix)-L2x1)+L2y1;
    }
    if(verbose) cout << "The intersection point is " << (*Ix) << "," << (*Iy) << endl;
   
    //check if x is within line segments
    double L1xmin = 0,L1xmax = 0;
    double L2xmin = 0,L2xmax = 0;
    double L1ymin = 0,L1ymax = 0;
    double L2ymin = 0,L2ymax = 0;
    if(L1x1 < L1x2){
      L1xmin = L1x1;      L1xmax = L1x2;
    }
    else{
      L1xmin = L1x2;      L1xmax = L1x1;
    }
    if(L2x1 < L2x2){
      L2xmin = L2x1;      L2xmax = L2x2;
    }
    else{
      L2xmin = L2x2;      L2xmax = L2x1;
    }
    if(L1y1 < L1y2){
      L1ymin = L1y1;      L1ymax = L1y2;
    }
    else{
      L1ymin = L1y2;      L1ymax = L1y1;
    }
    if(L2y1 < L2y2){
      L2ymin = L2y1;      L2ymax = L2y2;
    }
    else{
      L2ymin = L2y2;      L2ymax = L2y1;
    }
    if(verbose)cout << "For line 1, x ranges from " << L1xmin << " to " << L1xmax << " and y values range from " << L1ymin << " to " <<  L1ymax << endl;
    if(verbose)cout << "For line 2, x ranges from " << L2xmin << " to " << L2xmax << " and y values range from " << L2ymin << " to " <<  L2ymax<< endl;

    //check inf slope for L1x intersections
    if(fabs(L1x2-L1x1) < 1e-8){
      //check if L2x's are on either side
      double L1x = (L1x1+L1x2)/2.0;
      if((L2x1 < L1x && L2x2 > L1x) ||
	 (L2x2 < L1x && L2x1 > L1x)){
	//Check if y is in range
	if((*Iy) >= L1ymin && (*Iy) <= L1ymax &&
	   (*Iy) >= L2ymin && (*Iy) <= L2ymax)returnval =  1;
      }
      }
    else if(fabs(L2x2-L2x1) < 1e-8){
      //check if L1x's are on either side
      double L2x = (L2x1+L2x2)/2.0;
      if((L1x1 < L2x && L1x2 > L2x) ||
	 (L1x2 < L2x && L1x1 > L2x)){
	//Check if y is in range
	if((*Iy) >= L1ymin && (*Iy) <= L1ymax &&
	   (*Iy) >= L2ymin && (*Iy) <= L2ymax)returnval =  1;
      }
    }
    else if(fabs(L1ymin-L1ymax) < 1e-8){//slope = 0
      //check if l2y is in rane
      if((*Iy) >= L2ymin && (*Iy) <= L2ymax){
	//check if x is in range
	if((*Ix) >= L1xmin && (*Ix) <= L1xmax &&
	   (*Ix) >= L2xmin && (*Ix) <= L2xmax)returnval=1;
      }
    }
    else if(fabs(L2ymin-L2ymax) < 1e-8){//slope = 0
      //check if l1y is in rane
      if((*Iy) >= L1ymin && (*Iy) <= L1ymax){
	//check if x is in range
	if((*Ix) >= L1xmin && (*Ix) <= L1xmax &&
	   (*Ix) >= L2xmin && (*Ix) <= L2xmax)returnval=1;
      }
    }
    else if((*Ix) >= L1xmin && (*Ix) <= L1xmax &&
	    (*Ix) >= L2xmin && (*Ix) <= L2xmax &&
	    (*Iy) >= L1ymin && (*Iy) <= L1ymax &&
	    (*Iy) >= L2ymin && (*Iy) <= L2ymax)returnval =  1;
      if(verbose)cout<< "returnval is " << returnval;
      cout.flush();
  }
  return returnval;
}

////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////

double mod(double a, double b)
{ return fmod((fmod(a,b)+b),b); }


////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////





// int imod(int a, int b)
// { return ((a%b)+b)%b; }




// double DistanceFromPointToLineSegment(double px0, double py0, double Lx1, double Ly1, double Lx2, double Ly2){
//   //See if point of normal intersection exists on line segment
//   double m = (Ly2-Ly1)/(Lx2-Lx1);
//   double xint = (Lx1*m + px0/m + py0 - Ly1)/(1/m + m); 
//   if((Lx1 > xint && Lx2 > xint) ||(Lx1 < xint && Lx2 < xint)){
//     double dist1 = Distancexy(px0,py0,Lx1,Ly1);
//     double dist2 = Distancexy(px0,py0,Lx2,Ly2);
//     if(dist1 < dist2)return dist1;
//     else return dist2;
//   }
//   else return fabs(((Ly2-Ly1)*px0-(Lx2-Lx1)*py0 +Lx2*Ly1 - Ly2*Lx1)/Distancexy(Lx1,Ly1,Lx2,Ly2));
// }

// void CalculateMeanAndStdev(vector<double> * data, double* mean, double* stdev){
//   *mean = 0;
//   *stdev = 0;
//   //calculate mean
//   for(int i = 0; i < (*data).size();i++)
//     *mean += (*data)[i];
//   *mean /= (*data).size();

//   for(int i = 0; i < (*data).size();i++)
//     (*stdev)+= pow((*data)[i]-(*mean),2.0);
//   (*stdev) /= (double)((*data).size()-1);
//   return;
// }

// void FitExponential(vector<vector<double> > * data, double* exponent, double* constant){
//   //Algorithm from mathworld.wolfram.com/LeastSquaresFittingExponential.html

//   double sumx2y=0,sumylny=0,sumxy=0,sumxylny=0,sumy=0;
//   for(int i = 0; i < (*data).size();i++){
//     sumx2y   += (*data)[i][0]*(*data)[i][0]*(*data)[i][1];
//     sumylny  += (*data)[i][1]*log((*data)[i][1]);
//     sumxy    += (*data)[i][0]*(*data)[i][1];
//     sumxylny += (*data)[i][0]*(*data)[i][1]*log((*data)[i][1]);
//     sumy     += (*data)[i][1];
//   }
//   //Solving for a where A=exp(a) and A is the constant
//   double anum = sumx2y*sumylny-sumxy*sumxylny;
//   double aden = sumy*sumx2y-pow(sumxy,2.0);
//   *constant = exp(anum/aden);

//   //Solving for b where B=b and B is the exponent
//   double bnum = sumy*sumxylny-sumxy*sumylny;
//   double bden = sumy*sumx2y-pow(sumxy,2.0);
//   *exponent = bnum/bden;
// }

// void CalculatePeriodogram(vector<vector<double> > * data, vector<vector<double> > *results){
//   int N = (*data).size();
//   if((*data).size()%2==0)N--;
  
//   //Calculate mean (a0)
//   double a0=0;
//   for(int i = 0; i < N;i++)a0+=(*data)[i][1];
//   a0=a0/(double)N;

//   for(int i = 0; i < (N-1)/2; i++){
//     double ai = 0,bi = 0;
//     double fi =(double)i/(double)N;
//     for(int t = 0; t < (*data).size();t++){
//       ai+=(*data)[t][1]* cos(2*PI*fi*t);
//       bi+=(*data)[t][1]* sin(2*PI*fi*t);
//     }
//     ai=ai*2.0/(double)N;
//     bi=bi*2.0/(double)N;
//     (*results)[i][0]=fi;
//     (*results)[i][1]=(double)N/2.0*(pow(ai,2)+pow(bi,2));
//   }
//   return;
// }

// void CalculateSpectrum(vector<vector<double> > * data, vector<vector<double> > *results,double f_interval){
//   int N = (*data).size();
//   int nf= (*results).size();
//   double df=f_interval/(double)(nf-1);
  
//   //Calculate mean of data
//   double mean = 0;
//   for(int i = 0; i < (*data).size();i++)mean += (*data)[i][1];
//   mean = mean/(double)(*data).size();

//  //Calculate c0
//   double c0 = 0;
//   for(int t = 0; t < (*data).size();t++)
//     c0+=pow(((*data)[t][1]-mean),2.0);
//   c0 = c0/(double)(*data).size();

//  //Calculate ck*cos
//   for(int f=0;f<nf;f++){
//     (*results)[f][0]=f*df;
//     (*results)[f][1] = c0;
//     double sumtot = 0;
//     for(int k = 1; k < (*data).size()-1;k++){
//       //Calculate ck
//       double ck = 0;
//       for(int t = 1; t< (*data).size()-k; t++)
// 	ck +=((*data)[t][1]-mean)*((*data)[t+k][1]-mean);
//       ck = ck/(double)N;
//       sumtot+= ck*cos(2*PI*f*df*k);
//     }
//     (*results)[f][1]+=2*sumtot;
//   }


//   return;
// }




// void CalculateAutoCorrelation(vector<vector<double> > * data, vector<double> *results){
  
//   //Calculate mean of data
//   double mean = 0;
//   for(int i = 0; i < (*data).size();i++)mean += (*data)[i][1];
//   mean = mean/(double)(*data).size();

//   //Calculate c0
//   double c0 = 0;
//   for(int t = 0; t < (*data).size();t++)
//     c0+=pow(((*data)[t][1]-mean),2.0);
//   c0 = c0/(double)(*data).size();
    
//   //Calculate ck
//   (*results)[0]=1;
//   for(int k = 1; k < (*data).size()/2;k++){
//     for(int t = 0; t< (*data).size()-k; t++)
//       (*results)[k]=((*data)[t][1]-mean)*((*data)[t+k][1]-mean);
//     (*results)[k] = (*results)[k]/(double)(*data).size()/c0;
//   }

//   return;
// }
      


// double CalculatePearsonCorrelation(vector<double>* Data1, vector<double>* Data2){
//   //Calculate means for Data1 and Data2;vectors must be the same size 
//   double mean1 = 0, mean2 = 0;
//   cout << "Data 1: " << (*Data1).size() << " entries and Data 2 " << (*Data2).size() << " entries " ;cout.flush();
//   for(int i = 0; i < (*Data1).size();i++){
//     mean1 += (*Data1)[i];
//     mean2 += (*Data2)[i];
//   }
//   mean1 = mean1/(double)(*Data1).size();
//   mean2 = mean2/(double)(*Data2).size();

//   double sum1 = 0, sum2 = 0, sum3 = 0;
//   for(int i = 0; i < (*Data1).size(); i++){
//     sum1 += ((*Data1)[i]-mean1)*((*Data2)[i]-mean2);
//     sum2 += pow((*Data1)[i]-mean1,2);
//     sum3 += pow((*Data2)[i]-mean2,2);
//   }
//   return sum1/sqrt((sum2*sum3));
// }

// double  CalculateKullback_Leibler(vector<double>* DataDist, vector<double>* RefDist, double binsize){
//   double sum = 0;
//   bool verbose = 0;
//   for(int i = 0; i < (*DataDist).size();i++){
//     if(verbose)cout << "Datadist[" << i << " is " << (*DataDist)[i] << endl;
//     if(verbose)cout << "Refdist[" << i << " is " << (*RefDist)[i] << endl;
//     if((*DataDist)[i]>1e-6)
//       sum += (*DataDist)[i]*log((*DataDist)[i]/(*RefDist)[i])*binsize;
//   } 
//   if(verbose)cout << "sum is " << sum << endl;
//   return sum;
// }
  
// void CreateDistribution (vector<double>* Data, vector<double>* Dist, double mindata,double maxdata){
  
//   int nbins = (*Dist).size();
//   double binsize = (maxdata-mindata)/(double)nbins;
//   vector<int> Hist(nbins,0);
  
//   int ndata = 0;
//   for(int i = 0;i < (*Data).size();i++){
//     if((*Data)[i] == maxdata)Hist[Hist.size()-1]++;
//     else Hist[(int)(((*Data)[i]-mindata)/binsize)]++;
//     ndata++;
//   }
  
//   //Create Distribution
//   for(int i = 0; i < (*Dist).size();i++)
//     (*Dist)[i] = (double)Hist[i]/(double)ndata;

//   //Remove Vector
//   Hist.clear();vector<int> ().swap(Hist);
//   return ;
// }

// void Shuffle2d(vector<vector<int> >*Vin, vector<vector<int> > *Vout){
  
//   //Create Array that is (#rows)*(#columns)
//   vector<int> Varray;
//   for(int i = 0; i < (*Vin).size();i++){
//     for(int j = 0; j < (*Vin)[i].size(); j++){
//       Varray.push_back((*Vin)[i][j]);
//     }
//   }

//   //Shuffle Varray
//   random_shuffle(Varray.begin(),Varray.end());

//   //Printout to Vout
//   int s = 0;
//   for(int i = 0; i < (*Vout).size();i++){
//     for(int j = 0; j < (*Vout)[i].size();j++){
//       (*Vout)[i][j] = Varray[s];
//       s++;
//     }
//   }

//   Varray.clear(); vector<int> ().swap(Varray);
//   return;
// }

// double Distance(vector<double>& p1,vector<double>& p2){
//   return sqrt(pow(p2[0]-p1[0],2.0)+pow(p2[1]-p1[1],2));
// }



// bool TwoLinesIntersect(double L1x1, double L1y1,double L1x2, double L1y2,double L2x1, double L2y1,double L2x2, double L2y2){
//   bool returnval;
//   //Find slopes of two lines
//   double m1 = 100000;
//   if(fabs(L1x2-L1x1) > 1e-4){
//     m1 = (L1y2-L1y1)/(L1x2-L1x1);
//   }
//   double m2 = 100000;
//   if(fabs(L2x2-L2x1) > 1e-4){
//     m2 = (L2y2-L2y1)/(L2x2-L2x1);
//   }
//   if(fabs(m1-m2)<1e-4)return 0;
//   else{
//     double x = (L1y1-L2y1+m2*L2x1-m1*L1x1)/(m2-m1);
   
//     //check if x is within line segments
//     double L1xmin = 0,L1xmax = 0;
//     double L2xmin = 0,L2xmax = 0;
//     if(L1x1 < L1x2){
//       L1xmin = L1x1;
//       L1xmax = L1x2;
//     }
//     else{
//       L1xmin = L1x2;
//       L1xmax = L1x1;
//     }
//     if(L2x1 < L2x2){
//       L2xmin = L2x1;
//       L2xmax = L2x2;
//     }
//     else{
//       L2xmin = L2x2;
//       L2xmax = L2x1;
//     }
//     if(x >= L1xmin && x <= L1xmax &&x >= L2xmin && x <= L2xmax)returnval =  1;
//     else returnval =  0;
//   }
//   return returnval;
// }







  
