#include "./GraphNN_VersionForPaper.h"

#define PI 3.14159

using namespace boost;
void CalcGofr(vector<IsletGraph>* Gp,vector<vector<double> > *Gofr, string celltype,int Gno){ 
  std::pair<vertex_iter, vertex_iter>vp,vp1,vp2;
  
  //if(Gno != 99999)cout << "Calculating Gofr for Islet " << Gno << " out of " << (*Gp).size() << "...";//algorithm can be run for all Islets in a set using Gno = 99999
  //else cout << "Calculating Gofr ..." << endl;
  int vectorcount = 0;
  bool verbose=0,verbosemem=0;
  //Calculate longest distance between cells
  double maxdist = 0;
  for(int i = 0; i < num_vertices((*Gp)[Gno].allGraph);i++){
    for(int j = i+1; j < num_vertices((*Gp)[Gno].allGraph);j++){
      double distance = Distancexy((*Gp)[Gno].allGraph[i].x, (*Gp)[Gno].allGraph[i].y,(*Gp)[Gno].allGraph[j].x, (*Gp)[Gno].allGraph[j].y);
      if(distance > maxdist)maxdist = distance;
    }
  }
  double dr = 0.5;
  int nbins = maxdist/dr + 4.0;
  
  vector<vector<double> > gk_norm((*Gp).size(),vector<double>(nbins,0.0));vectorcount++;if(verbosemem)cout << "Adding gk_norm: vectorcount = " << vectorcount << endl; cout.flush();
  int numIslets = 0;
  int Gbeg = 0, Gend = (*Gp).size();
  if(Gno!=99999){
    Gbeg = Gno;
    Gend = Gno+1;
  }

  for(int i = Gbeg; i < Gend;i++){
    if(verbose &&  i%100 == 0){
      cout << i << " " ;
      cout.flush();
    }
    
    //Determine the (#vertices)^2
    int numvertsq = 0;
    if (celltype == "abd")  numvertsq = (int)pow(num_vertices((*Gp)[i].allGraph),2.0);
    else if(celltype == "aa")numvertsq = (int)pow((*Gp)[i].nalpha,2.0);
    else if(celltype == "bb")numvertsq = (int)pow((*Gp)[i].nbeta,2.0);
    else if(celltype == "dd")numvertsq = (int)pow((*Gp)[i].ndelta,2.0); 
    else if(celltype == "adm")numvertsq = (int)pow((*Gp)[i].nalpha+(*Gp)[i].ndelta,2.0);
    else if(celltype == "ab")numvertsq = (*Gp)[i].nalpha*(*Gp)[i].nbeta;
    else if(celltype == "ad")numvertsq = (*Gp)[i].nalpha*(*Gp)[i].ndelta;
    else if(celltype == "bd")numvertsq = (*Gp)[i].ndelta*(*Gp)[i].nbeta;
    else if(celltype == "adb")numvertsq = ((*Gp)[i].nalpha+(*Gp)[i].ndelta)*(*Gp)[i].nbeta;
    
    if((celltype == "abd" && num_vertices((*Gp)[i].allGraph) > 1)||
       (celltype == "aa" && (*Gp)[i].nalpha > 1)||
       (celltype == "bb" && (*Gp)[i].nbeta > 1)||
       (celltype == "dd" && (*Gp)[i].ndelta > 1)||
       (celltype == "adm" && (*Gp)[i].nalpha+(*Gp)[i].ndelta> 1)||
       (celltype == "ab" && (*Gp)[i].nalpha > 0 && (*Gp)[i].nbeta > 0) ||
       (celltype == "ad" && (*Gp)[i].nalpha > 0 && (*Gp)[i].ndelta > 0) ||
       (celltype == "bd" && (*Gp)[i].nbeta > 0 && (*Gp)[i].ndelta > 0) ||
       (celltype =="adb" && ((*Gp)[i].nalpha > 0 ||(*Gp)[i].ndelta >0) && (*Gp)[i].nbeta > 0)){
      numIslets++;
      
      //Determine bounding box for islet area
      double IsletArea = 0;
      double xlength = 0,ylength = 0;
      double xmin = 100000, xmax = 0, ymin = 100000,ymax = 0;
      for(vp = vertices((*Gp)[i].allGraph);vp.first!=vp.second;vp.first++){
	if(((*Gp)[i].allGraph)[*vp.first].x >xmax)xmax = ((*Gp)[i].allGraph)[*vp.first].x;
	if(((*Gp)[i].allGraph)[*vp.first].y >ymax)ymax = ((*Gp)[i].allGraph)[*vp.first].y;
	if(((*Gp)[i].allGraph)[*vp.first].x <xmin)xmin = ((*Gp)[i].allGraph)[*vp.first].x;
	if(((*Gp)[i].allGraph)[*vp.first].y <ymin)ymin = ((*Gp)[i].allGraph)[*vp.first].y;
      }
      xlength = xmax-xmin;
      ylength = ymax-ymin;
      IsletArea = xlength*ylength;
      
      if(verbose)cout << "For Islet " << (*Gp)[i].IsletNum <<"\t" <<  xmin << "\t" << xmax << "\t" << xlength << "\t" <<  ymin << "\t" << ymax << "\t" << ylength << endl;
      
      //calculate the number of possible edges
      for(vp1=vertices((*Gp)[i].allGraph);vp1.first!=vp1.second;vp1.first++){
	for(vp2=vertices((*Gp)[i].allGraph);vp2.first!=vp2.second;vp2.first++){
	  if((*vp1.first != *vp2.first)&&
	     (celltype == "abd"|| 
	      (celltype == "aa" && (*Gp)[i].allGraph[*vp1.first].type == "a" && (*Gp)[i].allGraph[*vp2.first].type == "a") || 
	      (celltype == "bb" && (*Gp)[i].allGraph[*vp1.first].type == "b" && (*Gp)[i].allGraph[*vp2.first].type == "b") || 
	      (celltype == "dd" && (*Gp)[i].allGraph[*vp1.first].type == "d" && (*Gp)[i].allGraph[*vp2.first].type == "d") || 
	      (celltype == "ab" && (*Gp)[i].allGraph[*vp1.first].type == "a" && (*Gp)[i].allGraph[*vp2.first].type == "b") ||
	      (celltype == "ad" && (*Gp)[i].allGraph[*vp1.first].type == "a" && (*Gp)[i].allGraph[*vp2.first].type == "d") ||
	      (celltype == "bd" && (*Gp)[i].allGraph[*vp1.first].type == "d" && (*Gp)[i].allGraph[*vp2.first].type == "b") ||
	      (celltype == "adb" && ((*Gp)[i].allGraph[*vp1.first].type == "a" ||(*Gp)[i].allGraph[*vp1.first].type == "d") &&(*Gp)[i].allGraph[*vp2.first].type == "b")  ||
	      (celltype == "adm" && 
	       ((*Gp)[i].allGraph[*vp1.first].type == "a" || (*Gp)[i].allGraph[*vp1.first].type == "d")  &&
	       ((*Gp)[i].allGraph[*vp2.first].type == "a" || (*Gp)[i].allGraph[*vp2.first].type == "d")))) {
	    
	    double x1 = ((*Gp)[i].allGraph)[*vp1.first].x, y1 = ((*Gp)[i].allGraph)[*vp1.first].y;
	    double x2 = ((*Gp)[i].allGraph)[*vp2.first].x, y2 = ((*Gp)[i].allGraph)[*vp2.first].y;
	      
	    double radius = Distancexy(x1,y1,x2,y2);
	    int arc = 0;
	    if(verbose)cout << "For vertex "<<(*Gp)[i].allGraph[*vp1.first].IsletCellNum << " (" << (*Gp)[i].allGraph[*vp1.first].x << 
	      "," << (*Gp)[i].allGraph[*vp1.first].y << ") with type " << (*Gp)[i].allGraph[*vp1.first].type << "  - " << 
	      (*Gp)[i].allGraph[*vp2.first].IsletCellNum << " (" << (*Gp)[i].allGraph[*vp2.first].x << "," << (*Gp)[i].allGraph[*vp2.first].y << 
	      ") with type " << (*Gp)[i].allGraph[*vp2.first].type << " with radius " << radius << " :" << endl;
	    //Determine how much of the circumference of circle, with center = cell 1 and radius = distance between cell 1 and cell2, is within the bounding box
	    //First quadrant
	    for(int deg = 0;deg < 90;deg++){
	      double theta = (double)deg*PI/180.0;
	      double xs = x2+radius*cos(theta);
	      double ys = y2+radius*sin(theta);
	      if(verbose)cout << " for source deg " << deg << ": (x,y) is (" << xs << "," << ys << ")" << endl;
	      if(xs <= xmax && xs >= xmin && ys <= ymax && ys >= ymin)
		arc++;
	      if(verbose)cout << "  Adding one to arc" << endl;
	    }
	    if(verbose)cout << ".....After 1st quad: arc is " << arc <<  endl;
	    //2nd quadrant
	    for(int deg = 0;deg < 90;deg++){
	      double theta = (double)deg*PI/180.0;
	      double xs = x2-radius*cos(theta);
	      double ys = y2+radius*sin(theta);
	      if(verbose)cout << " for source deg " << deg << ": (x,y) is (" << xs << "," << ys << ")" << endl;
	      if(xs <= xmax && xs >= xmin && ys <= ymax && ys >= ymin)
		arc++;
	      if(verbose)cout << " Adding one to arc" << endl;
	      
	      if(verbose)cout << " for deg " << deg << ": (x,y) is (" << xs << "," << ys << ")" << endl;
	    }
	    if(verbose)cout << ".....After 2ND quad: arc is " << arc << endl;
	    //3rd quadrant
	    for(int deg = 0;deg < 90;deg++){
	      double theta = (double)deg*PI/180.0;
	      double xs = x2-radius*cos(theta);
	      double ys = y2-radius*sin(theta);
	      if(verbose)cout << " for source deg " << deg << ": (x,y) is (" << xs << "," << ys << ")" << endl;
	      if(xs <= xmax && xs >= xmin && ys <= ymax && ys >= ymin)
		arc++;
	      if(verbose)cout << " Adding one to arc" << endl;
	    }
	    if(verbose)cout << ".....After 3rd quad: arc is " << arc <<  endl;
	    //4th quadrant
	    for(int deg = 0;deg < 90;deg++){
	      double theta = (double)deg*PI/180.0;
	      double xs = x2+radius*cos(theta);
	      double ys = y2-radius*sin(theta);
	      if(verbose)cout << " for source deg " << deg << ": (x,y) is (" << xs << "," << ys << ")" << endl;
	      if(xs <= xmax && xs >= xmin && ys <= ymax && ys >= ymin)
		arc++;
	      if(verbose)cout << " Adding one to arc" << endl;
	    }
	    if(verbose)cout << ".....After 4th quad: arc is " << arc <<   endl;
	    double arcdegrees = (double)arc*PI/180.0;
	    if(verbose)cout << "...Arcdegrees = " << arcdegrees <<  endl;
	    
	    if(arc !=0){
	      if((int)(radius/dr)>=gk_norm[i].size()-1){
		        if(verbose)cout << "radius/dr= " << radius/dr << "gk_norm.size is " << gk_norm[i].size()-1 << endl;
	      }
	      else gk_norm[i][(int)(radius/dr)]+=IsletArea/(arcdegrees*radius*dr*numvertsq);
	    }
	  }
	}//closes vp2 
      }//closes vp1
      if(verbose){
	        //Printout gknorm values
	        cout << "For islet " << i << ": gknorm values are " ;
	        for(int j = 0; j < gk_norm[i].size();j++){
	          cout << gk_norm[i][j] << " ";
	        }
	        cout << endl;
      }
    }
  }//closes Graphs
  
  for(int j = 0; j < nbins;j++){
    vector<double> avggk(2,0);vectorcount++;if(verbosemem)cout << "Adding avggk: vectorcount = " << vectorcount << endl; cout.flush();
    avggk[0]=dr*j;
    for(int i = 0; i < (*Gp).size();i++){
      avggk[1]+= gk_norm[i][j];//i represents the islet number; if single islet, i=0
    }
    if(numIslets!=0)avggk[1] = avggk[1]/(double)numIslets;
    else avggk[1]=0;
    (*Gofr).push_back(avggk);
    RemoveVector(&avggk);vectorcount--;if(verbosemem)cout << "Removing avggk: vectorcount = " << vectorcount << endl; cout.flush();
  }
  RemoveVector2d(&gk_norm);vectorcount--;if(verbosemem)cout << "Removing gknorm: vectorcount = " << vectorcount << endl; cout.flush();
  //if(vectorcount > 0)cout << "Problem with vectors! vectorcount = " << vectorcount << " in CalcGofr()" << endl;
  //cout << " Finished" << endl;
  return ;
}

//*********************************************************************************************//
//*********************************************************************************************//
//*********************************************************************************************//

void DetermineThreshold(vector<IsletGraph>* I,string type, bool printgr, string fig_directory, string fileprefix){
  int vectorcount = 0;
  bool verbose=0,verbosemem=0;
  
  int numIslets = 0;
  int fignum=0;
  
  //check sorting of Islets
  //check if biggest first
  bool biggestFirst = 1;
  for(int i = 0;i < (*I).size()-1;i++){
    if(num_vertices((*I)[i+1].allGraph) <num_vertices((*I)[i+1].allGraph)){
      biggestFirst = 0;
      break;
    }
  }
  //check if smallest first
  bool smallestFirst = 1;
  for(int i = 0;i < (*I).size()-1;i++){
    if(num_vertices((*I)[i+1].allGraph) >num_vertices((*I)[i+1].allGraph)){
      smallestFirst = 0;
      break;
    }
  }

  //sort islets based on IsletNumber to segregate subjects
  sort((*I).begin(), (*I).end(), compareByIsletNumSmallestFirst);
  
  int i = 0;
  int subjno = 1;
  while(i < (*I).size()){
    //find islets of each subject
    int numIslets = 0;
    int subjbegin = i;
    vector<double> MinBetweenPeaks;vectorcount++;if(verbosemem)cout << "Adding MinBetweenPeaks: vectorcount = " << vectorcount << endl; cout.flush();

    //while((*I)[i].IsletNum/10000 == subjno){

      if((type =="bb" && (*I)[i].nbeta > 5) || 
	 (type == "adm" &&  (*I)[i].nalpha+(*I)[i].ndelta > 5) ||
	 (type == "adb" && ((*I)[i].nbeta > 5 && (*I)[i].nalpha+(*I)[i].ndelta > 5))){
	numIslets++;
	
	//Calculate G(r) for each individual islet in the subject
	vector<vector<double> > Gofr; vectorcount++;if(verbosemem)cout << "Adding Gofr: vectorcount = " << vectorcount << endl; cout.flush();
	
	CalcGofr(I,&Gofr,type,i);
	
	//Smooth vector twice
	vector<vector<double> > Smooth(Gofr.size(),vector<double>(2,0));vectorcount++;if(verbosemem)cout << "Adding Smooth: vectorcount = " << vectorcount << endl; cout.flush();
	vector<vector<double>  >Smooth2(Gofr.size(),vector<double>(2,0));vectorcount++;if(verbosemem)cout << "Adding Smooth2: vectorcount = " << vectorcount << endl; cout.flush();
	SmoothVector5(&Gofr,&Smooth2); SmoothVector5(&Smooth2,&Smooth);
	
	//Calculate Peaks
	vector<vector<double>  >Peaki;vectorcount++;if(verbosemem)cout << "Adding Peaki: vectorcount = " << vectorcount << endl; cout.flush();
	
	PeakFinder(&Smooth,&Peaki);
	
	if(verbose)cout << "before MinPeaks: There are " << Peaki.size() << " peaks" << endl; 
	if(verbose){
	  cout << "for type " << type << " with " << (*I)[i].nalpha << " alphas, " << (*I)[i].nbeta << " betas, and " << (*I)[i].ndelta << " deltas " << endl;
	  
	}
	
	double min2_3 = MinAfter2ndHighestPeak(&Smooth,&Peaki);
	if(Peaki.size() > 2)MinBetweenPeaks.push_back(min2_3);
	if(type == "bb")(*I)[i].bbthresh = min2_3;
	else if(type == "adm")(*I)[i].admthresh = min2_3;
	else(*I)[i].adbthresh = min2_3;

	double totalcells = (double)((*I)[i].nalpha+(*I)[i].nbeta+(*I)[i].ndelta);
	if(printgr){
	    PrintoutGnuplotGofr4x5(&Smooth, &((*I)[i]), fig_directory, fileprefix, type, min2_3);
	}

      //Remove vectors
      RemoveVector2d(&Gofr);vectorcount--;if(verbosemem)cout << "Removing Gofr: vectorcount = " << vectorcount << endl; cout.flush();
      RemoveVector2d(&Smooth);vectorcount--;if(verbosemem)cout << "Removing Smooth: vectorcount = " << vectorcount << endl; cout.flush();
      RemoveVector2d(&Smooth2);vectorcount--;if(verbosemem)cout << "Removing Smooth2: vectorcount = " << vectorcount << endl; cout.flush();
      RemoveVector2d(&Peaki);vectorcount--;if(verbosemem)cout << "Removing Peaki: vectorcount = " << vectorcount << endl; cout.flush();
      }
      i++;

    //}
    
    double thresh = 0;
    if(MinBetweenPeaks.size() > 0){
      //Calculate mean 
      for(int p = 0; p < MinBetweenPeaks.size();p++)thresh += MinBetweenPeaks[p];
      //Printout individual islet's threshold
      thresh/=(double)MinBetweenPeaks.size();
    }
    //check subject threshold values: if 0, replace with average value or 20
    for(int check_i = subjbegin; check_i < i+1; check_i++){
      if(type =="adm" && (*I)[check_i].admthresh < 1e-4){//islet hasn't been set
	if(thresh > 0)
	  (*I)[check_i].admthresh = thresh;
	else (*I)[check_i].admthresh = 20;
      }
      if(type == "bb" && (*I)[check_i].bbthresh < 1e-4){
	if(thresh > 0)
	  (*I)[check_i].bbthresh = thresh;
	else (*I)[check_i].bbthresh = 20;
      }
      if(type =="adb" && (*I)[check_i].adbthresh < 1e-4){
	if(thresh > 0)
	  (*I)[check_i].adbthresh = thresh;
	else (*I)[check_i].adbthresh = 20;
      }
    }
    RemoveVector(&MinBetweenPeaks);vectorcount--;
    subjno++;
  }
  
  if(verbose)cout << "Finished all peaks for " << numIslets << " islets \n";

  
  if(vectorcount != 0)cout << "Problem with vectors! vectorcount = " << vectorcount << " in DetermineThreshold()" << endl;
  return;

}

//*********************************************************************************************//
//*********************************************************************************************//
//*********************************************************************************************//

void FindCellShadows(IsletGraph * I, vertex_t uall, vertex_t vall, vector<double>*shadowRange, double shadowthresh){
  bool verbose=0,verbosemem=0;
  
  if(verbose)cout << endl << endl << "Inside FindCellShadows" << endl;

  //Calculate angle of line from uall to vall
  double Angle = GetAngleOfLine((*I).allGraph[uall].x,(*I).allGraph[uall].y,(*I).allGraph[vall].x,(*I).allGraph[vall].y);
  double distance = Distancexy((*I).allGraph[uall].x,(*I).allGraph[uall].y,(*I).allGraph[vall].x,(*I).allGraph[vall].y);
  double AnglePlus = mod(Angle + atan(shadowthresh/distance), 2*PI);
  double AngleMinus = mod(Angle - atan(shadowthresh/distance), 2*PI);

  if(verbose){
    cout << "The shadowrange before is " << endl;
    cout.flush();
    for(int s = 0; s < (*shadowRange).size();s=s+2)
      cout << "[" << (*shadowRange)[s]*180.0/PI << "," << (*shadowRange)[s+1]*180.0/PI << "] ";
    cout << endl;
    cout.flush();
    cout << "The angle between vertex " << (*I).allGraph[uall].IsletCellNum << " and " << (*I).allGraph[vall].IsletCellNum << " is " << Angle*360.0/(2*PI) << endl;
    cout << "The distance is " << distance << endl;
    cout << " and angleminus and angleplus are " << AngleMinus*180.0/PI << "," << AnglePlus*180.0/PI << endl;
  }

  
  //Remove area from shadowRange
  if(AngleMinus<AnglePlus){//segment doesn't wrap
    int ampos = 0;
    while(ampos < (*shadowRange).size() && AngleMinus > (*shadowRange)[ampos] )ampos++;
    int appos = ampos;
    while(appos < (*shadowRange).size() && AnglePlus > (*shadowRange)[appos]  )appos++;
    if(verbose) cout << "ampos and appos are " << ampos << " and " << appos << endl;

    //Remove entries from ampos to and including appos
    (*shadowRange).erase((*shadowRange).begin()+ampos, (*shadowRange).begin()+appos);
    if(verbose){
      cout << "The shadowrange after removing in between ampos and appos is " << endl;
      cout.flush();
      for(int s = 0; s < (*shadowRange).size();s=s+2)
  	cout << "[" << (*shadowRange)[s]*180.0/PI << "," << (*shadowRange)[s+1]*180.0/PI << "] ";
      cout << endl;
      cout.flush();
    }
    if(ampos%2 == 1){
      (*shadowRange).insert((*shadowRange).begin()+ampos, AngleMinus);
      ampos++;
    }
    if(appos%2 == 1)(*shadowRange).insert((*shadowRange).begin()+ampos, AnglePlus);
  }
  else{
    int appos = 0;
    while(appos < (*shadowRange).size() && AnglePlus > (*shadowRange)[appos])appos++;
    int ampos = appos;
    while(ampos < (*shadowRange).size() && AngleMinus > (*shadowRange)[ampos])ampos++;
    if(verbose) cout << "ampos and appos are " << ampos << " and " << appos << endl;

    //Remove entries from ampos to end and from beginning to appos
    (*shadowRange).erase((*shadowRange).begin()+ampos, (*shadowRange).end());
    if(ampos%2 == 1)(*shadowRange).insert((*shadowRange).begin()+ampos, AngleMinus);
    if(verbose){
      cout << "The shadowrange after updating ampos is " << endl;
      cout.flush();
      for(int s = 0; s < (*shadowRange).size();s=s+2)
  	cout << "[" << (*shadowRange)[s]*180.0/PI << "," << (*shadowRange)[s+1]*180.0/PI << "] ";
      cout << endl;
      cout.flush();
    }
    (*shadowRange).erase((*shadowRange).begin(), (*shadowRange).begin()+appos);
    if(appos%2 == 1)(*shadowRange).insert((*shadowRange).begin(), AnglePlus);
    if(verbose){
      cout << "The shadowrange after updating appos is " << endl;
      cout.flush();
      for(int s = 0; s < (*shadowRange).size();s=s+2)
  	cout << "[" << (*shadowRange)[s]*180.0/PI << "," << (*shadowRange)[s+1]*180.0/PI << "] ";
      cout << endl;
      cout.flush();
    }
  }	

  //  if(verbose){
  //   cout << "The shadowrange after is " << endl;
  //   cout.flush();
  //   for(int s = 0; s < (*shadowRange).size();s=s+2)
  //     cout << "[" << (*shadowRange)[s]*180.0/PI << "," << (*shadowRange)[s+1]*180.0/PI << "] ";
  //   cout << endl;
  //   cout.flush();
    
  // }
   //no vectors,graphs,paths or cells 
   return;
}
//*********************************************************************************************//
//*********************************************************************************************//
//*********************************************************************************************//


void FindInteriorPath(IsletGraph* Ip, Graph* Gp ,vertex_t startv, vector<vertex_t>*Path , vertex_t *nextv , int pathnum,double cellx,double celly){
   std::pair<out_edge_iter, out_edge_iter> oep;
   vertex_t u;
  bool verbose=0,verbosemem=0;
  int vectorcount = 0;
  if(verbose)PrintoutGnuplot(Gp,"tempFindIntPath");
  if(verbose)cout << "Finding Interior Path" << endl;
  
  //let (cellx,celly) be the beta cell inside the loop (or startv from the previous path) and startv be the first ad cell in ring
  (*Path).push_back(startv);
  vector <CellAngle> CAngles;vectorcount++;if(verbosemem)cout << "Adding CAngles: vectorcount = " << vectorcount << endl; cout.flush();
  if(verbose)cout << "Cell " << (*Gp)[startv].IsletCellNum << " has " << out_degree(startv,(*Gp)) <<" neighbor" <<  endl;cout.flush();

  MakeGraphLocallyPlanar(Ip,Gp,startv);

  for(oep = out_edges(startv,(*Gp));oep.first!=oep.second; oep.first++){
    CellAngle temp;
    temp.Cell = target(*oep.first,(*Gp));
    temp.Angle = GetAngleOfLine((*Gp)[startv].x,(*Gp)[startv].y, (*Gp)[target(*oep.first,(*Gp))].x,(*Gp)[target(*oep.first,(*Gp))].y);
    CAngles.push_back(temp);
  }
  sort(CAngles.begin(),CAngles.end(),compareByAngle);
  if(verbose){
    cout << "After sorting: cell angles are " << endl;
    for(int i =0; i < CAngles.size();i++){
      cout << (*Gp)[CAngles[i].Cell].IsletCellNum << "\t" << CAngles[i].Angle << endl;
    }
  }
  //Find cell located after angle between startv and nextv
  //if(verbose)cout << "Before finding b_angle:The points are (" << endl;
  double b_angle = GetAngleOfLine((*Gp)[startv].x,(*Gp)[startv].y,cellx,celly);
  if(verbose)cout << "Beta angle is " << b_angle << endl;
  int CApos = 0; //if b_angle is gt all CellAngles then next is CellAngle[0]
  for(int i = 0; i < CAngles.size();i++){
    if(verbose)cout << "CAngles[" << i << "].Cell is " <<(*Gp)[CAngles[i].Cell].IsletCellNum << endl;
    cout.flush();
    if(CAngles[i].Angle > b_angle){
      CApos = i;
      break;
    }
  }

  // nextv and Path[-1] are the same
  (*Path).push_back(CAngles[CApos].Cell);
  (*nextv) = CAngles[CApos].Cell;


  //Remove CAngles
  RemoveVector(&CAngles);vectorcount--;if(verbosemem)cout << "Removing CAngles: vectorcount = " << vectorcount << endl; cout.flush();
  
  double MidAngle = 0,MidAngleBegin = 0,MidAngleEnd = 99999;
  bool foundpath = 0;

  while(!foundpath){

    // Run this loop till you come back to starting point
    while((*Path).size() == 1 ||(*Path)[(*Path).size()-1]!=startv){
      // ||((*Path)[(*Path).size()-1]==startv && PathMidAngle[0] !=PathMidAngle[PathMidAngle.size()-1]) ){
      FindNextCellInPath(Ip,Gp, Path, &MidAngle,pathnum,nextv);

      // Why is this there? -- Manu
      if((*Path).size()==3) MidAngleBegin = MidAngle;

      if(verbose){
	        cout << "The newPath around the cluster is " ;
	        for(int i = 0; i < (*Path).size();i++)cout << (*Gp)[(*Path)[i]].IsletCellNum <<" ";
	        cout << endl;
      }
    }

    FindNextCellInPath(Ip,Gp, Path, &MidAngle,pathnum,nextv);
    FindNextCellInPath(Ip,Gp, Path, &MidAngle,pathnum,nextv);
    if(MidAngle == MidAngleBegin){
      foundpath = 1;
      (*Path).pop_back();
      (*Path).pop_back();
    }

  }

 
  
  if(verbose){
    cout << "The newPath around the cluster is " ;
    for(int i = 0; i < (*Path).size();i++)cout << (*Gp)[(*Path)[i]].IsletCellNum <<" ";
    cout << endl;
    cout << "nextv is cell " << (*Gp)[*nextv].CellNum << endl;
  }
  if(vectorcount != 0) cout << "Problem with vectors! vectorcount = " << vectorcount << " in FindInteriorEdges" << endl;
  return;
}

//*********************************************************************************************//
//*********************************************************************************************//
//*********************************************************************************************//

void FindPartialLoops( IsletGraph * I, double shadowthresh,string fprefix){
  
  bool verbose=0,verbosemem=0;
  std::pair<vertex_iter, vertex_iter> vp;
  std::pair<out_edge_iter, out_edge_iter> oep; 
  std::pair<edge_iter, edge_iter> ep; 
  int graphcount = 0, vectorcount = 0, Igraphcount = 0;; 
  if(verbose)cout << "vectorcount = " << vectorcount << endl;
  if(verbose)PrintoutGnuplotExteriorPaths(I,fprefix+"orig.exteriorpaths");
  if(verbose)PrintoutGnuplotSetLoops(I,0 ,fprefix+".2xgr.setLoops");
  if(verbose)cout << "Inside FindPartialLoops for Islet " << (*I).IsletNum << " with " << num_vertices((*I).allGraph) << " vertices" << endl;
    

  verbose=0;
  //Find ranges of admissable edges
  vector<string> types;vectorcount++;if(verbosemem)cout << "Adding types: vectorcount = " << vectorcount << endl;
  types.push_back("b"); types.push_back( "ad");
  
  vector<vector<vertex_t> > *Cp, *Sp;
  vector<vector<vertex_t> > * Ep1, *Ep2,*Lp, *Mp;
  Graph * Gp1, * Gp2 ;
  vector<int> * Mtp;// *Mtncomp;
  vector<vector<double> > * MtPerc;
  edge_t e; bool b;

  for(int t = 0; t < types.size();t++){
    if(types[t] == "b"){
      Ep1 = &((*I).bExteriorPaths);
      Ep2 = &((*I).adExteriorPaths);
      Gp1 = &((*I).bGraph);
      Gp2 = &((*I).adGraph);
      Mp = &((*I).bMantle);
      Lp = &((*I).betasetLoops);
      Cp = &((*I).bComponents);
      Mtp = &((*I).bMantle_type);
      Sp = &((*I).betasets);
      MtPerc = &((*I).bMantlePercent);
      // Mtncomp = &((*I).bMantle_ncomp);
    }
     
    if(types[t] == "ad"){
      Ep1 = &((*I).adExteriorPaths);
      Ep2 = &((*I).bExteriorPaths);
      Gp1 = &((*I).adGraph);
      Gp2 = &((*I).bGraph);
      Mp = &((*I).adMantle);
      Lp = &((*I).adsetLoops);
      Cp = &((*I).adComponents);
      Mtp = &((*I).adMantle_type);
      Sp = &((*I).adsets);
      MtPerc = &((*I).adMantlePercent);
      //Mtncomp = &((*I).adMantle_ncomp);
    }

    //Loop through components and add in original betasets
    vector<int> vertsInSets;vectorcount++;if(verbosemem)cout << "Adding vertsInSets: vectorcount = " << vectorcount << endl; cout.flush();
    for(int c1r = 0; c1r < (*Cp).size();c1r++){
      if(verbose){
  	      cout << "Component " << c1r << " contains " <<(*Cp)[c1r].size() << " cells and is " << endl;
  	      for(int cc  = 0;cc < (*Cp)[c1r].size();cc++)
  	        cout << (*Gp1)[(*Cp)[c1r][cc]].IsletCellNum << " ";
  	      cout << endl;
  	      cout << "The sets contain: " << endl;
  	      for(int s = 0; s < (*Sp).size();s++){
  	        cout << s << ": "; 
  	        for(int ss = 0; ss < (*Sp)[s].size();ss++)cout << (*Gp1)[(*Sp)[s][ss]].IsletCellNum << " ";
  	        cout << endl;
  	      }
      }
      //Check if c1r is in an original set
      bool inset = 0;
      if(verbose)cout << "Sp size is " << (*Sp).size() << endl;
      for(int i = 0; i < (*Sp).size();i++){
  	      if(verbose)cout << "i = " << i << " and Sp i .size is " << (*Sp)[i].size() << endl;
  	      for(int j = 0; j < (*Sp)[i].size();j++){
  	        //check if vert is already in vertsInSets
  	        bool vinSet = 0;
  	        for(int vv = 0; vv < vertsInSets.size();vv++){
  	          if((*Sp)[i][j] == vertsInSets[vv]){
  	            vinSet = 1;
  	            break;
  	          }
  	        }
  	        if(!vinSet)vertsInSets.push_back((*Sp)[i][j]);
  	        //if (verbose) cout<< "  j = " << j << " and Cp[c1r] size is " << (*Cp)[c1r].size() << endl;
  	        for(int ci = 0; ci < (*Cp)[c1r].size();ci++){
  	          //if (verbose)cout << "  ci = " << ci << endl;
  	          if((*Sp)[i][j] == (*Cp)[c1r][ci]){
  	            if(verbose)cout << "Cell " << (*Gp1)[(*Sp)[i][j]].IsletCellNum << " is = to " <<  (*Gp1)[(*Cp)[c1r][ci]].IsletCellNum  << endl;
  	            inset = 1;
  	            (*Mp)[c1r].insert((*Mp)[c1r].begin(),(*Lp)[i].begin(),(*Lp)[i].end());
  	            (*Mtp)[c1r] = 1;
	            
  	            for(int cc = 0; cc < (*Ep1)[c1r].size(); cc++){
		                if((*Gp1)[(*Ep1)[c1r][cc]].IsletCellNum < 100000)(*MtPerc)[c1r].push_back(1.0);
		                else (*MtPerc)[c1r].push_back(0.0);
	            }
	           
	            //  	      (*Mtncomp)[c1r] = 1;
  	            if(verbose){
  	      	        if (verbose)cout << "Component " << c1r << " is in an original set.   " << endl;
  	      	        cout << "Its Mantle consists of " << endl;cout.flush();
  	      	        for(int m = 0; m < (*Mp)[c1r].size();m++){
  	      	          //if((*Mp)[c1r][m] < 100000)
		                cout << (*Gp2)[(*Mp)[c1r][m]].IsletCellNum << " ";cout.flush();
		                  //else
  	      	            //cout << (*Mp)[c1r][m] << " " ;
		                
  	      	        }
  	      	        cout << endl << endl;
  	            }
	            
  	            break;
  	          }
  	          //if (verbose)cout << "After break: i = " << i << " j = " << j << " ci = " << ci << endl;
  	        }
  	        if(inset)j = (*Sp)[i].size();
  	        //if (verbose) cout << "After break2: i = " << i << " j = " << j <<  endl;
  	      }
  	      if(inset)i = (*Sp).size();
  	      //if (verbose)	cout << "After break3: i = " << i <<  endl;
      }
    }
   
    //Check if all vertices are in sets
    bool allInSet = 0;
    if(types[t] == "b" && vertsInSets.size() == (*I).nbeta)allInSet = 1;
    if(types[t] == "ad" && vertsInSets.size() == (*I).nalpha + (*I).ndelta)allInSet = 1;
    RemoveVector(&vertsInSets);vectorcount--;if(verbosemem)cout << "Removing vertsInSets: vectorcount = " << vectorcount << endl; cout.flush();

     if(!allInSet){
      //Make copy of IsletGraph with threshold 2*thresh for mantle
      IsletGraph Icopy; Igraphcount++;if(verbosemem)cout << "Adding Icopy: IGraphcount = " << Igraphcount << endl;
      CopyIsletGraph(I,&Icopy);
      bool AddEdgesFailed = 0;
       if(types[t] == "ad")  AddEdgesFailed = AddEdges(&Icopy.allGraph,"bb",2*(*I).bbthresh,shadowthresh);
      else if(types[t] == "b")AddEdgesFailed = AddEdges(&Icopy.allGraph,"adm",2*(*I).admthresh,shadowthresh);
      //AddEdgesFailed = AddEdges(&Icopy.allGraph,"adb",2*(*I).adbthresh,shadowthresh);
      MakeGraphPlanarByRemovingLongerEdges(&Icopy.allGraph,1);
      //RemoveVerticesAwayFromSets(&Icopy.allGraph,types[t],2*(*I).adbthresh);
      Create_adAnd_bGraphs(&Icopy);
      //PrintoutGnuplot(&Icopy.allGraph,"Icopyall");
      //PrintoutGnuplot(&Icopy.bGraph,"Icopyb");
      //PrintoutGnuplot(&Icopy.adGraph,"Icopyadm");
      //if(types[t] == "b")PrintoutGnuplotIslet(I,fprefix);
      if(verbose) PrintoutGnuplotIslet(&Icopy,fprefix+".2xgr."+types[t]);
      if(types[t] == "ad")SetComponentsAndLabelExteriorPaths(&Icopy,"b",fprefix, "2xgr");
      if(types[t] == "b")SetComponentsAndLabelExteriorPaths(&Icopy,"ad",fprefix,"2xgr");
      if(verbose) cout << "Finished setting components" << endl;cout.flush();
      
      //string sIsletNum = static_cast<ostringstream*>( &(ostringstream() << (*I).IsletNum) )->str();
      if(verbose)PrintoutGnuplotExteriorPaths(&Icopy,fprefix+".2xgr.exteriorpaths_"+types[t]);
    
      FindAD_BetaSets(&Icopy,fprefix,"2xgr", types[t]);
      verbose=0;
      if(verbose)PrintoutGnuplotSetLoops(&Icopy,0 ,fprefix+".2xgr.setLoops."+types[t]);
      verbose=0;
      cout << "Running Partial Loops on " << types[t] << " components" << endl;
      cout.flush();

      vector<vector<vertex_t> > * Spcopy;
      vector<vector<vertex_t> > * Lpcopy;
      if(types[t] == "b"){
  	    Spcopy = &(Icopy.betasets);
  	    Lpcopy = &(Icopy.betasetLoops);
      }
      if(types[t] == "ad"){
  	    Spcopy = &(Icopy.adsets);
  	    Lpcopy = &(Icopy.adsetLoops);
      }
      
      if(verbose){
  	    //cout << "At beginning:Mp size is " << (*Mp).size() << endl;
  	    //for(int i = 0; i < (*Mp).size();i++)cout << " i = " << i << " and Mp[i].size is " << (*Mp)[i].size() << endl;
  	    cout << "There are " << (*Cp).size() << " components" << endl;cout.flush();
      }
      for(int c1r = 0; c1r < (*Cp).size();c1r++){
  	    if((*Mtp)[c1r] == 0){
  	      bool inset = 0;
  	      //check if c1r is in set with increased threshold (Icopy)
  	      for(int i = 0; i < (*Spcopy).size();i++){
  	        for(int j = 0; j < (*Spcopy)[i].size();j++){
  	          for(int ci = 0; ci < (*Cp)[c1r].size();ci++){
  	    	if((*Spcopy)[i][j] == (*Cp)[c1r][ci]){
  	    	  inset = 1;
  	    	  (*Mtp)[c1r] = 2;
  	    	  for(int cc = 0; cc < (*Ep1)[c1r].size(); cc++){
		        if((*Gp1)[(*Ep1)[c1r][cc]].IsletCellNum < 100000)(*MtPerc)[c1r].push_back(1.0);
		        else (*MtPerc)[c1r].push_back(0.0);
		      }
  	    	  //(*MtPerc)[c1r] = 1.0;
  	    	  // //determine number of components
  	    	  // vector<int> compnums;vectorcount++;if(verbosemem)cout << "Adding compnums: vectorcount = " << vectorcount << endl;
  	    	  // for(int lp = 0; lp < (*Lpcopy)[i].size();lp++){
		      //   // if((*Lpcopy)[i][lp] < 100000){
		      //   if(verbose) cout << "Lpcopy[" << lp << "] is " << (*Gp2)[(*Lpcopy)[i][lp]].IsletCellNum << " with compnum " << (*Gp2)[(*Lpcopy)[i][lp]].CompNum <<  endl; 
		      //   bool incomp = 0;
		      //   for(int cn = 0; cn < compnums.size();cn++){
		      //     if((*Gp2)[(*Lpcopy)[i][lp]].CompNum ==compnums[cn]){
		      // 	incomp = 1;
		      // 	cn = compnums.size();
		      //     }
  	    	      
  	    	  //     if(!incomp)compnums.push_back((*Gp2)[(*Lpcopy)[i][lp]].CompNum);
  	    	  //   }
  	    	    
  	    	  // }
  	    	  // (*Mtncomp)[c1r] = compnums.size();
  	    	  // if(verbose)cout << "The loop is made up of " << (*Mtncomp)[c1r] << " components "<< endl;
  	    	  // RemoveVector(&compnums);vectorcount--;if(verbosemem)cout << "Removing compnums: vectorcount = " << vectorcount << endl;
  	    	  //update path based on original graph
  	    	  int ncellsAdded = 0;
  	    	  if(verbose){
  	    	    cout << "Its loop before consists of " << endl;
  	    	    for(int m = 0; m < (*Lpcopy)[i].size();m++){
  	    	      //if((*Lpcopy)[i][m] < 100000)
  	    		    cout << (*Gp2)[(*Lpcopy)[i][m]].IsletCellNum << "(" << (*Lpcopy)[i][m] <<")  ";
		    	      //else
  	    		    //cout << (*Lpcopy)[ci][m] << " " ;
  	    	    }
  	    	    cout << endl << endl;
  	    	  }
  	    	  for(int k = 0; k < (*Lpcopy)[i].size()-1;k++){
  	    	    //check if edge between k and k+1 in Gp1
  	    	    if((*Lpcopy)[i][k] < 100000 && (*Lpcopy)[i][k+1] < 100000){
  	    	      tie(e,b) = edge((*Lpcopy)[i][k],(*Lpcopy)[i][k+1], (*Gp2));
  	    	      if(!b){
  	    		        if(verbose)cout << "There is no edge between " << (*Gp2)[(*Lpcopy)[i][k]].IsletCellNum << " and " << (*Gp2)[(*Lpcopy)[i][k+1]].IsletCellNum << " with compnums " << (*Gp2)[(*Lpcopy)[i][k]].CompNum << " and " <<(*Gp2)[(*Lpcopy)[i][k+1]].CompNum <<  endl;
  	    		        //check if in same component
  	    		        if((*Gp2)[(*Lpcopy)[i][k]].CompNum == (*Gp2)[(*Lpcopy)[i][k+1]].CompNum){
		    	          
  	    		          //New code to replace dijkstra
  	    		          vector<vertex_t> tpath;vectorcount++;if(verbosemem)cout << "Adding tpath: vectorcount = " << vectorcount << endl; cout.flush();
  	    		          bool nopath = ShortestWeightedPathBetweenTwoVertices(Gp2,(*Lpcopy)[i][k],(*Lpcopy)[i][k+1],&tpath);
  	    		          //check if members of path are already in loop: if so, do not add path
  	    		          bool inpath = 0;
  	    		          for(int tp = 1; tp < tpath.size()-1;tp++){
  	    		            for(int lp = 0; lp < (*Lpcopy)[i].size();lp++){
  	    		              if(tpath[tp] == (*Lpcopy)[i][lp]){
  	    		        	inpath = 1; 
  	    		        	tp = tpath.size();
  	    		        	break;
  	    		              }
  	    		            }
  	    		          }
  	    		          if(!inpath)(*Lpcopy)[i].insert((*Lpcopy)[i].begin()+k,tpath.begin()+1, tpath.end()-1);
  	    		          RemoveVector(&tpath);vectorcount--;if(verbosemem)cout << "Removing tpath: vectorcount = " << vectorcount << endl; cout.flush();
		    	          
  	    		        }
  	    	      }
  	    	    }
  	    	  }
  	    	  if(verbose){
  	    	    cout << "Its loop after consists of " << endl;
  	    	    for(int m = 0; m < (*Lpcopy)[i].size();m++){
  	    	      if((*Lpcopy)[i][m] < 100000) cout << (*Gp2)[(*Lpcopy)[i][m]].IsletCellNum << " ";
  	    	      else cout << (*Lpcopy)[i][m] << " " ;
  	    	    }
  	    	    cout << endl << endl;
  	    	  }
  	    	  (*Mp)[c1r].insert((*Mp)[c1r].begin(),(*Lpcopy)[i].begin(),(*Lpcopy)[i].end());
		      
		      
		      
  	    	  if(verbose){
  	    	    cout << "Component " << c1r << " is in a saturated set.  " << endl;
  	    	    cout << "Its Mantle consists of " << endl;
  	    	    for(int m = 0; m < (*Mp)[c1r].size();m++){
  	    	      if((*Mp)[c1r][m] < 100000) cout << (*Gp2)[(*Mp)[c1r][m]].IsletCellNum << " ";
  	    	      else cout << (*Mp)[c1r][m] << " " ;
		          
  	    	    }
  	    	    cout << endl << endl;
  	    	  }
  	    	  break;
  	    	}
  	          }
  	          if(inset)j = (*Spcopy)[i].size(); 
  	        }
  	        if(inset)i = (*Spcopy).size();
  	      }
  	    }
      }
      RemoveIsletGraph(&Icopy);Igraphcount--;if(verbosemem)cout << "Removing Icopy: IGraphcount = " << Igraphcount << endl;
     }
     
     //Check if all mantles have been found
     bool mantleNotFound = 0;
     for(int i = 0; i < (*Mp).size();i++){
       if((*Mp)[i].size() == 0){
	        mantleNotFound = 1;
	        break;
       }
     }
     if(mantleNotFound){
       //Make copy of allGraph and add adb edges 
       Graph Gpcopy;graphcount++;if(verbosemem)cout << "Adding graph Gpcopy: graphcount = " << graphcount << endl;cout.flush();
       CopyGraphVertices(&(*I).allGraph,&Gpcopy,"abd");
       bool AddEdgesFailed = AddEdges(&Gpcopy,"adb",(*I).adbthresh,shadowthresh);
       //CopyGraphEdges(&(*I).allGraph,&Gpcopy,"adb");
       //PrintoutGnuplot(&Gpcopy, "temp.abdedges");
       
       //Remove edges of Gpcopy that are not in exterior
       vector<vector<int> > RemoveEdges;vectorcount++;if(verbosemem)cout << "Adding RemoveEdges: vectorcount = " << vectorcount << endl; cout.flush();
       vector<int> InteriorVertices; vectorcount++;if(verbosemem)cout << "Adding InteriorVertices: vectorcount = " << vectorcount << endl; cout.flush();
       if(verbose)cout <<"Checking edges that are not in exterior from Gpcopy" << endl;
       for(ep = edges(Gpcopy); ep.first!=ep.second; ep.first++){
	        //Set up CompPathNums
	        for(int i = 0; i < 2;i++){
	          Gpcopy[*ep.first].adCompPathNums.push_back(99999);
	          Gpcopy[*ep.first].bCompPathNums.push_back(99999);
	        }
	        
	        vertex_t uall = source(*ep.first,Gpcopy);
	        vertex_t vall = target(*ep.first,Gpcopy);
	        vertex_t u = (*I).allGraph[uall].vother;
	        vertex_t v = (*I).allGraph[vall].vother;
	        
	        if(verbose)cout << "Checking edge (" << Gpcopy[uall].IsletCellNum << "," << Gpcopy[vall].IsletCellNum << ")" << endl;
	        
	        bool inInteriorVerts = 0;
	        for(int i = 0; i < InteriorVertices.size();i++){
	          if(InteriorVertices[i]==uall){
	            inInteriorVerts = 1;
	            break;
	          }
	          if(InteriorVertices[i]==vall){
	            inInteriorVerts = 1;
	            break;
	          }
	        }
	        if(!inInteriorVerts){
	          if(verbose)cout << "Not in InteriorVerts" << endl;
	          //Check if uall and vall are in exterior paths
	          vector<vertex_t>  * Ep1,*Ep2;
	          Graph * Gp1,*Gp2;
	          int c1 = 0, c2 = 0;
	          bool uInExt = 1, vInExt = 1;
	          if(Gpcopy[uall].type =="b" ){
	            Gp1 = &((*I).bGraph);
	            c1 = (*Gp1)[u].CompNum;
	            Ep1 = &((*I).bExteriorPaths[c1]);
	            Gp2 = &((*I).adGraph);
	            c2 = (*Gp2)[v].CompNum;
	            Ep2 = &((*I).adExteriorPaths[c2]);
	          }
             
	          else{
	            Gp1 = &((*I).adGraph);
	            c1 = (*Gp1)[u].CompNum;
	            Ep1 = &((*I).adExteriorPaths[c1]);
	            Gp2 = &((*I).bGraph);
	            c2= (*Gp2)[v].CompNum;
	            Ep2 = &((*I).bExteriorPaths[c2]);
	          }
	          if(verbose) cout << "Set pointers for uall" << endl;
	          //Check uall
	          vector<int> ep1;vectorcount++;if(verbosemem)cout << "Adding ep1: vectorcount = " << vectorcount << endl; cout.flush();
	          if(verbose)cout <<  " The size of ep1 is " << (*Ep1).size() << endl;cout.flush();
	          for(int i = 0; i < (*Ep1).size()-1;i++){
	            if(verbose)cout << "Checking ep1[i] " << (*Ep1)[i] << endl; cout.flush();
	            if((*Ep1)[i]<100000){
	              if(verbose)cout << "Checking " << (Gpcopy)[uall].IsletCellNum  << " with " << Gpcopy[(*Gp1)[(*Ep1)[i]].vother].IsletCellNum << endl;
	              cout.flush();
	              if(uall ==(*Gp1)[(*Ep1)[i]].vother)
	          ep1.push_back(i);
	            }
	          }
	          if(verbose){
	            cout << "External Path 1 consists of " ;cout.flush();
	            for(int e = 0; e < (*Ep1).size();e++){
	              if((*Ep1)[e]<100000) cout << (*Gp1)[(*Ep1)[e]].IsletCellNum << " ";
	            }
	            cout << endl;
	            cout << "Vertex " << Gpcopy[uall].IsletCellNum << " is at locations " ;cout.flush();
	            for(int e = 0; e < ep1.size();e++)
	              cout << ep1[e] << " ";cout.flush();
	            cout << endl;
	          }
	          if(ep1.size() ==0){
	            uInExt = 0;
	            if(verbose)cout << "removing edges from vertex" << (Gpcopy)[vall].IsletCellNum << endl;cout.flush();
	            InteriorVertices.push_back(uall);
	          }
	          //Check vall
	          vector<int> ep2;vectorcount++;if(verbosemem)cout << "Adding ep2: vectorcount = " << vectorcount << endl; cout.flush();
	          for(int i = 0; i < (*Ep2).size()-1;i++){
	            if((*Ep2)[i]<100000){
	              //cout << "Checking " << (Gpcopy)[*vp.first].IsletCellNum  << " with " << Gpcopy[(*Gp)[(*Ep)[i]].vother].IsletCellNum << endl;
	              cout.flush();
	              if(vall ==(*Gp2)[(*Ep2)[i]].vother) ep2.push_back(i);
	            }
	          }
	          if(verbose){
	            cout << "External Path 2 consists of " ;
	            for(int e = 0; e < (*Ep2).size();e++){
	              if((*Ep2)[e]<100000) cout << (*Gp2)[(*Ep2)[e]].IsletCellNum << " ";
	            }
	            cout << endl;
	            cout << "Vertex " << Gpcopy[vall].IsletCellNum << " is at locations " ;
	            for(int e = 0; e < ep2.size();e++)
	              cout << ep2[e] << " ";
	            cout << endl;
	          }
	          if(ep2.size() ==0){
	            vInExt = 0;
	            if(verbose)cout << "removing edges from vertex" << (Gpcopy)[vall].IsletCellNum << endl;
	            InteriorVertices.push_back(vall);
	          }
             
	          //Check if uall and vall are in each other's shadows
	          if(ep1.size() > 0 && ep2.size() > 0){
	            bool edgeexists = 0;
	            //loop thru ep1
	            for(int e1 = 0; e1 < ep1.size();e1++){
	              vector<double> shadowRange1;vectorcount++;if(verbosemem)cout << "Adding shadowRange1(1): vectorcount = " << vectorcount << endl; cout.flush();
	              string ty; if(Gpcopy[uall].type == "b")ty = "b"; else ty = "ad";
	              if(verbose)cout << "Before FindShadowRange with c1 = " << c1 << "ep1[e1] = " << ep1[e1] << " type = " << Gpcopy[uall].type  <<  endl;cout.flush();
	              FindShadowRange(I,c1,ep1[e1],ty,shadowthresh,&shadowRange1);
	              if(verbose){
	                  cout << "Shadow range for c1 vertex " << (*Gp1)[(*Ep1)[ep1[e1]]].IsletCellNum << " is ";cout.flush();
	                  for(int s = 0; s < shadowRange1.size();s=s+2)cout << "[" << shadowRange1[s] << "," << shadowRange1[s+1] << "]" ;
	                  cout << endl;cout.flush();
	              }
	              for(int e2 = 0; e2 < ep2.size();e2++){
	                  //check if e2 is in shadow range of e1
	                  double edgeAngle1 = GetAngleOfLine(Gpcopy[uall].x,Gpcopy[uall].y,Gpcopy[vall].x,Gpcopy[vall].y);
	                  if(verbose)cout << "edgeAngle1 is " << edgeAngle1 << endl;
	                  for(int s = 0; s < shadowRange1.size();s=s+2){
	                    if(edgeAngle1 >= shadowRange1[s] && edgeAngle1 <= shadowRange1[s+1]){
	                      if(verbose)cout << "  edgeangle1 is within shadow range" << endl;
	                      bool inRange = 0;
	                      vector<double> shadowRange2;vectorcount++;if(verbosemem)cout << "Adding shadowrange2(1):vectorcount = " << vectorcount << endl;
	                      if(Gpcopy[uall].type == "b")ty = "ad"; else ty = "b";
	                      FindShadowRange(I,c2,ep2[e2],ty,shadowthresh,&shadowRange2);
	                      if(verbose){
	                        cout << "Shadow range for c2 vertex " << (*Gp2)[(*Ep2)[e2]].IsletCellNum << " is ";
	                        for(int s = 0; s < shadowRange2.size();s=s+2)cout << "[" << shadowRange2[s] << "," << shadowRange2[s+1] << "]" ;
	                        cout << endl;
	                      }
	                      double edgeAngle2 = mod(edgeAngle1+PI,2*PI);
	                      if(verbose)cout << "Edge angle 2 is " << edgeAngle2 << endl;
	                      for(int s2 = 0; s2 < shadowRange2.size();s2=s2+2){
	                        if(edgeAngle2 >= shadowRange2[s2] && edgeAngle2 <= shadowRange2[s2+1]){
	         	                  inRange = 1;
	         	                  if(Gpcopy[uall].type == "b"){
	         	                    Gpcopy[*ep.first].bCompPathNums[0] = c1;
	         	                    Gpcopy[*ep.first].bCompPathNums[1] = ep1[e1];
	         	                    Gpcopy[*ep.first].adCompPathNums[0] = c2;
	         	                    Gpcopy[*ep.first].adCompPathNums[1] = ep2[e2];
	         	                  }
	         	                  else{
	         	                    Gpcopy[*ep.first].adCompPathNums[0] = c1;
	         	                    Gpcopy[*ep.first].adCompPathNums[1] = ep1[e1];
	         	                    Gpcopy[*ep.first].bCompPathNums[0] = c2;
	         	                    Gpcopy[*ep.first].bCompPathNums[1] = ep2[e2];
	         	                  }
	         	                  break;
	                        }
	                      }
	                      if(inRange){
	                        edgeexists = 1;
	                        s=shadowRange1.size();
	                      }
	                      RemoveVector(&shadowRange2);vectorcount--;if(verbosemem)cout << "Removing shadowRange2(1): vectorcount = " << vectorcount << endl; cout.flush();
	                    }
	                  }
	                  if(edgeexists)e2 = ep2.size();
	           
	              }
	              RemoveVector(&shadowRange1);vectorcount--;if(verbosemem)cout << "Removing shadowRange1(1): vectorcount = " << vectorcount << endl; cout.flush();
	              if(edgeexists)e1 = ep1.size();
	            }
	            if(!edgeexists){
	              vector<int> temp; vectorcount++;
	              temp.push_back(Gpcopy[uall].CellNum);
	              temp.push_back(Gpcopy[vall].CellNum);
	              RemoveEdges.push_back(temp);
	              RemoveVector(&temp);vectorcount--;
	            }
	          }
	          RemoveVector(&ep2);vectorcount--;if(verbosemem)cout << "Removing ep1: vectorcount = " << vectorcount << endl; cout.flush();
	          RemoveVector(&ep1);vectorcount--;if(verbosemem)cout << "Removing ep2: vectorcount = " << vectorcount << endl; cout.flush();
	        }
       }
 
       //Clear out interior vertices
       if(verbose)cout << "Clearing vertices " ;
       for(int i = 0; i < InteriorVertices.size();i++){
	        if(verbose)cout << Gpcopy[InteriorVertices[i]].IsletCellNum << " ";
	        for(oep = out_edges(InteriorVertices[i],Gpcopy);oep.first!=oep.second;oep.first++)ClearEdge(&Gpcopy,*oep.first);
	        clear_vertex(InteriorVertices[i],Gpcopy);//removes edges but leave vertex
       }
       if(verbose) cout << endl << endl;
       RemoveVector(&InteriorVertices);vectorcount--;if(verbosemem)cout << "Removing InteriorVertices: vectorcount = " << vectorcount << endl; cout.flush();
  
       //Remove edges from Gpcopy
       if(verbose)cout << "Removing edges " ;
       for(int i = 0; i < RemoveEdges.size();i++){
	        if(verbose)cout << "[" << Gpcopy[RemoveEdges[i][0]].IsletCellNum << "," << Gpcopy[RemoveEdges[i][1]].IsletCellNum << "] " ;
	        edge_t e; bool b;
	        tie(e,b) = edge(RemoveEdges[i][0], RemoveEdges[i][1],Gpcopy);
	        ClearEdge(&Gpcopy,e);
	        remove_edge(RemoveEdges[i][0], RemoveEdges[i][1],Gpcopy);
       }
       if(verbose) cout << endl << endl;
       RemoveVector2d(&RemoveEdges);vectorcount--;if(verbosemem)cout << "Removing RemoveEdges: vectorcount = " << vectorcount << endl; cout.flush();
       //PrintoutGnuplot(&Gpcopy, "temp.abdedges2");

       //Add mantles to unassigned sets
       for(int c1r = 0; c1r < (*Cp).size();c1r++){
	        if((*Mtp)[c1r]==0){
	          (*Mtp)[c1r] = 3;
	          vector<edge_t>  mantle;vectorcount++; if(verbosemem)cout << "Adding mantle: vectorcount = " << vectorcount << endl;cout.flush();
	          vector<vector<int> > VertexComponentPathNum;vectorcount++;if(verbosemem)cout << "Adding VertexComponentPathNum: vectorcount = " << vectorcount << endl; cout.flush();
	          //iterate through c1r external path
	          if(verbose){
	            cout << "Component " << c1r << " is not in a set.   " << endl;
	            cout << "Its external path is " ;cout.flush();
	            for(int ep1 = 0;ep1 < (*Ep1)[c1r].size();ep1++){
	              //if((*Ep1)[c1r][ep1]< 100000)
	              cout << (*Gp1)[(*Ep1)[c1r][ep1]].IsletCellNum << " ";
	            }
	            cout << endl;
	          }
	          for(int  ep1 = 0;ep1 < (*Ep1)[c1r].size()-1;ep1++){
	            if((*Ep1)[c1r][ep1] >99999)
	              (*MtPerc)[c1r].push_back(0);
	            else{
	              vertex_t u = (*Ep1)[c1r][ep1]; //associated with adGraph
	              vertex_t uall = (*Gp1)[(*Ep1)[c1r][ep1]].vother;
	              //check out_degree of uall wrt Gpcopy to see if there are any edges with mantle
	              if(out_degree(uall,Gpcopy) == 0){
	          (*MtPerc)[c1r].push_back(0);
	          vector<double> shadowRange_orig;vectorcount++;if(verbosemem)cout << "Adding shadowrange_orig:vectorcount = " << vectorcount << endl;cout.flush();
	          FindShadowRange(I,c1r,ep1,types[t],shadowthresh,&shadowRange_orig);
	          //if(verbose)PrintoutGnuplotPartialLoopsShadowRange(I,c1r,ep1,&shadowRange_orig,fprefix+"."+types[t],"./",types[t]);
	          //if(verbose)PrintoutGnuplotPartialLoopsShadowRange(I,c1r,ep1,&shadowRange_orig,fprefix+"."+types[t]+".left","./",types[t],"blue");
	          RemoveVector(&shadowRange_orig);vectorcount--;if(verbosemem)cout << "Removing shadowrange_orig:vectorcount = " << vectorcount << endl;cout.flush();
	          if(verbose)cout << "Vertex " << Gpcopy[uall].IsletCellNum << " has no neighbors in Gpcopy" << endl;
	              } 
	              else{
	          vector<double> shadowRange_orig;vectorcount++;if(verbosemem)cout << "Adding shadowrange_orig:vectorcount = " << vectorcount << endl;cout.flush();
	          vector<double> shadowRange_left;vectorcount++;if(verbosemem)cout << "Adding shadowrange_left:vectorcount = " << vectorcount << endl;cout.flush();
	          FindShadowRange(I,c1r,ep1,types[t],shadowthresh,&shadowRange_orig);
	          //copy shadowrange_orig to shadowrange_left
	          for(int s = 0; s < shadowRange_orig.size();s++)shadowRange_left.push_back(shadowRange_orig[s]);
	          if(verbose){
	            cout << "Shadow range for c1 vertex " << (*Gp1)[u].IsletCellNum << " is ";
	            for(int s = 0; s < shadowRange_orig.size();s=s+2)cout << "[" << shadowRange_orig[s]*180.0/PI << "," << shadowRange_orig[s+1]*180.0/PI << "]" ;
	            cout << endl;
	          }
	         
	          double origShadowLength = 0;
	          for(int s = 0; s < shadowRange_orig.size();s+=2)origShadowLength += shadowRange_orig[s+1]-shadowRange_orig[s];
	          if(verbose){
	            if(origShadowLength > 0){
	              //string sc1r = static_cast<ostringstream*>( &(ostringstream() << c1r) )->str();
	              //if(verbose)PrintoutGnuplotPartialLoopsShadowRange(I,c1r,ep1,&shadowRange_orig,fprefix+"."+types[t],"./",types[t]);
	            }
	          }
	          if(shadowRange_orig.size()>0){
	            //add edges to mantle and update shadow range
	            for(oep = out_edges(uall,Gpcopy);oep.first!=oep.second; oep.first++){
	              vertex_t vall = target(*oep.first,Gpcopy);
	              //add oep to mantle 
	              bool inmantle = 0;
	              for(int m = 0;m < mantle.size();m++){
	                if(source(mantle[m],Gpcopy) == vall && target(mantle[m],Gpcopy) == uall){
	         	 inmantle = 1;
	         	 break;
	                }
	              }
	              if(!inmantle){
	                mantle.push_back(*oep.first);
	                //Check if in VertexComponentPathnum 
	                bool inVCP = 0;
	               
	                for(int vcp = 0; vcp < VertexComponentPathNum.size();vcp++){
	         	 if(verbose)cout << "comparing (" << VertexComponentPathNum[vcp][0] << "," << VertexComponentPathNum[vcp][1] << "," << VertexComponentPathNum[vcp][2] << ") with (" << (*I).allGraph[vall].vother << "," << Gpcopy[*oep.first].adCompPathNums[0] << "," << Gpcopy[*oep.first].adCompPathNums[1] << ")" << endl;
	         	 if(types[t] == "b" && VertexComponentPathNum[vcp][0] == (*I).allGraph[vall].vother &&
	         	    VertexComponentPathNum[vcp][1] == Gpcopy[*oep.first].adCompPathNums[0] &&
	         	    VertexComponentPathNum[vcp][2] == Gpcopy[*oep.first].adCompPathNums[1] ){
	         	   inVCP = 1;
	         	   break;
	         	 }
	         	 else if(types[t] == "ad" && VertexComponentPathNum[vcp][0] == (*I).allGraph[vall].vother &&
	         		 VertexComponentPathNum[vcp][1] == Gpcopy[*oep.first].bCompPathNums[0] &&
	         		 VertexComponentPathNum[vcp][2] == Gpcopy[*oep.first].bCompPathNums[1] ){
	         	   inVCP = 1;
	         	   break;
	         	 }
	                }
	                if(!inVCP){
	         	 vector<int> temp;vectorcount++;
	         	 temp.push_back((*I).allGraph[vall].vother);
	         	 if(types[t] == "b"){
	         	   temp.push_back(Gpcopy[*oep.first].adCompPathNums[0]);
	         	   temp.push_back(Gpcopy[*oep.first].adCompPathNums[1]);
	         	 }
	         	 else{
	         	   temp.push_back(Gpcopy[*oep.first].bCompPathNums[0]);
	         	   temp.push_back(Gpcopy[*oep.first].bCompPathNums[1]);
	         	 }
	         	 VertexComponentPathNum.push_back(temp);
	         	 if(verbose)cout << "for edge [" << Gpcopy[source(*oep.first,Gpcopy)].IsletCellNum << "," << Gpcopy[target(*oep.first,Gpcopy)].IsletCellNum << " adding " << (*Gp2)[temp[0]].IsletCellNum << "," << temp[1] << "," << temp[2] << ") to VertexComponentPathNum" << endl;
	         	 RemoveVector(&temp);vectorcount--;
	                }
	              }
	              //update shadowrange_left
	              double shadthresh = 0;
	              if(types[t]=="b")shadthresh = (*I).admthresh / 2.0;
	              else shadthresh = (*I).bbthresh / 2.0;
	              if(verbose)cout << "Before find cell shadows" << endl;
	              cout.flush();
	              if(shadowRange_left.size() > 0)FindCellShadows(I,  uall, vall, &shadowRange_left,shadthresh);
	              if(verbose)cout << "After find cell shadows" << endl;cout.flush();
	            }
	          }
	          double leftShadowLength = 0;
	          for(int s = 0; s < shadowRange_left.size();s+=2)leftShadowLength += shadowRange_left[s+1]-shadowRange_left[s];
	          if(verbose){
	            //PrintoutGnuplotPartialLoopsShadowRange(I,c1r,ep1,&shadowRange_left,fprefix+"."+types[t]+".left","./",types[t],"blue");
	          }
	          if(origShadowLength > 0)(*MtPerc)[c1r].push_back(1.0-leftShadowLength/origShadowLength);
	          else (*MtPerc)[c1r].push_back(0);
	          if(verbose)cout << "The percentage of vertex " << (Gpcopy)[uall].IsletCellNum << " covered by a mantle is " << (*MtPerc)[c1r][ep1] << endl;
	          RemoveVector(&shadowRange_orig);vectorcount--;if(verbosemem)cout << "Removing shadowRange_orig: vectorcount = " << vectorcount << endl; cout.flush();
	          RemoveVector(&shadowRange_left);vectorcount--;if(verbosemem)cout << "Removing shadowRange_left: vectorcount = " << vectorcount << endl; cout.flush();
	              }
	            }
	          }
	          RemoveVector(&mantle);vectorcount--;if(verbosemem)cout << "Removing mantle: vectorcount = " << vectorcount << endl; cout.flush();
	          //Sort mantle by component and pathnum
	          sort(VertexComponentPathNum.begin(),VertexComponentPathNum.end(),compareByLastEntry);
	          sort(VertexComponentPathNum.begin(),VertexComponentPathNum.end(),compareBySecondEntry);
	          if(verbose){
	            cout << "The sorted components and pathnums in the mantle are:" << endl;
	            for(int m = 0; m < VertexComponentPathNum.size();m++)
	              cout << (*Gp2)[VertexComponentPathNum[m][0]].IsletCellNum << "\t" << VertexComponentPathNum[m][1] << "\t" << VertexComponentPathNum[m][2] << endl;
	            cout << endl;
	          }
	          // //count number of components in mantle
	          // int cnum = 0;
	          // if(VertexComponentPathNum.size() > 0){
	          //   VertexComponentPathNum[0][1];(*Mtncomp)[c1r]++;
	          //   for(int p = 0; p < VertexComponentPathNum.size();p++){
	          //     if(VertexComponentPathNum[p][1] !=cnum){
	          // 	(*Mtncomp)[c1r]++;
	          // 	cnum =VertexComponentPathNum[p][1];
	          //     }
	          //   }
	          // }
	          // if(verbose) cout << "There are " << (*Mtncomp)[c1r] << "components in the mantle " << endl;cout.flush();
	           
	          //upate ordering per component so that vertices are in a straight line
	          if(VertexComponentPathNum.size()>2){
	            int begpos = 0,pos = 0, endpos = 0, compnum = VertexComponentPathNum[0][1];
	            if(verbose){
	              cout << "The component being checked has external path " ;
	              for(int p = 0; p < (*Ep2)[compnum].size();p++){
	                  if((*Ep2)[compnum][p]<100000) cout << (*Gp2)[(*Ep2)[compnum][p]].IsletCellNum << " ";
	              }
	              cout << endl << endl;
	            }
	            while(pos < VertexComponentPathNum.size()){
	             
	              if(verbose)cout << "Inside while" << endl;cout.flush();
	              while(pos <VertexComponentPathNum.size() && VertexComponentPathNum[pos][1] == compnum)pos++;
	              if(pos <VertexComponentPathNum.size()-1)endpos = pos;
	              else endpos = pos-1;
	              if(verbose)cout << "At beginning: begpos = " << begpos << " pos = " << pos << " and endpos = " << endpos << endl;cout.flush();
	              bool noloop = 0;
	              bool edgeexists = 0;
	              tie(e,edgeexists) = edge(VertexComponentPathNum[begpos][0],VertexComponentPathNum[endpos][0], *Gp2);
	              if(edgeexists){//check if loop exists in mantle
	                  for(int p = begpos; p < endpos;p++){
	                    bool b = 0;
	                    tie(e,b) = edge(VertexComponentPathNum[p][0],VertexComponentPathNum[p+1][0], *Gp2);
	                    if(!b){ 
	                      noloop = 1;
	                      break;
	                    }
	                  }
	              }
	              while(edgeexists && endpos - begpos > 1 && noloop){
	                 if(verbose)cout << " An edge exists between " << (*Gp2)[VertexComponentPathNum[begpos][0]].IsletCellNum << " and " << (*Gp2)[VertexComponentPathNum[endpos][0]].IsletCellNum << endl;cout.flush();
	                 //remove vertex from beginpos and add to pos
	                 vector<int> temp;vectorcount++; 
	                 for(int vv = 0; vv < 3;vv++)temp.push_back(VertexComponentPathNum[begpos][vv]);
	                 VertexComponentPathNum.erase(VertexComponentPathNum.begin()+begpos);
	                 VertexComponentPathNum.insert(VertexComponentPathNum.begin()+endpos,temp);
	                 RemoveVector(&temp);vectorcount--;
	                 tie(e,edgeexists) = edge(VertexComponentPathNum[begpos][0],VertexComponentPathNum[endpos][0], *Gp2);
	              }
	              begpos = pos;
	              if(pos <VertexComponentPathNum.size())compnum = VertexComponentPathNum[pos][1];
	              if(verbose)cout << "After iteration: begpos = " << begpos << " pos = " << pos << " and endpos = " << endpos << " and compnum = " << compnum <<  endl;cout.flush();
	            }
	          }
	          if(verbose){
	            cout << "The sorted lines in the mantle are:" << endl;cout.flush();
	            for(int m = 0; m < VertexComponentPathNum.size();m++)
	              cout << (*Gp2)[VertexComponentPathNum[m][0]].IsletCellNum << "\t" << VertexComponentPathNum[m][1] << "\t" << VertexComponentPathNum[m][2] << endl;cout.flush();
	            cout << endl;cout.flush();
	          }
	          for(int m = 0; m < VertexComponentPathNum.size();m++)(*Mp)[c1r].push_back(VertexComponentPathNum[m][0]);
	          RemoveVector2d(&VertexComponentPathNum);vectorcount--;if(verbosemem)cout << "Removing VertexComponentPathNum: vectorcount = " << vectorcount << endl; cout.flush();
	        }
       }
       if(verbose)cout << "finished setting mantle" << endl;cout.flush();
       RemoveGraph(&Gpcopy);graphcount--;if(verbosemem)cout << "Removing graph Gpcopy:graphcount = " << graphcount << endl; cout.flush();
     }
 
     if(verbose){
       cout << "MantleSet consists of " << endl;cout.flush();
       for(int m = 0; m < (*Mp).size();m++){
	        for(int mm = 0; mm < (*Mp)[m].size();mm++){
	          if((*Mp)[m][mm] < 100000)
	            cout << (*Gp2)[(*Mp)[m][mm]].IsletCellNum << " ";
	          else cout << (*Mp)[m][mm] << " ";
	        }
	        cout << endl;
       }
     }
  }
  RemoveVector(&types);vectorcount--;if(verbosemem)cout << "Removing types: vectorcount = " << vectorcount << endl; cout.flush();
  
  cout << "Before removing IsletGraph" << endl;cout.flush();
  
  if(vectorcount != 0 || graphcount != 0 || Igraphcount != 0) cout << "Problem with vectors and/or graphs ! vectorcount = " << vectorcount << " and graphcount = " << graphcount <<" and Igraphcount = " << Igraphcount <<" in FindPartialLoops" << endl;
  cout << "Finished FindPartialLoops" <<endl;cout.flush();
  return;
}
        
//*********************************************************************************************//
//*********************************************************************************************//
//*********************************************************************************************//

void FindShadowRange(IsletGraph * I,  int c, int pathnum, string type, double shadowthresh, vector<double> * shadowRange){

  vector<vector<vertex_t> > * Cp;
  vector<vector<vertex_t> > * Ep;
  Graph * Gp;
  bool verbose=0,verbosemem=0;

  if(type == "b"){
    Cp = &((*I).bComponents);
    Ep = &((*I).bExteriorPaths);
    Gp = &((*I).bGraph);
  }
  if(type == "ad"){
    Cp = &((*I).adComponents);
    Ep = &((*I).adExteriorPaths);
    Gp = &((*I).adGraph);
  }
  vertex_t u = (*Ep)[c][pathnum];	     
  if(verbose)cout << "Checking pathnum = " << pathnum  << ". The component size is " << (*Cp)[c].size() << " and out degree of u is " << out_degree(u,(*Gp)) << endl;
  //if component is an isolated cell
  if((*Cp)[c].size() == 1){
    if(verbose)cout << "Cell " << (*Gp)[u].IsletCellNum << " is isolated "<< endl;
    (*shadowRange).push_back(0); (*shadowRange).push_back(2*PI);
  }
  
  //if cell only has 1 neighbor
  else if(out_degree(u,(*Gp)) == 1){//pendant end: subtract out vertex edge shadow from 2PI
    if(verbose)cout << "Cell " << (*Gp)[u].IsletCellNum << " is a pendant" << endl;
    double vx = 0, vy = 0;
    //    if((*Ep)[c][pathnum+1] <100000){
      vertex_t v = (*Ep)[c][pathnum+1]; //associated with bGraph
      vx = (*Gp)[v].x;
      vy = (*Gp)[v].y;
      //}
    double edgeAngle = GetAngleOfLine((*Gp)[u].x,(*Gp)[u].y,vx,vy);
    double distance = Distancexy((*Gp)[u].x,(*Gp)[u].y,vx,vy);
    double shadowAngle = tan(shadowthresh/distance);
    double shadowbegin = mod(edgeAngle-shadowAngle, 2*PI);
    double shadowend = mod(edgeAngle+shadowAngle, 2*PI);
    if (verbose) cout << " The edge angle is " << edgeAngle/(2*PI)*360 << " and the shadow angle is " << shadowAngle << ". The edge's shadow begins and ends at " << shadowbegin/(2*PI)*360 << " and " << shadowend/(2*PI)*360 << endl;
    if(shadowbegin < shadowend){
      (*shadowRange).push_back(0); (*shadowRange).push_back(shadowbegin);
      (*shadowRange).push_back(shadowend); (*shadowRange).push_back(2*PI);
    }
    else{
      (*shadowRange).push_back(shadowend); (*shadowRange).push_back(shadowbegin);
    }
  }
  //cell is part of path
 
  else{
    if(verbose)cout << "Cell " << (*Gp)[u].IsletCellNum << " is a part of a path" << endl; cout.flush();
    //Find previous cell in path
    double prevvx = 0, prevvy = 0;
    int prevpos = 0;
    if(pathnum == 0) prevpos = (*Ep)[c].size()-2;
    else prevpos = pathnum-1;
    //    if((*Ep)[c][prevpos] < 100000){
      if(verbose)cout << "  Previous vertex is " << (*Gp)[(*Ep)[c][prevpos]].IsletCellNum << endl;cout.flush();
      vertex_t prevv;
      prevv=(*Ep)[c][prevpos];
      prevvx = (*Gp)[prevv].x; prevvy = (*Gp)[prevv].y;
      
      //}
    
    //Find next cell in path
    double nextvx = 0, nextvy = 0;
    //if((*Ep)[c][pathnum+1] < 100000){
      vertex_t nextv = (*Ep)[c][pathnum+1];
      nextvx = (*Gp)[nextv].x; nextvy = (*Gp)[nextv].y;
      if(verbose)cout << "  Next vertex is " << (*Gp)[(*Ep)[c][pathnum+1]].IsletCellNum << endl;
      //}
    
    double edgeAngle1 = GetAngleOfLine((*Gp)[u].x,(*Gp)[u].y,prevvx,prevvy);
    double distance1 = Distancexy((*Gp)[u].x,(*Gp)[u].y,prevvx,prevvy);
    double edgeAngle2 = GetAngleOfLine((*Gp)[u].x,(*Gp)[u].y,nextvx,nextvy);
    double distance2 = Distancexy((*Gp)[u].x,(*Gp)[u].y,nextvx,nextvy);
    double shadowbegin = mod(edgeAngle1+tan(shadowthresh/distance1),2*PI);
    double shadowend = mod(edgeAngle2-tan(shadowthresh/distance2),2*PI);
    if (verbose) cout << " Edge angle 1 is " << edgeAngle1/(2*PI)*360 << " and " << endl 
		      << " Edge angle 2 is " << edgeAngle2/(2*PI)*360 << endl
		      << " the edge's shadow begins and ends at " << shadowbegin/(2*PI)*360 << " and " << shadowend/(2*PI)*360 << endl;
    if(edgeAngle1 < edgeAngle2){
      if(shadowbegin< shadowend){
	if(verbose) cout << "Case 1: not checked" <<endl;
	(*shadowRange).push_back(shadowbegin); (*shadowRange).push_back(shadowend);
      }
    }
    else{
      if(shadowbegin > edgeAngle1){
	if(shadowend < edgeAngle2 ){ //case 1 
	  if(verbose) cout << "Case 2: not checked" <<endl;
	  (*shadowRange).push_back(0); (*shadowRange).push_back(shadowend);
	  (*shadowRange).push_back(shadowbegin); (*shadowRange).push_back(2*PI);
	}
	else if(shadowbegin < shadowend){
	  if(verbose) cout << "Case 3: not checked" <<endl;
	  (*shadowRange).push_back(shadowbegin); (*shadowRange).push_back(shadowend);
	}
      }
      else if(shadowend> shadowbegin){
	if(verbose) cout << "Case 4: not checked" <<endl;
	(*shadowRange).push_back(shadowbegin); (*shadowRange).push_back(shadowend);
      }
    }
  }
  if(verbose){
    cout << "The shadowrange is " << endl;
    cout.flush();
    for(int s = 0; s < (*shadowRange).size();s=s+2)
      cout << "[" << (*shadowRange)[s] << "," << (*shadowRange)[s+1] << "] ";
    cout << endl;
    cout.flush();
  }  
  //no vectors, graphs, or paths, or cells
  return;
}

//*********************************************************************************************//
//*********************************************************************************************//
//*********************************************************************************************//


void MakeGraphPlanarByRemovingLongerEdges(Graph * Gp,bool opposite){
  std::pair<edge_iter, edge_iter> ep1,ep2;
  vertex_t u1, v1, u2, v2, u;
  edge_t e;
  bool verbose=0,verbosemem=0;

  if(verbose || verbosemem)cout << "*****Inside MakeGraphPlanarByRemovingLongerEdges" << endl;
  if(verbose && opposite)cout << "using opposites" << endl;
  
  bool intersections = 0;
  double Ix = 0,Iy = 0;
  edge_t ei, ej;
  for(ep1 = edges((*Gp));ep1.first!=ep1.second;ep1.first++){
    bool oppositepassed = 0;
    ep2 = edges((*Gp));
    while(ep2.first!=ep1.first)ep2.first++;
    ep2.first++;
    for(; ep2.first!=ep2.second;ep2.first++){
      oppositepassed = 0;
      u1 = source(*ep1.first, (*Gp));
      v1 = target(*ep1.first, (*Gp));
      u2 = source(*ep2.first, (*Gp));
      v2 = target(*ep2.first, (*Gp));
      
      //Check if edge is of opposite type
      int bcell1 = 0,adcell1 = 0;
      if((*Gp)[u1].type == "b")bcell1++;
      else adcell1++;
      if((*Gp)[v1].type == "b")bcell1++;
      else adcell1++;
      int bcell2=0,adcell2=0;
      if((*Gp)[u2].type == "b")bcell2++;
      else adcell2++;
      if((*Gp)[v2].type == "b")bcell2++;
      else adcell2++;
      if((bcell1==2 &&adcell2==2)|| (adcell1==2 && bcell2==2))oppositepassed =1;
      if(!opposite || (opposite && oppositepassed)){
	if(verbose)cout << "inside opposite check1: u1,v1 are of type "  << (*Gp)[u1].type <<  "," << (*Gp)[v1].type<< " and u2,v2 is of type " << (*Gp)[u2].type <<  "," << (*Gp)[v2].type<<endl;
	bool intersect = 0;
	verbose=0;
	intersect = IntersectionPointBetweenLines((*Gp)[u1].x, (*Gp)[u1].y,
						  (*Gp)[v1].x, (*Gp)[v1].y,
						  (*Gp)[u2].x, (*Gp)[u2].y,
						  (*Gp)[v2].x, (*Gp)[v2].y,
						  &Ix, &Iy,verbose);
	verbose=0;
	if(verbose)cout << "intersect = " << intersect <<  endl;
	if(intersect){
	  verbose=0;
	  if(verbose)cout << "There is an intersection between (" << (*Gp)[u1].IsletCellNum << "," << (*Gp)[v1].IsletCellNum << ") and (" << (*Gp)[u2].IsletCellNum << "," << (*Gp)[v2].IsletCellNum << ") with IP = (" << Ix << "," << Iy << endl;
	  verbose=0;
	  intersections = 1;
	  ei = *ep1.first; ej = *ep2.first;
	  while(ep1.first!= ep1.second)ep1.first++;
	  ep1.first--;
	  break;
	}
      }
    }
  }
  while(intersections){
    verbose=0;
    if(verbose)cout << "inside intersections" << endl;
    //Remove longer edge
    double dist1=Distancexy((*Gp)[u1].x,(*Gp)[u1].y,(*Gp)[v1].x,(*Gp)[v1].y);
    double dist2=Distancexy((*Gp)[u2].x,(*Gp)[u2].y,(*Gp)[v2].x,(*Gp)[v2].y);
    
    if(verbose)cout << "dist1 = " << dist1 << " and dist2 = " << dist2 << endl;
    if (dist1>=dist2){
      verbose=0;
      if(verbose)cout << "There is an intersection between (" << (*Gp)[u1].IsletCellNum << "," << (*Gp)[v1].IsletCellNum << ") and (" << (*Gp)[u2].IsletCellNum << "," << (*Gp)[v2].IsletCellNum << ") with IP = (" << Ix << "," << Iy << endl;
      if(verbosemem)cout << "Removing edge " << (*Gp)[u1].IsletCellNum<< "," << (*Gp)[v1].IsletCellNum  << ")" << endl;
      ClearEdge(Gp,ei);
      remove_edge(ei,(*Gp));
      verbose=0;
    }
    else{
      verbose=0;
      if(verbose)cout << "There is an intersection between (" << (*Gp)[u1].IsletCellNum << "," << (*Gp)[v1].IsletCellNum << ") and (" << (*Gp)[u2].IsletCellNum << "," << (*Gp)[v2].IsletCellNum << ") with IP = (" << Ix << "," << Iy << endl;
      if(verbosemem)cout << "Removing edge " << (*Gp)[u2].IsletCellNum<< "," << (*Gp)[v2].IsletCellNum  << ")" << endl;
      ClearEdge(Gp,ej);
      remove_edge(ej,(*Gp));
      verbose=0;
    }
    
    //Check if any other intersections
    if(verbose)cout << "There are " << num_edges((*Gp)) << " edges left" << endl;
    if(num_edges((*Gp)) == 1)intersections =0;
    else{
      for(ep1 = edges((*Gp));ep1.first!=ep1.second;ep1.first++){
	bool oppositepassed = 0;
	ep2 = edges((*Gp));
	while(ep2.first!=ep1.first)ep2.first++;
	ep2.first++;
	for(; ep2.first!=ep2.second;ep2.first++){
	  oppositepassed = 0;
	  u2 = source(*ep2.first, (*Gp));
	  v2 = target(*ep2.first, (*Gp));
	  if(verbose)cout << "Checking for new intersections" << endl;
	  int bcell1 = 0,adcell1 = 0;
	  u1 = source(*ep1.first, (*Gp));
	  v1 = target(*ep1.first, (*Gp));
	  if((*Gp)[u1].type == "b")bcell1++;
	  else adcell1++;
	  if((*Gp)[v1].type == "b")bcell1++;
	  else adcell1++;
	  if(verbose) cout << "For edge (" << (*Gp)[u1].IsletCellNum << "," << (*Gp)[v1].IsletCellNum << ") bcell1 = " << bcell1 << " and adcell1 = " << adcell1 << endl;
	  bool intersect = 0;
	  intersections = 0;
	  int bcell2=0,adcell2=0;
	  if((*Gp)[u2].type == "b")bcell2++;
	  else adcell2++;
	  if((*Gp)[v2].type == "b")bcell2++;
	  else adcell2++;
	  if(verbose) cout << "  For edge2 (" << (*Gp)[u2].IsletCellNum << "," << (*Gp)[v2].IsletCellNum << ") bcell2 = " << bcell2 << " and adcell2 = " << adcell2 << endl;
	  if((bcell1==2 &&adcell2==2)|| (adcell1==2 && bcell2==2))oppositepassed=1;
	  if(!opposite || (opposite && oppositepassed)){
	    if(verbose)cout << "inside opposite check: u1,u2 are of type "  << (*Gp)[u1].type <<  "," << (*Gp)[u2].type<< " and v1,v2 is of type " << (*Gp)[v1].type <<  "," << (*Gp)[v2].type<<endl;
	    verbose=0;
	    intersect = IntersectionPointBetweenLines((*Gp)[u1].x, (*Gp)[u1].y,
						      (*Gp)[v1].x, (*Gp)[v1].y,
						      (*Gp)[u2].x, (*Gp)[u2].y,
						      (*Gp)[v2].x, (*Gp)[v2].y,
						      &Ix, &Iy,verbose);
	    verbose=0;
	    if(verbose)cout<< "    intersect = " << intersect << endl;
	    if(intersect){
	      verbose=0;
	      if(verbose)cout << "There is an intersection between (" << (*Gp)[u1].IsletCellNum << "," << (*Gp)[v1].IsletCellNum << ") and (" << (*Gp)[u2].IsletCellNum << "," << (*Gp)[v2].IsletCellNum << ") with IP = (" << Ix << "," << Iy << ")" << endl;
	      verbose=0;
	      intersections = 1;
	      ei = *ep1.first; ej = *ep2.first;
	      while(ep1.first!= ep1.second)ep1.first++;
	      ep1.first--;
	      break;
	    }
	  }
	}
      }
    }
  }
  
  if(verbose || verbosemem)cout << "*****Finished MakeGraphPlanarByRemovingLongerEdge" << endl;
  cout.flush();
  //no vectors
  return;
}
//*********************************************************************************************//
//*********************************************************************************************//
//*********************************************************************************************//

double MinAfter2ndHighestPeak(vector<vector<double> > * vec, vector<vector<double> > * Peaks){
  int max1pos = 0,max2pos = 0;
  double maxpeak1 = 0, maxpeak2 = 0;
  double minx = 0;
  bool verbose=0;
  
  
  if(verbose){
    cout << "The peaks after cleaning are at" << endl;
    for(int i = 0; i < (*Peaks).size();i++)
      cout << "(" << (*Peaks)[i][0] << "," << (*Peaks)[i][1] << "), ";
    cout << endl;
  }
  //Find first 2 highest peaks
  if((*Peaks).size() > 2){
    for(int i = 0; i < (*Peaks).size();i++){
      if((*Peaks)[i][0] > 30)break;
      if((*Peaks)[i][1] > maxpeak1){
	maxpeak2 = maxpeak1;
	max2pos = max1pos;
	maxpeak1 = (*Peaks)[i][1];
	max1pos = i;
      }
      else if((*Peaks)[i][1] > maxpeak2){
	maxpeak2 = (*Peaks)[i][1];
	max2pos = i;
      }
      
    }
    //find min between peak 1 and 2
    int peak1xi = 0;
    while((*vec)[peak1xi][0] < (*Peaks)[0][0])peak1xi++;
    int peak2xi = peak1xi;
    while((*vec)[peak2xi][0] < (*Peaks)[1][0])peak2xi++;
    double min0x = 0, min0y = 100000;
    for(int i = peak1xi; i < peak2xi; i++){
      if((*vec)[i][1] < min0y){
	min0y = (*vec)[i][1];
	min0x = (*vec)[i][0];
      }
    }
    if(verbose)cout << "min between peaks1 and 2 is at (" << min0x << "," << min0y << ")" << endl;

      
    int cutoffPeakpos = max2pos;
    double cutoffPeak = maxpeak2;
    if(max2pos < max1pos){
      cutoffPeakpos = max1pos;
      cutoffPeak = maxpeak1;
    }
    if(verbose)cout << "cutoff peak is at (" << (*Peaks)[cutoffPeakpos][0] << "," << cutoffPeak << ") with pos " << cutoffPeakpos << " and Peaks.size() is " << (*Peaks).size() <<  endl; 
    //find cutoff peak in vec
    if(cutoffPeakpos == (*Peaks).size()-1){
      minx = (*Peaks)[cutoffPeakpos][0]+2;
    }
    else if((*Peaks).size() > 2){
      //find peak2 in original vector
      int peak2xi = 0; 
      while((*vec)[peak2xi][0] < (*Peaks)[cutoffPeakpos][0])peak2xi++;
      if(verbose)cout << "peak2 is at (" << (*vec)[peak2xi][0] << "," << (*vec)[peak2xi][1] << ")" << endl;
      //find peak3 in original vector
      int peak3xi = peak2xi;
      while ((*vec)[peak3xi][0] < (*Peaks)[cutoffPeakpos+1][0])peak3xi++;
      if(verbose)cout << "peak3 is at (" << (*vec)[peak3xi][0] << "," << (*vec)[peak3xi][1] << ")" << endl;
      //find min between peak 2 and 3
      double min1y = 1000000;
      for(int i = peak2xi; i < peak3xi;i++){
	if((*vec)[i][1] < min1y){
	  min1y = (*vec)[i][1];
	  minx = (*vec)[i][0];
	}
      }
      if(verbose) cout << "the min between 2 and 3 is " << min1y << " at " << minx << endl;

      if((*Peaks).size() > cutoffPeakpos+2){
	if ((*Peaks)[cutoffPeakpos+2][0]<40){
	  //find peak 4
	  int peak4xi = peak3xi;
	  while ((*vec)[peak4xi][0] < (*Peaks)[cutoffPeakpos+2][0])peak4xi++;
	  if(verbose)cout << "peak4 is at (" << (*vec)[peak4xi][0] << "," << (*vec)[peak4xi][1] << ")" << endl;
	  //find min between peaks 3 and 4
	  double min2y = 100000;
	  double min2x = 100000;
	  for(int i = peak3xi;i<peak4xi;i++){
	    if((*vec)[i][1] < min2y ){
	      min2y = (*vec)[i][1];
	      min2x = (*vec)[i][0];
	    }
	  }
	  if(verbose) cout << "the min between 3 and 4 is " << min2y << " at " << min2x << endl;
	  //if min between peaks 3 and 4 is less than min between peaks 2 and 3, update location of min
	  if(min2y < min1y )minx = min2x;
	
	}
      }
    }
    else if((*Peaks).size() > 0){
      int peak1xi = 0; 
      while((*vec)[peak1xi][0] < (*Peaks)[cutoffPeakpos][0])peak1xi++;
      int min3xi = peak1xi;
      while ((*vec)[min3xi][1] > 1.0)min3xi++;
      minx = min3xi;
    }
  }
  else minx = 0;
  return minx;
}




// //*********************************************************************************************//
// //*********************************************************************************************//
// //*********************************************************************************************//



void PeakFinder(vector<vector<double> > * vec, vector<vector<double> >* Peaks){
  bool verbose=0,verbosemem=0;
  int vectorcount = 0;
  //find all peaks
  for(int i = 1; i < (*vec).size()-1;i++){
    if((*vec)[i][1]> (*vec)[i-1][1] && (*vec)[i][1]> (*vec)[i+1][1] && (*vec)[i][1] > 0.5 && (*vec)[i][0]> 5)
      (*Peaks).push_back((*vec)[i]);
  }
  if(verbose){
    cout << "The peaks before cleaning are at" << endl;
    for(int i = 0; i < (*Peaks).size();i++)
      cout << "(" << (*Peaks)[i][0] << "," << (*Peaks)[i][1] << "), ";
    cout << endl;
  }
  
  //find peaks to remove: if peaks are less than 3units apart remove lesser-valued peak 
  sort((*Peaks).begin(),(*Peaks).end(),compareByLastEntryBiggestFirst);
  for(int i = 0; i < (*Peaks).size();i++){
    for(int j = i+1;j < (*Peaks).size();j++){
      if(fabs((*Peaks)[i][0]-(*Peaks)[j][0])<3){
	(*Peaks)[j].clear();vector<double>().swap((*Peaks)[j]);
	(*Peaks).erase((*Peaks).begin()+j);
	j--;
      }
    }
  }
  sort((*Peaks).begin(),(*Peaks).end(),compareByFirstEntry);

  //clean peaks of local plateaus: compare peak with min between nearest peaks
  vector<int>PeaksToRemove;vectorcount++;
  for(int i =(*Peaks).size()-1; i > 0; i--){
    //find minpeak
    int minpeak = i;
    if((*Peaks)[i-1][1] < (*Peaks)[minpeak][1])minpeak = i-1;
   
    //find min between peak i and peak i-1
    int peakim1x = 0;
    while((*vec)[peakim1x][0] < (*Peaks)[i-1][0])peakim1x++;
    int peakix = peakim1x;
    while((*vec)[peakix][0] < (*Peaks)[i][0])peakix++;
    double minval=100000, xminval = 0;
    for(int j = peakim1x; j < peakix;j++){
      if((*vec)[j][1] < minval){
	minval = (*vec)[j][1];
	xminval = (*vec)[j][0];
      }
    }
    if(verbose) cout << "For peak " << "(" << (*Peaks)[i-1][0] << "," << (*Peaks)[i-1][1] << ") and peak (" << (*Peaks)[i][0] << "," << (*Peaks)[i][1] << ") the minval is " << minval << " at pos " << xminval << endl;
    // cout << "minpeak is (" << (*Peaks)[minpeak][0] << "," << (*Peaks)[minpeak][1] << ")" << endl;
    if(((*Peaks)[minpeak][1] -minval < 0.05) ){//testing smallest between peaks i and i-1
      //check if minpeak is in Peaks to remove
      bool inRemove = 0;
      for(int p = 0; p < PeaksToRemove.size();p++){
	if(PeaksToRemove[p] == minpeak){
	  inRemove = 1;
	  break;
	}
      }
      if(!inRemove){
	if (verbose)cout << "Removing peak (" << (*Peaks)[minpeak][0] << "," << (*Peaks)[minpeak][1] << ")" <<endl; 
	PeaksToRemove.push_back(minpeak);
      }
    }
  }
      
  //Remove peaks
  sort(PeaksToRemove.begin(),PeaksToRemove.end());
  reverse(PeaksToRemove.begin(),PeaksToRemove.end());
  for(int i = 0;i < PeaksToRemove.size();i++){
    if(verbose)cout << "Removing peak (" << (*Peaks)[PeaksToRemove[i]][0] << "," << (*Peaks)[PeaksToRemove[i]][1] << ") " << endl;
    RemoveVector(&((*Peaks)[PeaksToRemove[i]]));
    (*Peaks).erase((*Peaks).begin()+PeaksToRemove[i]);
  }
  RemoveVector(&PeaksToRemove);vectorcount--;

  
  if(verbose){
    cout << "The peaks after cleaning are at" << endl;
    for(int i = 0; i < (*Peaks).size();i++)
      cout << "(" << (*Peaks)[i][0] << "," << (*Peaks)[i][1] << "), ";
    cout << endl;
  }
  if(vectorcount!=0)cout << "Problem in Peakfinder():vectorcount = " << vectorcount << endl;
  return;
}

//*********************************************************************************************//
//*********************************************************************************************//
//*********************************************************************************************//


bool setGrThresholds(vector<IsletGraph>* Islets, bool runGr, bool printgr, string gr_directory, string fig_directory, string fileprefix){

  int vectorcount = 0;
  bool verbose=0,verbosemem=0;
  if(verbosemem || verbose)cout << "*****Inside setGrThresholds" << endl;cout.flush();
  //See if Thresholds file exists
  //get first part of file header
  //string group;
  //group = getGroup(fprefix);
  //string ssubjno = static_cast<ostringstream*>( &(ostringstream() << subjno) )->str();
  //string tfile = gr_directory + "/"+group+".large."+ssubjno+".smooth2.thresh.minBetPeak2And3.dat";
  int subjno = 0;
  string ssubjno,tfile;

  tfile = gr_directory + "/" + fileprefix + ".smooth2.thresh.minBetPeak2And3.dat";

  //*****************************************************************************//
  //             Try to read in existing Islet threshold values                  //
  //*****************************************************************************//
  if(!runGr){
    bool addedAllThresh = 0;

    int ino = 0;

    while (!addedAllThresh){
      while(ino != (*Islets).size() && (*Islets)[ino].admthresh !=0  )ino++;
      if(verbose)cout << "ino = " << ino << " and (*Islets).size is " << (*Islets).size() << endl;
      if(ino == (*Islets).size())addedAllThresh = 1;
      else{

      	//subjno = (*Islets)[ino].IsletNum/10000;
      	//ssubjno = static_cast<ostringstream*>( &(ostringstream() << subjno) )->str();

      	if(file_exists(tfile)){
	        //cout << "Adding g(r) thresholds for subject " << group << " " << ssubjno << endl;
      	  if(verbose) cout << "Islet threshold file is " << tfile << endl; cout.flush();

      	  ifstream tin(tfile.c_str());

      	  while(tin.good()){
      	    int isletnum = 0; 
      	    double admt = 0.0, bbt = 0.0, adbt = 0.0;
      	    tin >> isletnum >> admt >> bbt >> adbt;
      	    if(tin.eof())break;
      	    //Set thresh
      	    isletnum += subjno*10000;
      	    for(int i = 0; i < (*Islets).size();i++){
      	      if((*Islets)[i].IsletNum == isletnum){
      	      	if(verbose)cout << "for islet " << (*Islets)[i].IsletNum << " adding threshes " << admt << " " << bbt << " " << adbt << endl;
      	       	(*Islets)[i].admthresh = admt;
		
      	       	(*Islets)[i].bbthresh = bbt;
      	       	(*Islets)[i].adbthresh = adbt;
		i = (*Islets).size();
      	       }
      	    }
      	  }
      	}
      	if(!file_exists(tfile) || (*Islets)[ino].admthresh == 0){
      	  if(!file_exists(tfile))cout << "...G(r) file " << tfile << " not found! Rerun G(r)!!!!" << endl;
      	  if((*Islets)[ino].admthresh == 0)cout << "The islet wasn't located in the file " << tfile << ". Rerun G(r)!!! "<< endl; cout.flush();
      	  return 1;
      	}
      }
    }	  
  }
  //*****************************************************************************//
  //                          Calculate G(r) thresholds                          //
  //*****************************************************************************//
  else{
    //set thresholds to 0
    for(int i = 0; i < (*Islets).size(); i++){
      (*Islets)[i].admthresh = 0;
      (*Islets)[i].bbthresh = 0;
      (*Islets)[i].adbthresh = 0;
    }
    //cout << endl;
    //cout << "...G(r) file " << tfile << " not found!  Calculating now (this could take a while!)..." ; cout.flush();
    
    ofstream throut;
    if(verbose)cout<< "created file " << tfile << endl;cout.flush();
     
    if(verbose)cout << "Determining threshold for adm..." << endl;
    DetermineThreshold(Islets, "adm", printgr, fig_directory, fileprefix);//printgr specifies whether to graph g(r)

    if(verbose)cout << "Determining threshold for bb..." << endl;
    DetermineThreshold(Islets, "bb", printgr, fig_directory, fileprefix);//printgr specifies whether to graph g(r)

    if(verbose)cout << "Determining threshold for adm..." << endl;
    DetermineThreshold(Islets, "adb", printgr, fig_directory, fileprefix);//printgr specifies whether to graph g(r)

    //cout << "Done." << endl;

    //Printout Thresholds
    
    //int old_subjno = 0;
    
    for(int i = 0; i < (*Islets).size(); i++){

      throut.open(tfile.c_str());
      throut << (*Islets)[i].IsletNum << "\t" << (*Islets)[i].admthresh << "\t" <<(*Islets)[i].bbthresh << "\t" << (*Islets)[i].adbthresh << endl; cout.flush();

      throut.close();
    }

  }
  
  if(vectorcount != 0)
    cout << "Problem with vectors! vectorcount = " << vectorcount << " in setGrThresholds() " << endl;
  if(verbosemem)cout << "*****Finished setGrThresholds: vectorcount = " << vectorcount << endl;
  return 0;
  //cout << "Finished\n" ;

}

//*********************************************************************************************//
//*********************************************************************************************//
//*********************************************************************************************//

bool ShortestWeightedPathBetweenTwoVertices(Graph * Gp, vertex_t u, vertex_t v, vector<vertex_t> * path){
  bool verbose=0,verbosemem=0;
  int vectorcount = 0;
  if(verbosemem)cout << "*****Inside ShortestWeightedPathBetweenTwoVertices" << endl;
  //check if edge exists between u and v
  edge_t e; bool b;
  tie(e,b) = edge(u,v,*Gp);
  if(b){
    if(verbose)cout << "An edge exists between " << (*Gp)[u].IsletCellNum << " and " << (*Gp)[v].IsletCellNum << ": finished" << endl;cout.flush();
    (*path).push_back(u);
    (*path).push_back(v);
    if(verbose)cout << "vectorcount = " << vectorcount << endl;
    return 0;
  }
 
  if(verbose)cout << "Inside ShortestWeightedPathBetweenTwoVertices...u = " << (*Gp)[u].IsletCellNum << " and " << (*Gp)[v].IsletCellNum <<  endl;cout.flush();
  //string su = static_cast<ostringstream*>( &(ostringstream() <<(*Gp)[u].IsletCellNum ) )->str();
  //string sv = static_cast<ostringstream*>( &(ostringstream() <<(*Gp)[v].IsletCellNum ) )->str();
  //if(verbose)PrintoutGnuplot(Gp, "temp.Shortest"+su+"_"+sv);

  vector<vertexListWeight> que1, que2;vectorcount+=2;if(verbosemem)cout << "Adding que1, qu2: vectorcount = " << vectorcount << endl;cout.flush();
  //Setup que1
  vertexListWeight temp;
  temp.vertex = u; temp.Parent = -1;temp.Weight = 0;
  que1.push_back(temp);
  std::pair<out_edge_iter, out_edge_iter> oep;
  bool foundPathIn1 = 0, foundPathIn2 = 0;
  vector<VertexDistance> Neighbors;vectorcount++;if(verbosemem)cout << "Adding Neighbors: vectorcount = " << vectorcount << endl;cout.flush();
  for(oep=out_edges(u,*Gp); oep.first!=oep.second; oep.first++){
    VertexDistance u1;
    u1.vertex = target(*oep.first,*Gp);
    u1.distance = (*Gp)[*oep.first].distance;
    Neighbors.push_back(u1);
  }
  sort(Neighbors.begin(), Neighbors.end(),compareByDistance);
  for(int i = 0; i < Neighbors.size();i++){
    temp.vertex = Neighbors[i].vertex; temp.Parent = 0;temp.Weight =Neighbors[i].distance; 
    que1.push_back(temp);
  }
  Neighbors.clear();
  
  //Setup que2
  
  temp.vertex = v;temp.Parent = -1;temp.Weight = 0; 
  que2.push_back(temp);
  for(oep=out_edges(v,*Gp); oep.first!=oep.second; oep.first++){
    VertexDistance v1;
    v1.vertex = target(*oep.first,*Gp);
    v1.distance = (*Gp)[*oep.first].distance;
    Neighbors.push_back(v1);
  }
  sort(Neighbors.begin(), Neighbors.end(),compareByDistance);
  for(int i = 0; i < Neighbors.size();i++){
    temp.vertex = Neighbors[i].vertex; temp.Parent = 0;temp.Weight =Neighbors[i].distance; 
    que2.push_back(temp);
  }
  Neighbors.clear();
  
  if(verbose){
    cout << "que1 consists of " << endl;
    for(int i = 0; i < que1.size(); i++)
      cout << (*Gp)[que1[i].vertex].IsletCellNum << "\t" << que1[i].Parent <<"\t" << que1[i].Weight <<  endl;
    cout.flush();
    cout << "que2 consists of " << endl;
    for(int i = 0; i < que2.size(); i++)
      cout << (*Gp)[que2[i].vertex].IsletCellNum << "\t" << que2[i].Parent <<"\t" << que2[i].Weight << endl;
    cout.flush();
  }
  int i1 = 1, i2 = 1;
  while(i1 < que1.size() && i2 < que2.size() && que1[i1].vertex != v && que2[i2].vertex != u){
    if(verbose)cout << "*** i1 = " << i1 << " with vertex " << (*Gp)[que1[i1].vertex].IsletCellNum << " and i2 = " << i2 << " with vertex " << (*Gp)[que2[i2].vertex].IsletCellNum <<" ***" << endl;cout.flush();
    
    for(oep=out_edges(que1[i1].vertex,*Gp); oep.first!=oep.second; oep.first++){
      VertexDistance u1;
      u1.vertex = target(*oep.first,*Gp);
      u1.distance = (*Gp)[*oep.first].distance;
      Neighbors.push_back(u1);
      if(verbosemem)cout << "Adding " << (*Gp)[u1.vertex].IsletCellNum << " to Neighbors " << endl;cout.flush();
    }
    sort(Neighbors.begin(), Neighbors.end(),compareByDistance);
    for(int i = 0; i < Neighbors.size();i++){
      //Check if neighbor is in que1
      bool inque1 = 0;
      for(int j = 0; j < que1.size();j++){
	        if(que1[j].vertex == Neighbors[i].vertex){
	          //check if Neighbor distance is greater than que distance
	          if(que1[j].Weight <Neighbors[i].distance+que1[i1].Weight)
	            inque1 = 1;  
	          else que1.erase(que1.begin()+j);
	          break;
	        }
      }
      if(!inque1){
	        temp.vertex = Neighbors[i].vertex;
	        temp.Parent = i1;
	        temp.Weight =Neighbors[i].distance+que1[i1].Weight;
	        if(temp.Weight > que1[que1.size()-1].Weight)que1.push_back(temp);
	        else{
	          int j = i1;
	          while(temp.Weight > que1[j].Weight)j++;
	         if(verbose)cout << "before insert" << endl;cout.flush();
	          que1.insert(que1.begin()+j,temp);
	         if(verbose)cout << "after insert" << endl;cout.flush();
	        }
      }
    }
    if(verbose){
      cout << "que1 consists of " << endl;
      for(int i = 0; i < que1.size(); i++)
	cout << (*Gp)[que1[i].vertex].IsletCellNum << "\t" << que1[i].Parent <<"\t" << que1[i].Weight <<  endl;
      cout.flush();
    }
    Neighbors.clear();
    
    //que2
     for(oep=out_edges(que2[i2].vertex,*Gp); oep.first!=oep.second; oep.first++){
      VertexDistance v1;
      v1.vertex = target(*oep.first,*Gp);
      v1.distance = (*Gp)[*oep.first].distance;
      Neighbors.push_back(v1);
      if(verbosemem)cout << "Adding " << (*Gp)[v1.vertex].IsletCellNum << " to Neighbors " << endl;cout.flush();
    }
    sort(Neighbors.begin(), Neighbors.end(),compareByDistance);
    for(int i = 0; i < Neighbors.size();i++){
      //Check if neighbor is in que2
      bool inque2 = 0;
      for(int j = 0; j < que2.size();j++){
	        if(que2[j].vertex == Neighbors[i].vertex){
	          //check if Neighbor distance is greater than que distance
	          if(que2[j].Weight <Neighbors[i].distance+que2[i2].Weight)
	            inque2 = 1;  
	          else que2.erase(que2.begin()+j);
	          break;
	        }
      }
      if(!inque2){
	        temp.vertex = Neighbors[i].vertex;
	        temp.Parent = i2;
	        temp.Weight =Neighbors[i].distance+que2[i2].Weight;
	        if(temp.Weight > que2[que2.size()-1].Weight)que2.push_back(temp);
	        else{
	          int j = i2;
	          while(temp.Weight > que2[j].Weight)j++;
	         if(verbose)cout << "before insert" << endl;cout.flush();
	          que2.insert(que2.begin()+j,temp);
	         if(verbose)cout << "after insert" << endl;cout.flush();
	        }
      }
    }
    if(verbose){
      cout << "que2 consists of " << endl;
      for(int i = 0; i < que2.size(); i++)
	cout << (*Gp)[que2[i].vertex].IsletCellNum << "\t" << que2[i].Parent <<"\t" << que2[i].Weight <<  endl;
      cout.flush();
    }
    Neighbors.clear();
    
    
    i1++;i2++;
    
    if(verbose){
      cout << "que1 consists of " << endl;
      for(int i = 0; i < que1.size(); i++)
	cout << (*Gp)[que1[i].vertex].IsletCellNum << "\t" << que1[i].Parent <<"\t" << que1[i].Weight <<  endl;
      cout.flush();
      cout << "que2 consists of " << endl;
      for(int i = 0; i < que2.size(); i++)
	cout << (*Gp)[que2[i].vertex].IsletCellNum << "\t" << que2[i].Parent <<"\t" << que2[i].Weight << endl;
      cout.flush();
    }
  }
  RemoveVector(&Neighbors);vectorcount--;if(verbosemem)cout << "Removing Neighbors: vectorcount = " << vectorcount << endl;cout.flush();

  if(que1.size() == i1 || que2.size() == i2){
    if(verbose && que1.size() == i1)cout << "que1.size = i1: no path" << endl;
    if(verbose && que2.size() == i2)cout << "que2.size = i2: no path" << endl;
    RemoveVector(&que1);vectorcount--;if(verbosemem)cout << "Removing que1: vectorcount = " << vectorcount << endl;cout.flush();
    RemoveVector(&que2);vectorcount--;if(verbosemem)cout << "Removing que2: vectorcount = " << vectorcount << endl;cout.flush();
    if(vectorcount!=0)cout << "Problem with vectorcount! vectorcount = " << vectorcount << endl;
    return 1;
  }
  else if( que1[i1].vertex == v){
    if(verbose){
      cout << "que1 consists of " << endl;
      for(int i = 0; i < que1.size(); i++)
	cout << (*Gp)[que1[i].vertex].IsletCellNum << "\t" << que1[i].Parent << endl;
    }
    vector<vertex_t> temp2;vectorcount++;if(verbosemem)cout << "Adding temp2: vectorcount = " << vectorcount << endl;cout.flush();
    int i = i1;
    while(i !=-1){
      temp2.push_back(que1[i].vertex);
      if(verbosemem)cout << "Adding " << (*Gp)[que1[i].vertex].IsletCellNum << " to path " << endl;
      i = que1[i].Parent;
    }
    //Reverse temp2
    for(int j = temp2.size()-1;j > -1; j--)
      (*path).push_back(temp2[j]);
    RemoveVector(&temp2);vectorcount--;if(verbosemem)cout << "Removing temp2: vectorcount = " << vectorcount << endl;cout.flush();
  }
  else{
    int i = i2;
    while(i !=-1){
      (*path).push_back(que2[i].vertex);
      i = que2[i].Parent;
    }
  }
  if(verbose){
    cout << "The path from vertex " << (*Gp)[u].IsletCellNum << " to " << (*Gp)[v].IsletCellNum << " is " << endl;cout.flush();
    for(int i = 0; i < (*path).size();i++)
      cout << (*Gp)[(*path)[i]].IsletCellNum << " ";
    cout << endl;
  }
  RemoveVector(&que1);vectorcount--;if(verbosemem)cout << "Removing que1: vectorcount = " << vectorcount << endl;cout.flush();
  RemoveVector(&que2);vectorcount--;if(verbosemem)cout << "Removing que2: vectorcount = " << vectorcount << endl;cout.flush();
  if(vectorcount!=0)cout << "Problem in ShortestWeightedPathBetweenTwoVertices with vectorcount! vectorcount = " << vectorcount << endl;
  
  return 0;
}

//*********************************************************************************************//
//*********************************************************************************************//
//*********************************************************************************************//

void SmoothVector5(vector<vector<double> > * vec, vector<vector<double> > * smvec){
  for(int j = 0; j < 2; j++){
    for(int i = 0; i < 2; i++)
      (*smvec)[i][j]=(*vec)[i][j];
    for(int i = (*vec).size()-2; i < (*vec).size(); i++)
      (*smvec)[i][j]=(*vec)[i][j];
  }
    
  for(int i = 2; i < (*vec).size()-2;i++){
    (*smvec)[i][0]=(*vec)[i][0];
    (*smvec)[i][1]=((*vec)[i-2][1]+(*vec)[i-1][1]+(*vec)[i][1]+(*vec)[i+1][1]+(*vec)[i+2][1])/5.0;
  }
  //no vectors
  return;
}
//*********************************************************************************************//
//*********************************************************************************************//
//*********************************************************************************************//

      
// void RotateIsletMaximalDistanceXaxis(IsletGraph *RI, string type){ 
 
//   int vectorcount = 0;
//   bool verbose=0,verbosemem=0;
//   if(verbose)PrintoutGnuplot(&(*RI).allGraph, "BeforeRotating");
//   //RI is originally a copy of I
//   //Determine furthest cells
//   double maxdist = 0;
//   int jmax = 0, kmax = 0;
//   for(int j = 0;j < num_vertices((*RI).allGraph);j++){
//     if((type == "bb" && (*RI).allGraph[j].type == "b") ||
//        (type == "adm" && ((*RI).allGraph[j].type == "a" || (*RI).allGraph[j].type == "d")) ||
//        type == "adb"){ 
//       double x1 = (*RI).allGraph[j].x, y1 = (*RI).allGraph[j].y;
//       for(int k = j+1; k < num_vertices((*RI).allGraph);k++){
// 	if((type == "bb" && (*RI).allGraph[k].type == "b") ||
// 	   (type == "adm" && ((*RI).allGraph[k].type == "a" || (*RI).allGraph[k].type == "d")) ||
// 	   type == "adb"){
// 	  double x2 = (*RI).allGraph[k].x, y2 = (*RI).allGraph[k].y;
// 	  double radius = Distancexy(x1,y1,x2,y2);
// 	  if(radius > maxdist){
// 	    maxdist = radius;
// 	    jmax = j;
// 	    kmax = k;
// 	  }
// 	}
//       }
//     }
//   }
  
//   //Make p1(x1,y1) the first point from left to right
//   double x1 = ((*RI).allGraph)[jmax].x, y1 = ((*RI).allGraph)[jmax].y;
//   double x2 = ((*RI).allGraph)[kmax].x, y2 = ((*RI).allGraph)[kmax].y;
//   if(x1>x2){
//     double tempx = x1; x1 = x2; x2 = tempx;
//     double tempy = y1; y1 = y2; y2 = tempy;
//   }
//   //find angle
//   double angle = 0;
//   double dy =(y2-y1);
//   double dx = (x2-x1);
//   angle = atan2(dy,dx);

//   for(int i = 0; i < num_vertices((*RI).allGraph);i++){
//     double oldx = (*RI).allGraph[i].x, oldy =(*RI).allGraph[i].y;
//     (*RI).allGraph[i].x = (oldx - x1)*cos(angle) + (oldy-y1)*sin(angle);//see Friday 8/4/2017 notes
//     (*RI).allGraph[i].y = (oldy - y1)*cos(angle) - (oldx-x1)*sin(angle);
//   }
//   if(verbose)PrintoutGnuplot(&(*RI).allGraph, "AfterRotating");
//   return;
// }

  

// void RandomizeOppositeEdgesTest(IsletGraph * I, string fileprefix){
//   bool verbose=0;
//   int vectorcount = 0;
//   if(verbose)cout << "Inside RandomizeOppositeEdges" << endl;cout.flush();
//   if(verbose)cout << "There are " << num_vertices((*I).allGraph) << " vertices and " << num_edges((*I).allGraph) << " edges" << endl;
//   //vector<int>List_num_opposites,numCompsad,numLoopsad,numCompsb,numLoopsb;vectorcount+=5;
//   //vector<double>PercInLoopsad,PercMantleSurroundedad,PercInLoopsb,PercMantleSurroundedb;vectorcount+=4;
 
//   std::pair<edge_iter, edge_iter> ep;
//   int numiter = 500;
//   if(5*num_vertices((*I).allGraph) > numiter)numiter = 5*num_vertices((*I).allGraph);
//   for(int i = 0; i < numiter; i++){
//     int num_opposites = 0;
//     for(ep = edges((*I).allGraph);ep.first!=ep.second;ep.first++){
//       vertex_t v1 = source(*ep.first,(*I).allGraph);
//       vertex_t v2 = target(*ep.first,(*I).allGraph);
//       if(((*I).allGraph)[v1].type == "b"){
// 	if(((*I).allGraph)[v2].type == "a" || ((*I).allGraph)[v2].type == "d")num_opposites++;
//       }
//       else if(((*I).allGraph)[v2].type == "b")num_opposites++;
//     }
//     if(num_opposites == 0){
//       cout << "Graph has no opposites" << endl;
//       //RemoveVector(&List_num_opposites);vectorcount--;
//       return;
//     }
//     if(verbose)cout << " i =" << i << ": There are " << num_opposites << endl;cout.flush();
//     ofstream opp((fileprefix+".Num_opposite_type").c_str(), ios::app);
//     opp << num_opposites << endl;
//     opp.close();
//     // List_num_opposites.push_back(num_opposites);
//     //Pick a random edge
//     int redgenum = rand()%num_opposites;
//     if(verbose)cout << "the random edge is " << redgenum << endl;cout.flush();
//     int edgeiter = 0;
//     for(ep = edges((*I).allGraph);ep.first!=ep.second;ep.first++){
//       vertex_t v1 = source(*ep.first,(*I).allGraph);
//       vertex_t v2 = target(*ep.first,(*I).allGraph);
//       if(((*I).allGraph)[v1].type == "b"){
// 	if(((*I).allGraph)[v2].type == "a" || ((*I).allGraph)[v2].type == "d")
// 	  edgeiter++;
//       }
//       else if(((*I).allGraph)[v2].type == "b")edgeiter++;
//       if(redgenum < edgeiter) {
// 	//Switch cell types
// 	string temptype = ((*I).allGraph)[source(*ep.first,(*I).allGraph)].type;
// 	((*I).allGraph)[source(*ep.first,(*I).allGraph)].type = ((*I).allGraph)[target(*ep.first,(*I).allGraph)].type;
// 	((*I).allGraph)[target(*ep.first,(*I).allGraph)].type = temptype;
// 	Create_adAnd_bGraphs(I);
// 	ClearIsletGraphVectors(I);
// 	SetComponentsAndLabelExteriorPaths(I,"b",fileprefix);
// 	SetComponentsAndLabelExteriorPaths(I,"ad",fileprefix);
// 	//Printout numcomps
// 	opp.open((fileprefix+".Num_bcomponents").c_str(), ios::app);
// 	opp << (*I).bComponents.size() << endl;
// 	opp.close();
// 	opp.open((fileprefix+".Num_adcomponents").c_str(), ios::app);
// 	opp << (*I).adComponents.size() << endl;
// 	opp.close();
// 	//Printout num beta/ad sets
// 	FindAD_BetaSets(I,fileprefix);
// 	opp.open((fileprefix+".Num_betasets").c_str(), ios::app);
// 	opp << (*I).betasets.size() << endl;
// 	opp.close();
// 	opp.open((fileprefix+".Num_adsets").c_str(), ios::app);
// 	opp << (*I).adsets.size() << endl;
// 	opp.close();
// 	int numBetasinSets = 0;
// 	for(int s = 0; s < (*I).betasets.size();s++)
// 	  numBetasinSets += (*I).betasets[s].size();
// 	opp.open((fileprefix+".Perc_betaCellsInSets").c_str(), ios::app);
// 	opp << (double)numBetasinSets/(double)(*I).nbeta << endl;
// 	opp.close();
// 	int numadsinSets = 0;
// 	for(int s = 0; s < (*I).adsets.size();s++)
// 	  numadsinSets += (*I).adsets[s].size();
// 	opp.open((fileprefix+".Perc_adCellsInSets").c_str(), ios::app);
// 	opp << (double)numadsinSets/(double)((*I).nalpha+(*I).ndelta)  << endl;
// 	opp.close();
// 	//Printout Partial Loops
// 	FindPartialLoops(I,4,fileprefix);
// 	double totPerc = 0;
// 	for(int s = 0; s < (*I).bMantlePercent.size(); s++){
// 	  if((*I).bMantle_type[s] == 1 ||(*I).bMantle_type[s] == 2)totPerc+=1;
// 	  else{
// 	    double Percsum = 0;
// 	    for(int ss = 0; ss < (*I).bMantlePercent[s].size();ss++) Percsum +=(*I).bMantlePercent[s][ss];
// 	    totPerc+=Percsum/(double)(*I).bMantlePercent[s].size();
// 	  }
// 	}
// 	opp.open((fileprefix+".Perc_bMantleSurrounded").c_str(), ios::app);
// 	opp << totPerc/(double)(*I).bMantlePercent.size()  << endl;
// 	opp.close();
// 	totPerc = 0;
// 	for(int s = 0; s < (*I).adMantlePercent.size(); s++){
// 	  if((*I).adMantle_type[s] == 1 ||(*I).adMantle_type[s] == 2)totPerc+=1;
// 	  else{
// 	    double Percsum = 0;
// 	    for(int ss = 0; ss < (*I).adMantlePercent[s].size();ss++) Percsum +=(*I).adMantlePercent[s][ss];
// 	    totPerc+=Percsum/(double)(*I).adMantlePercent[s].size();
// 	  }
// 	}
// 	opp.open((fileprefix+".Perc_adMantleSurrounded").c_str(), ios::app);
// 	opp << totPerc/(double)(*I).adMantlePercent.size()  << endl;
// 	opp.close();
// 	// numCompsad.push_back((*I).adComponents.size());
// 	// numLoopsad.push_back((*I).adsets.size());
// 	// numCompsb.push_back((*I).bComponents.size());
// 	// numLoopsb.push_back((*I).betasets.size());
// 	// //Calculate number of cells in sets
// 	// int numBetasinSets = 0;
// 	// for(int s = 0; s < (*I).betasets.size();s++)
// 	//   numBetasinSets += (*I).betasets[s].size();
// 	// PercInLoopsb.push_back((double)numBetasinSets/(double)(*I).nbeta);
// 	// int numadsinSets = 0;
// 	// for(int s = 0; s < (*I).adsets.size();s++)
// 	//   numadsinSets += (*I).adsets[s].size();
// 	// PercInLoopsad.push_back((double)numadsinSets/(double)((*I).nalpha+(*I).ndelta));
// 	//Calculate average percentage of coverage
// 	// double totPerc = 0;
// 	// for(int s = 0; s < (*I).bMantlePercent.size(); s++){
// 	//   if((*I).bMantle_type[s] == 1 ||(*I).bMantle_type[s] == 2)totPerc+=1;
// 	//   else{
// 	//     double Percsum = 0;
// 	//     for(int ss = 0; ss < (*I).bMantlePercent[s].size();ss++) Percsum +=(*I).bMantlePercent[s][ss];
// 	//     totPerc+=Percsum/(double)(*I).bMantlePercent[s].size();
// 	//   }
// 	// }
// 	// PercMantleSurroundedb.push_back(totPerc/(double)(*I).bMantlePercent.size());
// 	// totPerc = 0;
// 	// for(int s = 0; s < (*I).adMantlePercent.size(); s++){
// 	//   if((*I).adMantle_type[s] == 1 ||(*I).adMantle_type[s] == 2)totPerc+=1;
// 	//   else{
// 	//     double Percsum = 0;
// 	//     for(int ss = 0; ss < (*I).adMantlePercent[s].size();ss++) Percsum +=(*I).adMantlePercent[s][ss];
// 	//     totPerc+=Percsum/(double)(*I).adMantlePercent[s].size();
// 	//   }
// 	// }
// 	// PercMantleSurroundedad.push_back(totPerc/(double)(*I).adMantlePercent.size());
// 	break;
//       }
//     }
//   }
//  //  //Printout Results
//  //  ofstream opp((fileprefix+".Num_opposite_type").c_str());
//  //  for(int i = 0; i < List_num_opposites.size();i++)
//  //    opp << List_num_opposites[i] << endl;
//  //  opp.close();
//  //  RemoveVector(&List_num_opposites);vectorcount--;
//  //  opp.open((fileprefix+".Num_bcomponents").c_str());
//  //  for(int i = 0; i < numCompsb.size();i++)
//  //    opp << numCompsb[i] << endl;
//  //  opp.close();
//  //  RemoveVector(&numCompsb);vectorcount--;
//  //  opp.open((fileprefix+".Num_adcomponents").c_str());
//  //  for(int i = 0; i < numCompsad.size();i++)
//  //    opp << numCompsad[i] << endl;
//  //  opp.close();
//  //  RemoveVector(&numCompsad);vectorcount--;
//  // opp.open((fileprefix+".Num_bsets").c_str());
//  //  for(int i = 0; i < numLoopsb.size();i++)
//  //    opp << numLoopsb[i] << endl;
//  //  opp.close();
//  //  RemoveVector(&numLoopsb);vectorcount--;
//  //  opp.open((fileprefix+".Num_adsets").c_str());
//  //  for(int i = 0; i < numLoopsad.size();i++)
//  //    opp << numLoopsad[i] << endl;
//  //  opp.close();
//  //  RemoveVector(&numLoopsad);vectorcount--;
//  //  opp.open((fileprefix+".Perc_bInSets").c_str());
//  //  for(int i = 0; i < PercInLoopsb.size();i++)
//  //    opp << PercInLoopsb[i] << endl;
//  //  opp.close();
//  //  RemoveVector(&PercInLoopsb);vectorcount--;
//  //  opp.open((fileprefix+".Perc_adInSets").c_str());
//  //  for(int i = 0; i < PercInLoopsad.size();i++)
//  //    opp << PercInLoopsad[i] << endl;
//  //  opp.close();
//  //  RemoveVector(&PercInLoopsad);vectorcount--;
//  // opp.open((fileprefix+".Perc_bMantleSurrounded").c_str());
//  //  for(int i = 0; i < PercMantleSurroundedb.size();i++)
//  //    opp << PercMantleSurroundedb[i] << endl;
//  //  opp.close();
//  //  RemoveVector(&PercMantleSurroundedb);vectorcount--;
//  //  opp.open((fileprefix+".Perc_adMantleSurrounded").c_str());
//  //  for(int i = 0; i < PercMantleSurroundedad.size();i++)
//  //    opp << PercMantleSurroundedad[i] << endl;
//  //  opp.close();
//  //  RemoveVector(&PercMantleSurroundedad);vectorcount--;

//   PrintoutGnuplot(&(*I).allGraph,fileprefix+".Randomized");

//   if(vectorcount!=0)cout << "Problem in RandomizeOppositeEdges! vectorcount = " << vectorcount << endl;
//   return;
// }
  
// void RandomizeOppositeEdges(Graph * G, string fileprefix){
//   bool verbose=0;
//   int vectorcount = 0;
//   if(verbose)cout << "Inside RandomizeOppositeEdges" << endl;cout.flush();
//   if(verbose)cout << "There are " << num_vertices(*G) << " vertices and " << num_edges(*G) << " edges" << endl;
//   vector<int>List_num_opposites;vectorcount++;
//   std::pair<edge_iter, edge_iter> ep;
//   int numiter = 50000;
//   //if(5*num_vertices(*G) > 500)numiter = 5*num_vertices(*G);
//   for(int i = 0; i < numiter; i++){
//     int num_opposites = 0;
//     for(ep = edges(*G);ep.first!=ep.second;ep.first++){
//       vertex_t v1 = source(*ep.first,*G);
//       vertex_t v2 = target(*ep.first,*G);
//       if((*G)[v1].type == "b"){
// 	if((*G)[v2].type == "a" || (*G)[v2].type == "d")num_opposites++;
//       }
//       else if((*G)[v2].type == "b")num_opposites++;
//     }
//     if(num_opposites == 0){
//       cout << "Graph has no opposites" << endl;
//       RemoveVector(&List_num_opposites);vectorcount--;
//       return;
//     }
//     if(verbose)cout << " i =" << i << ": There are " << num_opposites << endl;cout.flush();
//     List_num_opposites.push_back(num_opposites);
//     //Pick a random edge
//     int redgenum = rand()%num_opposites;
//     if(verbose)cout << "the random edge is " << redgenum << endl;cout.flush();
//     int edgeiter = 0;
//     for(ep = edges(*G);ep.first!=ep.second;ep.first++){
//       vertex_t v1 = source(*ep.first,*G);
//       vertex_t v2 = target(*ep.first,*G);
//       if((*G)[v1].type == "b"){
// 	if((*G)[v2].type == "a" || (*G)[v2].type == "d")
// 	  edgeiter++;
//       }
//       else if((*G)[v2].type == "b")edgeiter++;
//       if(redgenum < edgeiter) {
// 	//Switch cell types
// 	string temptype = (*G)[source(*ep.first,*G)].type;
// 	(*G)[source(*ep.first,*G)].type = (*G)[target(*ep.first,*G)].type;
// 	(*G)[target(*ep.first,*G)].type = temptype;
// 	break;
//       }
//     }
//   }
//   ofstream opp((fileprefix+".Num_opposite_type").c_str());
//   //ofstream opp(("sims/"+fileprefix+".Num_opposite_type").c_str());
//   for(int i = 0; i < List_num_opposites.size();i++)
//     opp << List_num_opposites[i] << endl;
//   opp.close();
//   RemoveVector(&List_num_opposites);vectorcount--;
//   //PrintoutGnuplot(G,fileprefix+".Randomized");

//   if(vectorcount!=0)cout << "Problem in RandomizeOppositeEdges! vectorcount = " << vectorcount << endl;
//   return;
// }
  
  

// int connected_components_DS(Graph *Gp, vector<int> *components){
//   int vectorcount = 0;
//   bool verbose=0;
//   //if(verbose)PrintoutGnuplot(Gp,"connectedcomponentsDS");

//   std::pair<vertex_iter, vertex_iter> vp;
//   vector<vertexListWeight> Assigned;vectorcount++; //contains vertex_t vertex, int Parent (which will be component number); and Weight (which will be distance);
//   int Lastcompnum = -1;
//   for(vp = vertices(*Gp);vp.first!=vp.second;vp.first++){
//     vertex_t u = *vp.first;
//     if(verbose)cout << "Finding component for " << (*Gp)[u].IsletCellNum <<endl;
//     //Calculate distance between u and Assigned vertices
//     for(int i = 0; i < Assigned.size();i++)
//       Assigned[i].Weight = Distancexy((*Gp)[u].x, (*Gp)[u].y, (*Gp)[Assigned[i].vertex].x,(*Gp)[Assigned[i].vertex].y);
//     //Sort Assigned by distance
//     sort(Assigned.begin(),Assigned.end(),compareByWeight);
//     if(verbose){
//       cout << "The sorted Assigned vertices are:" << endl;
//       for(int i = 0;i < Assigned.size();i++)
//     	cout << (*Gp)[Assigned[i].vertex].IsletCellNum << "\t" << Assigned[i].Parent << "\t" << Assigned[i].Weight << endl;
//     }
//     //see if path exists between u and Assigned
//     bool NoPathExists = 1;
//     for(int i = 0; i < Assigned.size();i++){
//       vector<vertex_t>path;vectorcount++;
//       NoPathExists = PathBetweenTwoVertices(Gp, u, Assigned[i].vertex, &path);
//       if(!NoPathExists){//Assign u to i's component
// 	if(verbose)cout << "A path exists between " << (*Gp)[u].IsletCellNum << " and " << (*Gp)[Assigned[i].vertex].IsletCellNum << endl;
// 	vertexListWeight vlw;
// 	vlw.vertex = u;
// 	vlw.Parent = Assigned[i].Parent;
// 	Assigned.push_back(vlw);
// 	RemoveVector(&path);vectorcount--;
// 	break;
//       }
//       RemoveVector(&path);vectorcount--;
//     }
//     if(NoPathExists){
//       vertexListWeight vlw;
//       vlw.vertex = u;
//       Lastcompnum++;
//       vlw.Parent = Lastcompnum;
//       Assigned.push_back(vlw);
//       if(verbose)cout << "Assigning " << (*Gp)[u].IsletCellNum << " to component " << Lastcompnum  << endl;
//     }
//   }
//   for(int i = 0; i < Assigned.size();i++)
//     (*components)[(*Gp)[Assigned[i].vertex].CellNum] = Assigned[i].Parent;
//   if(verbose){
//     sort(Assigned.begin(),Assigned.end(),compareByParent);
//     cout << "The components are "<< endl;
//     for(int i = 0; i < Assigned.size();i++)
//       cout << (*Gp)[Assigned[i].vertex].IsletCellNum << "\t" << Assigned[i].Parent << endl;
//   }
//   RemoveVector(&Assigned);vectorcount--;
//   cout << "Finished finding components" <<endl;
//   if(vectorcount>0)cout << "Problem with connected_components_DS! vectorcount = " << vectorcount << endl; 
//   return Lastcompnum+1;
// }
      
  




// bool PathBetweenTwoVertices(Graph * Gp, vertex_t u, vertex_t v, vector<vertex_t> * path){
//   //Determines if path between two vertices exists
//   //quicker than Shortest weighted path
//   //returns a path between the two vertices in *path

//   bool verbose=0,verbosemem=0;
//   int vectorcount = 0;
// if(verbosemem)cout << "*****Inside PathBetweenTwoVertices" << endl;
//   vector<vertexListWeight> que1, que2;vectorcount+=2;if(verbosemem)cout << "Adding que1,que2: vectorcount = " << vectorcount << endl;cout.flush();

  
//   if(verbose)cout << "Inside PathBetweenTwoVertices...u = " << (*Gp)[u].IsletCellNum << " and " << (*Gp)[v].IsletCellNum <<  endl;cout.flush();
//   string su = static_cast<ostringstream*>( &(ostringstream() <<(*Gp)[u].IsletCellNum ) )->str();
//   string sv = static_cast<ostringstream*>( &(ostringstream() <<(*Gp)[v].IsletCellNum ) )->str();
 
//   //if(verbose)PrintoutGnuplot(Gp, "temp.Path"+su+"_"+sv);
//   //Setup que1
//   vertexListWeight temp;
//   temp.vertex = u; temp.Parent = -1;
//   que1.push_back(temp);
//   std::pair<out_edge_iter, out_edge_iter> oep;
//   bool foundPathIn1 = 0, foundPathIn2 = 0;
//   for(oep=out_edges(u,*Gp); oep.first!=oep.second; oep.first++){
//     vertex_t u1 = target(*oep.first,*Gp);
//     temp.vertex = u1; temp.Parent = 0;
//     que1.push_back(temp);
//     if(u1 == v){
//       if(verbose)cout << "Found in path 1" << endl;
//       foundPathIn1 = 1;
//       break;
//     }
//   }
//   if(!foundPathIn1){
//     temp.vertex = v; temp.Parent = -1;
//     que2.push_back(temp);
//     std::pair<out_edge_iter, out_edge_iter> oep;
//     for(oep=out_edges(v,*Gp); oep.first!=oep.second; oep.first++){
//       vertex_t v1 = target(*oep.first,*Gp);
//       temp.vertex = v1; temp.Parent = 0;
//       que2.push_back(temp);
//       if(v1 == u){
// 	foundPathIn2 = 1;
// 	if(verbose)cout << "Found in path 2" << endl;
// 	break;
//       }
//     }
//   }
//   int i1 = 1, i2 = 1;
//   while(!foundPathIn1 && !foundPathIn2 && i1 < que1.size() && i2 < que2.size()){
//     if(verbose)cout << "*** i1 = " << i1 << " with vertex " << (*Gp)[que1[i1].vertex].IsletCellNum << " and i2 = " << i2 << " with vertex " << (*Gp)[que2[i2].vertex].IsletCellNum <<" ***" << endl;
//     for(oep=out_edges(que1[i1].vertex,*Gp); oep.first!=oep.second; oep.first++){
//       vertex_t u1 = target(*oep.first,*Gp);
//       //Check if u1 is in que1
//       bool inque1 = 0;
//       for(int i = 0; i < que1.size();i++){
// 	if(que1[i].vertex == u1){
// 	  inque1 = 1;
// 	  break;
// 	}
//       }
//       if(!inque1){
// 	temp.vertex = u1; temp.Parent = i1;
// 	que1.push_back(temp);
//       }
//       if(u1 == v){
// 	foundPathIn1 = 1;
// 	if(verbose)cout << "Found in path 1" << endl;
// 	break;
//       }
//     }
//     if(!foundPathIn1){
//       temp.vertex = v; temp.Parent = 0;
//       que2.push_back(temp);
//       std::pair<out_edge_iter, out_edge_iter> oep;
//       for(oep=out_edges(que2[i2].vertex,*Gp); oep.first!=oep.second; oep.first++){
// 	vertex_t v1 = target(*oep.first,*Gp);
// 	//Check if v1 is in que2
// 	bool inque2 = 0;
// 	for(int i = 0; i < que2.size();i++){
// 	  if(que2[i].vertex == v1){
// 	  inque2 = 1;
// 	  break;
// 	  }
// 	}
// 	if(!inque2){
// 	  temp.vertex = v1; temp.Parent = i2;
// 	  que2.push_back(temp);
// 	}
// 	if(v1 == u){
// 	  foundPathIn2 = 1;
// 	  if(verbose)cout << "Found in path 2" << endl;
// 	  break;
// 	}
//       }
//     }
//     if(!foundPathIn1 && !foundPathIn2){
//       i1++;i2++;
//     }
//   }
 
//   if(que1.size() == i1 || que2.size() == i2){
//     if(verbose && que1.size() == i1)cout << "que1.size = i1: no path" << endl;
//     else if(verbose && que2.size() == i2)cout << "que2.size = i2: no path" << endl;
//     RemoveVector(&que1);vectorcount--;if(verbosemem)cout << "Removing que1: vectorcount = " << vectorcount << endl;cout.flush();
//     RemoveVector(&que2);vectorcount--;if(verbosemem)cout << "Removingque2: vectorcount = " << vectorcount << endl;cout.flush();
//     if(vectorcount!=0)cout << "Problem with vectorcount! vectorcount = " << vectorcount << endl;
//     return 1;
//   }
//   else if(foundPathIn1){
//     if(verbose){
//       cout << "que1 consists of " << endl;
//       for(int i = 0; i < que1.size(); i++)
// 	cout << (*Gp)[que1[i].vertex].IsletCellNum << "\t" << que1[i].Parent << endl;
//       cout.flush();
//     }
//     vector<vertex_t> temp2;vectorcount++;if(verbosemem)cout << "Adding temp2: vectorcount = " << vectorcount << endl;cout.flush();
//     int i = que1.size()-1;
//     while(i !=-1){
//       temp2.push_back(que1[i].vertex);
//       if(verbosemem)cout << "Adding " << (*Gp)[que1[i].vertex].IsletCellNum << " to path " << endl;
//       i = que1[i].Parent;
//     }
//     //Reverse temp2
//     for(int j = temp2.size()-1;j > -1; j--)
//       (*path).push_back(temp2[j]);
//     RemoveVector(&temp2);vectorcount--;if(verbosemem)cout << "Removing temp2: vectorcount = " << vectorcount << endl;cout.flush();
//   }
//   else{
//     int i = que2.size()-1;
//     while(i !=-1){
//       (*path).push_back(que2[i].vertex);
//       i = que2[i].Parent;
//     }
//   }
//   if(verbose){
//     cout << "The path from vertex " << (*Gp)[u].IsletCellNum << " to " << (*Gp)[v].IsletCellNum << " is " << endl;cout.flush();
//     for(int i = 0; i < (*path).size();i++)
//       cout << (*Gp)[(*path)[i]].IsletCellNum << " ";cout.flush();
//     cout << endl;cout.flush();
//   }
//   RemoveVector(&que1);vectorcount--;if(verbosemem)cout << "Removing que1: vectorcount = " << vectorcount << endl;cout.flush();
//   RemoveVector(&que2);vectorcount--;if(verbosemem)cout << "Removing que2: vectorcount = " << vectorcount << endl;cout.flush();
//   if(vectorcount!=0)cout << "Problem in PathBetweenTwoVertices with vectorcount! vectorcount = " << vectorcount << endl;cout.flush();

//   return 0;
// }


// void CalcGofrUsingShadow(vector<IsletGraph>* I,vector<vector<double> > *Gofr, string celltype,int Gno){ 
//   std::pair<vertex_iter, vertex_iter>vp,vp1,vp2;
  
//   if(Gno != 99999)cout << "Calculating Gofr for Islet " << Gno << " out of " << (*I).size() << "...";
//   else cout << "Calculating Gofr ..." << endl;
//   int vectorcount = 0;
//   bool verbose=0,verbosemem=0;
//   double shadowthresh = 4.0;
//   //Calculate longest distance between cells
//   double maxdist = 0;
//   for(int i = 0; i < num_vertices((*I)[Gno].allGraph);i++){
//     for(int j = i+1; j < num_vertices((*I)[Gno].allGraph);j++){
//       double distance = Distancexy((*I)[Gno].allGraph[i].x, (*I)[Gno].allGraph[i].y,(*I)[Gno].allGraph[j].x, (*I)[Gno].allGraph[j].y);
//       if(distance > maxdist)maxdist = distance;
//     }
//   }
//   double dr = 0.5;
//   int nbins = maxdist/dr + 4.0;
  
//   vector<vector<double> > gk_norm((*I).size(),vector<double>(nbins,0.0));vectorcount++;if(verbosemem)cout << "Adding gk_norm: vectorcount = " << vectorcount << endl; cout.flush();
//   int numIslets = 0;
//   int Gbeg = 0, Gend = (*I).size();
//   if(Gno!=99999){
//     Gbeg = Gno;
//     Gend = Gno+1;
//   }
  
//   for(int i = Gbeg; i < Gend;i++){
//     if(verbose &&  i%100 == 0){
//       cout << i << " " ;
//       cout.flush();
//     }
    
//     //Determine the (#vertices)^2
//     int numvertsq = 0;
//     if (celltype == "abd")  numvertsq = (int)pow(num_vertices((*I)[i].allGraph),2.0);
//     else if(celltype == "aa")numvertsq = (int)pow((*I)[i].nalpha,2.0);
//     else if(celltype == "bb")numvertsq = (int)pow((*I)[i].nbeta,2.0);
//     else if(celltype == "dd")numvertsq = (int)pow((*I)[i].ndelta,2.0); 
//     else if(celltype == "adm")numvertsq = (int)pow((*I)[i].nalpha+(*I)[i].ndelta,2.0);
//     else if(celltype == "ab")numvertsq = (*I)[i].nalpha*(*I)[i].nbeta;
//     else if(celltype == "ad")numvertsq = (*I)[i].nalpha*(*I)[i].ndelta;
//     else if(celltype == "bd")numvertsq = (*I)[i].ndelta*(*I)[i].nbeta;
//     else if(celltype == "adb")numvertsq = ((*I)[i].nalpha+(*I)[i].ndelta)*(*I)[i].nbeta;
    
//     if((celltype == "abd" && num_vertices((*I)[i].allGraph) > 1)||
//        (celltype == "aa" && (*I)[i].nalpha > 1)||
//        (celltype == "bb" && (*I)[i].nbeta > 1)||
//        (celltype == "dd" && (*I)[i].ndelta > 1)||
//        (celltype == "adm" && (*I)[i].nalpha+(*I)[i].ndelta> 1)||
//        (celltype == "ab" && (*I)[i].nalpha > 0 && (*I)[i].nbeta > 0) ||
//        (celltype == "ad" && (*I)[i].nalpha > 0 && (*I)[i].ndelta > 0) ||
//        (celltype == "bd" && (*I)[i].nbeta > 0 && (*I)[i].ndelta > 0) ||
//        (celltype =="adb" && ((*I)[i].nalpha > 0 ||(*I)[i].ndelta >0) && (*I)[i].nbeta > 0)){
//       numIslets++;
      
//       //Determine bounding box for islet area
//       double IsletArea = 0;
//       double xlength = 0,ylength = 0;
//       double xmin = 100000, xmax = 0, ymin = 100000,ymax = 0;
//       for(vp = vertices((*I)[i].allGraph);vp.first!=vp.second;vp.first++){
// 	if(((*I)[i].allGraph)[*vp.first].x >xmax)xmax = ((*I)[i].allGraph)[*vp.first].x;
// 	if(((*I)[i].allGraph)[*vp.first].y >ymax)ymax = ((*I)[i].allGraph)[*vp.first].y;
// 	if(((*I)[i].allGraph)[*vp.first].x <xmin)xmin = ((*I)[i].allGraph)[*vp.first].x;
// 	if(((*I)[i].allGraph)[*vp.first].y <ymin)ymin = ((*I)[i].allGraph)[*vp.first].y;
//       }
//       xlength = xmax-xmin;
//       ylength = ymax-ymin;
//       IsletArea = xlength*ylength;
      
//       if(verbose)cout << "For Islet " << (*I)[i].IsletNum <<"\t" <<  xmin << "\t" << xmax << "\t" << xlength << "\t" <<  ymin << "\t" << ymax << "\t" << ylength << endl;
      
//       //calculate the number of possible edges
//       for(vp1=vertices((*I)[i].allGraph);vp1.first!=vp1.second;vp1.first++){
// 	//calculate distance to neigboring cells
// 	vector<VertexDistance> Dists;vectorcount++;if(verbosemem)cout << "Adding vector Dists: vectorcount = " << vectorcount << endl;cout.flush();
// 	for(vp2=vertices((*I)[i].allGraph);vp2.first!=vp2.second;vp2.first++){
// 	  if(*vp1.first != *vp2.first){
// 	    double x1 = ((*I)[i].allGraph)[*vp1.first].x, y1 = ((*I)[i].allGraph)[*vp1.first].y;
// 	    double x2 = ((*I)[i].allGraph)[*vp2.first].x, y2 = ((*I)[i].allGraph)[*vp2.first].y;
// 	    VertexDistance vd;
// 	    vd.vertex=*vp2.first;
// 	    vd.distance = Distancexy(x1,y1,x2,y2);
// 	    Dists.push_back(vd);
// 	  }
// 	}
// 	sort(Dists.begin(),Dists.end(),compareByDistance);
// 	//Remove intervals and add calcGofr
// 	bool intervalfull = 0;
// 	vector<vector<double> > interval; vectorcount++;if(verbosemem)cout << "Adding vector interval: vectorcount = " << vectorcount << endl;cout.flush();
// 	vector<double> temp(2,0);vectorcount++;if(verbosemem)cout << "Adding vector temp: vectorcount = " << vectorcount << endl;cout.flush();
// 	temp[0]=0;temp[1]=2*PI;
// 	interval.push_back(temp);
// 	int dit = 0;
// 	while(interval.size() > 0 && dit < Dists.size()){
// 	  //Determine if dit is able to be added
// 	  //Calculate Angles
// 	  float alpha2 = atan(shadowthresh/Dists[dit].distance);
// 	  float beta2 = 0.0;
// 	  double x1 = (*I)[i].allGraph[*vp1.first].x,y1 = (*I)[i].allGraph[*vp1.first].y;
// 	  double x2 = (*I)[i].allGraph[Dists[dit].vertex].x,y2 = (*I)[i].allGraph[Dists[dit].vertex].y;
// 	  if(x2>=x1)beta2 = atan((y2-y1)/(x2-x1));//Q1 and 4
// 	  else if(y2>y1)beta2 = atan((y2-y1)/(x2-x1))+ PI;//Q2
// 	  else beta2 = atan((y2-y1)/(x2-x1))- PI;//Q3
// 	  beta2 = mod(beta2,2*PI);
// 	  vector<double> beta2shad(2,0);vectorcount++;if(verbosemem)cout << "Adding vector beta2shad: vectorcount = " << vectorcount << endl;cout.flush();
// 	  beta2shad[0] = mod(beta2-alpha2,2*PI);beta2shad[1] = mod(beta2+alpha2,2*PI);
// 	  if(verbose)cout<< "Checking " << (*I)[i].allGraph[*vp1.first].IsletCellNum << " with " << (*I)[i].allGraph[Dists[dit].vertex].IsletCellNum << "\n";
// 	  if(verbose){
// 	    cout << "v2's shadow is (" << beta2shad[0] << "," <<  beta2shad[1]  << ")" << endl;
// 	    cout << "Intervals are "  ;
// 	    for(int sn = 0; sn < interval.size();sn++)
// 	      cout << "[" << interval[sn][0] << "," << interval[sn][1] << "]," ;
// 	  } 
// 	  bool betagood = 0;
// 	  if(beta2shad[0]<beta2shad[1]){//new segment doesn't wrap
// 	    //loop thru segments to see where beta2shad[0] falls
// 	    for(int sn = 0;sn< interval.size(); sn++){
// 	      if(beta2shad[0]> interval[sn][0] && beta2shad[0] < interval[sn][1]){// beta2shad0 is located in this segment
// 		if(beta2shad[1] < interval[sn][1]){// beta2shad1 is located in this segment
// 		  betagood=1;
// 		  temp[0] = beta2shad[1];
// 		  temp[1] = interval[sn][1];
// 		  interval[sn][1] = beta2shad[0];
// 		  interval.insert(interval.begin()+sn+1,temp);
// 		  sn = interval.size();
// 		  if(verbose)cout << "Case 1: not checked" << endl;
// 		}
// 		else{//beta2shad1 is located in another segment
// 		  interval[sn][1] = beta2shad[0];
// 		  while(sn < interval.size()-1 && beta2shad[1]>interval[sn+1][1]){
// 		    interval.erase(interval.begin()+sn+1);
// 		  }
// 		  if(sn < interval.size()-1 && beta2shad[1]>interval[sn+1][0]){
// 		    interval[sn+1][0] = beta2shad[1];
// 		  }
// 		  if(verbose)cout << "Case 2: not checked" << endl;
// 		}
// 		sn = interval.size();
// 	      }
	      
// 	      else if((interval.size() > 1 && sn < interval.size()-1)&& (beta2shad[0]> interval[sn][1] && beta2shad[0] < interval[sn+1][0])){//beta2shad0 is located between two segments
// 		while(sn < interval.size()-1 && beta2shad[1] > interval[sn+1][1]){//remove segments that are in shadow
// 		  interval.erase(interval.begin()+sn+1);
// 		}
// 		if(sn < interval.size()-1 && beta2shad[1] > interval[sn+1][0]){//update segment that contains beta2shad1
// 		  // cout << "check here" << endl;
// 		  interval[sn+1][0] = beta2shad[1];
// 		}
// 		if(verbose)cout << "Case 3: not checked" << endl;
// 	      }
	      
// 	      else if(beta2shad[0] < interval[sn][0] && sn == 0){//beta2shad0 is located before first segment
// 		while(beta2shad[1] > interval[sn][1] && sn < interval.size()){//remove segments that are in shadow
// 		  interval.erase(interval.begin()+sn);
// 		}
// 		if(beta2shad[1] > interval[sn][0] && beta2shad[1] < interval[sn][1]){
// 		  interval[sn][0] = beta2shad[1];
// 		}
// 		if(verbose)cout << "Case 4: not checked" << endl;
// 	      }
// 	    }
// 	  }
// 	  else{//new segment wraps
// 	    if(beta2shad[0] < interval[interval.size()-1][1] && beta2shad[0] > interval[interval.size()-1][0] && interval[interval.size()-1][1] == 2*PI && interval[0][0]==0 && beta2shad[1] < interval[0][1]){//good wrap
// 	      betagood=1;
// 	      interval[interval.size()-1][1]=beta2shad[0];
// 	      interval[0][0]=beta2shad[1];
// 	      if(verbose)cout << "Case 5: not checked" << endl;
// 	    }
// 	    else{
// 	      for(int sn = 0;sn < interval.size(); sn++){
// 		while(beta2shad[0]< interval[sn][0] && sn < interval.size()){//remove segments that are in shadow
// 		  interval.erase(interval.begin()+sn);
// 		  if(verbose)cout << "Case 6: not checked" << endl;
// 		}
// 		while(beta2shad[1] > interval[sn][1] && sn < interval.size()){
// 		  interval.erase(interval.begin()+sn);
// 		  if(verbose)cout << "Case 7: not checked" << endl;
// 		}
// 		if( sn < interval.size() && beta2shad[1]> interval[sn][0] && beta2shad[1] < interval[sn][1]){// beta2shad1 is located in this segment
// 		  interval[sn][0] = beta2shad[1];
// 		  if(verbose)cout << "Case 8: not checked" << endl;
// 		}
// 		if( sn < interval.size() && beta2shad[0]>interval[sn][0] && beta2shad[0] < interval[sn][1]){//beta2shad0 is located in this segment
// 		  interval[sn][1] = beta2shad[0];
// 		  if(verbose)cout << "Case 9: not checked" << endl;
// 		}
// 	      }
// 	    }
// 	  }
// 	  //check intervals 
// 	  for(int sn = 0; sn < interval.size();sn++){
// 	    if(interval[sn][0]>interval[sn][1])cout << "Problem with intervals"  << endl;
// 	  }
// 	  if(betagood){
// 	    vertex_t * vpt = &Dists[dit].vertex;
// 	    double radius = Dists[dit].distance;
// 	    if(celltype == "abd"|| 
// 	       (celltype == "aa" && (*I)[i].allGraph[*vp1.first].type == "a" && (*I)[i].allGraph[*vpt].type == "a") || 
// 	       (celltype == "bb" && (*I)[i].allGraph[*vp1.first].type == "b" && (*I)[i].allGraph[*vpt].type == "b") || 
// 	       (celltype == "dd" && (*I)[i].allGraph[*vp1.first].type == "d" && (*I)[i].allGraph[*vpt].type == "d") || 
// 	       (celltype == "ab" && (*I)[i].allGraph[*vp1.first].type == "a" && (*I)[i].allGraph[*vpt].type == "b") ||
// 	       (celltype == "ad" && (*I)[i].allGraph[*vp1.first].type == "a" && (*I)[i].allGraph[*vpt].type == "d") ||
// 	       (celltype == "bd" && (*I)[i].allGraph[*vp1.first].type == "d" && (*I)[i].allGraph[*vpt].type == "b") ||
// 	       (celltype == "adb" && ((*I)[i].allGraph[*vp1.first].type == "a" ||(*I)[i].allGraph[*vp1.first].type == "d") &&(*I)[i].allGraph[*vpt].type == "b")  ||
// 	       (celltype == "adm" && 
// 		((*I)[i].allGraph[*vp1.first].type == "a" || (*I)[i].allGraph[*vp1.first].type == "d")  &&
// 		((*I)[i].allGraph[*vpt].type == "a" || (*I)[i].allGraph[*vpt].type == "d"))) {
	      
// 	      int arc = 0;
// 	      if(verbose)
// 		cout << "For vertex "<<(*I)[i].allGraph[*vp1.first].IsletCellNum << " (" << (*I)[i].allGraph[*vp1.first].x << 
// 		  "," << (*I)[i].allGraph[*vp1.first].y << ") with type " << (*I)[i].allGraph[*vp1.first].type << "  - " << 
// 		  (*I)[i].allGraph[*vpt].IsletCellNum << " (" << (*I)[i].allGraph[*vpt].x << "," << 
// 		  (*I)[i].allGraph[*vpt].y << ") with type " << (*I)[i].allGraph[*vpt].type << " with radius " << 
// 		  radius << " :" << endl;
	      
// 	      //First quadrant
// 	      for(int deg = 0;deg < 90;deg++){
// 		double theta = (double)deg*PI/180.0;
// 		double xs = x2+radius*cos(theta);
// 		double ys = y2+radius*sin(theta);
// 		if(verbose)cout << " for source deg " << deg << ": (x,y) is (" << xs << "," << ys << ")" << endl;
// 		if(xs <= xmax && xs >= xmin && ys <= ymax && ys >= ymin)
// 		  arc++;
// 		if(verbose)cout << "  Adding one to arc" << endl;
// 	      }
// 	      if(verbose)cout << ".....After 1st quad: arc is " << arc <<  endl;
// 	      //2nd quadrant
// 	      for(int deg = 0;deg < 90;deg++){
// 		double theta = (double)deg*PI/180.0;
// 		double xs = x2-radius*cos(theta);
// 		double ys = y2+radius*sin(theta);
// 		if(verbose)cout << " for source deg " << deg << ": (x,y) is (" << xs << "," << ys << ")" << endl;
// 		if(xs <= xmax && xs >= xmin && ys <= ymax && ys >= ymin)
// 		  arc++;
// 		if(verbose)cout << " Adding one to arc" << endl;
		
// 		if(verbose)cout << " for deg " << deg << ": (x,y) is (" << xs << "," << ys << ")" << endl;
// 	      }
// 	      if(verbose)cout << ".....After 2ND quad: arc is " << arc << endl;
// 	      //3rd quadrant
// 	      for(int deg = 0;deg < 90;deg++){
// 		double theta = (double)deg*PI/180.0;
// 		double xs = x2-radius*cos(theta);
// 		double ys = y2-radius*sin(theta);
// 		if(verbose)cout << " for source deg " << deg << ": (x,y) is (" << xs << "," << ys << ")" << endl;
// 		if(xs <= xmax && xs >= xmin && ys <= ymax && ys >= ymin)
// 		  arc++;
// 		if(verbose)cout << " Adding one to arc" << endl;
// 	      }
// 	      if(verbose)cout << ".....After 3rd quad: arc is " << arc <<  endl;
// 	      //4th quadrant
// 	      for(int deg = 0;deg < 90;deg++){
// 		double theta = (double)deg*PI/180.0;
// 		double xs = x2+radius*cos(theta);
// 		double ys = y2-radius*sin(theta);
// 		if(verbose)cout << " for source deg " << deg << ": (x,y) is (" << xs << "," << ys << ")" << endl;
// 		if(xs <= xmax && xs >= xmin && ys <= ymax && ys >= ymin)
// 		arc++;
// 		if(verbose)cout << " Adding one to arc" << endl;
// 	      }
// 	      if(verbose)cout << ".....After 4th quad: arc is " << arc <<   endl;
// 	      double arcdegrees = (double)arc*PI/180.0;
// 	      if(verbose)cout << "...Arcdegrees = " << arcdegrees <<  endl;
	      
// 	      if(arc !=0){
// 		if((int)(radius/dr)>=gk_norm[i].size()-1){
// 		  if(verbose)cout << "radius/dr= " << radius/dr << "gk_norm.size is " << gk_norm[i].size()-1 << endl;
// 		}
// 		else gk_norm[i][(int)(radius/dr)]+=IsletArea/(arcdegrees*radius*dr*numvertsq);
// 	      }
// 	    }
// 	  }//closes vp2
// 	  dit++;
// 	}//closes vp1
// 	if(verbose){
// 	  //Printout gknorm values
// 	  cout << "For islet " << i << ": gknorm values are " ;
// 	  for(int j = 0; j < gk_norm[i].size();j++){
// 	    cout << gk_norm[i][j] << " ";
// 	  }
// 	  cout << endl;
// 	}
//       }
//     }
//   }//closes Graphs
  
//   for(int j = 0; j < nbins;j++){
//     vector<double> avggk(2,0);vectorcount++;if(verbosemem)cout << "Adding avggk: vectorcount = " << vectorcount << endl; cout.flush();
//     avggk[0]=dr*j;
//     for(int i = 0; i < (*I).size();i++){
//       avggk[1]+= gk_norm[i][j];
//     }
//     if(numIslets!=0)avggk[1] = avggk[1]/(double)numIslets;
//     else avggk[1]=0;
//     (*Gofr).push_back(avggk);
//     RemoveVector(&avggk);vectorcount--;if(verbosemem)cout << "Removing avggk: vectorcount = " << vectorcount << endl; cout.flush();
//   }
//   RemoveVector2d(&gk_norm);vectorcount--;if(verbosemem)cout << "Removing gknorm: vectorcount = " << vectorcount << endl; cout.flush();
//   if(vectorcount > 0)cout << "Problem with vectors! vectorcount = " << vectorcount << " in CalcGofrUsingShadow()" << endl;
//   cout << " Finished" << endl;
//   return ;
// }

// //*********************************************************************************************//
// //*********************************************************************************************//
// //*********************************************************************************************//



// void CalcGofrUsingHull(vector<IsletGraph>* Gp,vector<vector<double> > *Gofr, string celltype,int Gno){ 
//   std::pair<vertex_iter, vertex_iter>vp,vp1,vp2;
  
//   if(Gno != 99999)cout << "Calculating Gofr for Islet " << Gno << " out of " << (*Gp).size() << "...";
//   else cout << "Calculating Gofr ..." << endl;
//   int vectorcount = 0;
//   bool verbose=0,verbosemem=0;
//   int nbins = 1000;
//   double dr = 0.5;
//   vector<vector<double> > gk_norm((*Gp).size(),vector<double>(nbins,0.0));vectorcount++;if(verbosemem)cout << "Adding gk_norm: vectorcount = " << vectorcount << endl; cout.flush();
//   int numIslets = 0;
//   int Gbeg = 0, Gend = (*Gp).size();
//   if(Gno!=99999){
//     Gbeg = Gno;
//     Gend = Gno+1;
//   }

//   for(int i = Gbeg; i < Gend;i++){
//     if(verbose &&  i%100 == 0){
//       cout << i << " " ;
//       cout.flush();
//     }
    
//     //Determine the (#vertices)^2
//     int numvertsq = 0;
//     if (celltype == "abd")  numvertsq = (int)pow(num_vertices((*Gp)[i].allGraph),2.0);
//     else if(celltype == "aa")numvertsq = (int)pow((*Gp)[i].nalpha,2.0);
//     else if(celltype == "bb")numvertsq = (int)pow((*Gp)[i].nbeta,2.0);
//     else if(celltype == "dd")numvertsq = (int)pow((*Gp)[i].ndelta,2.0); 
//     else if(celltype == "adm")numvertsq = (int)pow((*Gp)[i].nalpha+(*Gp)[i].ndelta,2.0);
//     else if(celltype == "ab")numvertsq = (*Gp)[i].nalpha*(*Gp)[i].nbeta;
//     else if(celltype == "ad")numvertsq = (*Gp)[i].nalpha*(*Gp)[i].ndelta;
//     else if(celltype == "bd")numvertsq = (*Gp)[i].ndelta*(*Gp)[i].nbeta;
//     else if(celltype == "adb")numvertsq = ((*Gp)[i].nalpha+(*Gp)[i].ndelta)*(*Gp)[i].nbeta;
    
//     if((celltype == "abd" && num_vertices((*Gp)[i].allGraph) > 1)||
//        (celltype == "aa" && (*Gp)[i].nalpha > 1)||
//        (celltype == "bb" && (*Gp)[i].nbeta > 1)||
//        (celltype == "dd" && (*Gp)[i].ndelta > 1)||
//        (celltype == "adm" && (*Gp)[i].nalpha+(*Gp)[i].ndelta> 1)||
//        (celltype == "ab" && (*Gp)[i].nalpha > 0 && (*Gp)[i].nbeta > 0) ||
//        (celltype == "ad" && (*Gp)[i].nalpha > 0 && (*Gp)[i].ndelta > 0) ||
//        (celltype == "bd" && (*Gp)[i].nbeta > 0 && (*Gp)[i].ndelta > 0) ||
//        (celltype =="adb" && ((*Gp)[i].nalpha > 0 ||(*Gp)[i].ndelta >0) && (*Gp)[i].nbeta > 0)){
//       numIslets++;
      
//       //Determine area of hull
//       double IsletArea = PolygonArea(&((*Gp)[i].allGraph),&((*Gp)[i].hullPath));
      
//       //if(verbose)cout << "For Islet " << (*Gp)[i].IsletNum <<"\t" <<  xmin << "\t" << xmax << "\t" << xlength << "\t" <<  ymin << "\t" << ymax << "\t" << ylength << endl;
      
//       //calculate the number of possible edges
//       for(vp1=vertices((*Gp)[i].allGraph);vp1.first!=vp1.second;vp1.first++){
// 	for(vp2=vertices((*Gp)[i].allGraph);vp2.first!=vp2.second;vp2.first++){
// 	  if((*vp1.first != *vp2.first)&&
// 	     (celltype == "abd"|| 
// 	      (celltype == "aa" && (*Gp)[i].allGraph[*vp1.first].type == "a" && (*Gp)[i].allGraph[*vp2.first].type == "a") || 
// 	      (celltype == "bb" && (*Gp)[i].allGraph[*vp1.first].type == "b" && (*Gp)[i].allGraph[*vp2.first].type == "b") || 
// 	      (celltype == "dd" && (*Gp)[i].allGraph[*vp1.first].type == "d" && (*Gp)[i].allGraph[*vp2.first].type == "d") || 
// 	      (celltype == "ab" && (*Gp)[i].allGraph[*vp1.first].type == "a" && (*Gp)[i].allGraph[*vp2.first].type == "b") ||
// 	      (celltype == "ad" && (*Gp)[i].allGraph[*vp1.first].type == "a" && (*Gp)[i].allGraph[*vp2.first].type == "d") ||
// 	      (celltype == "bd" && (*Gp)[i].allGraph[*vp1.first].type == "d" && (*Gp)[i].allGraph[*vp2.first].type == "b") ||
// 	      (celltype == "adb" && ((*Gp)[i].allGraph[*vp1.first].type == "a" ||(*Gp)[i].allGraph[*vp1.first].type == "d") &&(*Gp)[i].allGraph[*vp2.first].type == "b")  ||
// 	      (celltype == "adm" && 
// 	       ((*Gp)[i].allGraph[*vp1.first].type == "a" || (*Gp)[i].allGraph[*vp1.first].type == "d")  &&
// 	       ((*Gp)[i].allGraph[*vp2.first].type == "a" || (*Gp)[i].allGraph[*vp2.first].type == "d")))) {
	    
// 	    double x1 = ((*Gp)[i].allGraph)[*vp1.first].x, y1 = ((*Gp)[i].allGraph)[*vp1.first].y;
// 	    double x2 = ((*Gp)[i].allGraph)[*vp2.first].x, y2 = ((*Gp)[i].allGraph)[*vp2.first].y;
	      
// 	    double radius = Distancexy(x1,y1,x2,y2);
// 	    int arc = 0;
// 	    if(verbose)cout << "For vertex "<<(*Gp)[i].allGraph[*vp1.first].IsletCellNum << " (" << (*Gp)[i].allGraph[*vp1.first].x << 
// 	      "," << (*Gp)[i].allGraph[*vp1.first].y << ") with type " << (*Gp)[i].allGraph[*vp1.first].type << "  - " << 
// 	      (*Gp)[i].allGraph[*vp2.first].IsletCellNum << " (" << (*Gp)[i].allGraph[*vp2.first].x << "," << (*Gp)[i].allGraph[*vp2.first].y << 
// 	      ") with type " << (*Gp)[i].allGraph[*vp2.first].type << " with radius " << radius << " :" << endl;
	    
// 	    //First quadrant
// 	    for(int deg = 0;deg < 90;deg++){
// 	      double theta = (double)deg*PI/180.0;
// 	      double xs = x2+radius*cos(theta);
// 	      double ys = y2+radius*sin(theta);
// 	      if(verbose)cout << " for source deg " << deg << ": (x,y) is (" << xs << "," << ys << ")" << endl;
// 	      int wn = WindingNumber(&((*Gp)[i].allGraph), &((*Gp)[i].hullPath),xs, ys); 
// 	      if(wn !=0)arc++;
// 	      if(verbose)cout << "  Adding one to arc" << endl;
// 	    }
// 	    if(verbose)cout << ".....After 1st quad: arc is " << arc <<  endl;
// 	    //2nd quadrant
// 	    for(int deg = 0;deg < 90;deg++){
// 	      double theta = (double)deg*PI/180.0;
// 	      double xs = x2-radius*cos(theta);
// 	      double ys = y2+radius*sin(theta);
// 	      if(verbose)cout << " for source deg " << deg << ": (x,y) is (" << xs << "," << ys << ")" << endl;
// 	      int wn = WindingNumber(&((*Gp)[i].allGraph), &((*Gp)[i].hullPath),xs, ys);
// 	      if(wn !=0)arc++;
// 	      if(verbose)cout << " Adding one to arc" << endl;
	      
// 	      if(verbose)cout << " for deg " << deg << ": (x,y) is (" << xs << "," << ys << ")" << endl;
// 	    }
// 	    if(verbose)cout << ".....After 2ND quad: arc is " << arc << endl;
// 	    //3rd quadrant
// 	    for(int deg = 0;deg < 90;deg++){
// 	      double theta = (double)deg*PI/180.0;
// 	      double xs = x2-radius*cos(theta);
// 	      double ys = y2-radius*sin(theta);
// 	      if(verbose)cout << " for source deg " << deg << ": (x,y) is (" << xs << "," << ys << ")" << endl;
// 	      int wn = WindingNumber(&((*Gp)[i].allGraph), &((*Gp)[i].hullPath),xs, ys);
// 	      if(wn !=0)arc++;
// 	      if(verbose)cout << " Adding one to arc" << endl;
// 	    }
// 	    if(verbose)cout << ".....After 3rd quad: arc is " << arc <<  endl;
// 	    //4th quadrant
// 	    for(int deg = 0;deg < 90;deg++){
// 	      double theta = (double)deg*PI/180.0;
// 	      double xs = x2+radius*cos(theta);
// 	      double ys = y2-radius*sin(theta);
// 	      if(verbose)cout << " for source deg " << deg << ": (x,y) is (" << xs << "," << ys << ")" << endl;
// 	      int wn = WindingNumber(&((*Gp)[i].allGraph), &((*Gp)[i].hullPath),xs, ys);
// 	      if(wn !=0)arc++;
// 	      if(verbose)cout << " Adding one to arc" << endl;
// 	    }
// 	    if(verbose)cout << ".....After 4th quad: arc is " << arc <<   endl;
// 	    double arcdegrees = (double)arc*PI/180.0;
// 	    if(verbose)cout << "...Arcdegrees = " << arcdegrees <<  endl;
	    
// 	    if(arc !=0){
// 	      if((int)(radius/dr)>=gk_norm[i].size()-1){
// 		if(verbose)cout << "radius/dr= " << radius/dr << "gk_norm.size is " << gk_norm[i].size()-1 << endl;
// 	      }
// 	      else gk_norm[i][(int)(radius/dr)]+=IsletArea/(arcdegrees*radius*dr*numvertsq);
// 	    }
// 	  }
// 	}//closes vp2 
//       }//closes vp1
//       if(verbose){
// 	//Printout gknorm values
// 	cout << "For islet " << i << ": gknorm values are " ;
// 	for(int j = 0; j < gk_norm[i].size();j++){
// 	  cout << gk_norm[i][j] << " ";
// 	}
// 	cout << endl;
//       }
//     }
//   }//closes Graphs
  
//   for(int j = 0; j < nbins;j++){
//     vector<double> avggk(2,0);vectorcount++;if(verbosemem)cout << "Adding avggk: vectorcount = " << vectorcount << endl; cout.flush();
//     avggk[0]=dr*j;
//     for(int i = 0; i < (*Gp).size();i++){
//       avggk[1]+= gk_norm[i][j];
//     }
//     if(numIslets!=0)avggk[1] = avggk[1]/(double)numIslets;
//     else avggk[1]=0;
//     (*Gofr).push_back(avggk);
//     RemoveVector(&avggk);vectorcount--;if(verbosemem)cout << "Removing avggk: vectorcount = " << vectorcount << endl; cout.flush();
//   }
//   RemoveVector2d(&gk_norm);vectorcount--;if(verbosemem)cout << "Removing gknorm: vectorcount = " << vectorcount << endl; cout.flush();
//   if(vectorcount > 0)cout << "Problem with vectors! vectorcount = " << vectorcount << " in CalcGofr()" << endl;
//   cout << " Finished" << endl;
//   return ;
// }

// //*********************************************************************************************//
// //*********************************************************************************************//
// //*********************************************************************************************//

// double PolygonArea(Graph *G,vector<vertex_t> * Path){
//  http://www.cdn.geeksforgeeks.org/area-of-a-polygon-with-given-n-ordered-vertices/
//   double area = 0; 
//   // Calculate value of shoelace formula
//   int j = (*Path).size()- 2;
//   for (int i = 0; i < (*Path).size()-1; i++){
//     area += ((*G)[(*Path)[j]].x + (*G)[(*Path)[i]].x) * ((*G)[(*Path)[j]].y - (*G)[(*Path)[i]].y);
//     j = i;  // j is previous vertex to i
//   }
  
//   // Return absolute value
//   return fabs(area / 2.0);
// }

// //*********************************************************************************************//
// //*********************************************************************************************//
// //*********************************************************************************************//

// bool CleanIntersectionPointsInPaths(Graph * Gp, vector<vertex_t> * path, vector<vertex_t> *LeaveIP, bool isSetLoop, vector<vector<double> >* set){
 
//   std::pair<edge_iter, edge_iter> ep1,ep2;
//   vertex_t u1, v1, u2, v2, u;
//   edge_t e;
//    std::pair<vertex_iter, vertex_iter> vp ;
//   bool verbose=1,verbosemem=0;
//   bool ContainsIPs = 1;
//   int graphcount = 0,vectorcount = 0,pathcount = 0;
  
//   if(verbose){
//     PrintoutGnuplotWithIPs(Gp, "temp");
//     cout << "**************Beginning of CleanIntersectionPointsInPaths **********************" << endl;
//   }
//   //Make copy of Graph
//   Graph Gpcopy; graphcount++;if(verbosemem)cout << "Adding Gpcopy (graph): graphcount = " << graphcount << endl; cout.flush();
//   //if(verbose)PrintoutGnuplotWithIPs(Gp,"temp.beforecopyMinusIP");
//   CopyGraphMinusIP(Gp,&Gpcopy);
//   int iter = 0;
//   while(ContainsIPs){
//     bool fixed = 0;
//     if(verbose)cout << "Inside CleanIntermediatePoints" << endl;
    
//     if(verbose){
//       cout << "The path contains " ;
//       cout.flush();
//       for(int i = 0; i < (*path).size(); i++)cout << (*Gp)[(*path)[i]].IsletCellNum << " ";cout.flush();
//       cout << endl;
//       cout.flush();
//     }
    
//     //Check if all vertices in paths are IPs
//     int numNonIPs = 0;
//     for(int i = 0; i < (*path).size()-1;i++){
//       if((*Gp)[(*path)[i]].CellNum <100000){
// 	numNonIPs++;
//       }
//     }
//     if(verbose)cout << "The number of nonIPs is " << numNonIPs << " out of " << (*path).size() << " cells in path " << endl;cout.flush();
//     int endpt1pos = 0,endpt2pos=0;
//     vertex_t endpt1,endpt2;
//     int ippos = 0;
//     vertex_t IP;
//     //bool deg1endpt1 = 0, deg1endpt2=0, deg2endpt = 0,pendant = 0;
//     if(numNonIPs<1){
//       endpt1pos = -1;
//       endpt2pos = (*path).size();
//       IP = (*path)[0];
//       ippos = 0;
//       if(verbose)cout << "Endpt1 is at pos " << endpt1pos << " and endpt2 is at pos "<< endpt2pos << endl;cout.flush();
//     }
      
//     bool allconnected= 1;
   
//     int ipn1 = 1;
//     //Make copy of Graph
//     Graph Gpcopy2; graphcount++;if(verbosemem)cout << "Adding Gpcopy2 (graph): graphcount = " << graphcount << endl; cout.flush();
//     CopyGraphMinusIP(Gp,&Gpcopy2);
//     //if(verbose)PrintoutGnuplot(&Gpcopy2,"temp.aftercopyMinusIP");
//     if(numNonIPs >=1){ //has nonIPs
//       if(verbose) cout << "path has nonIPs" << endl;
//       //if first entry is an IP: move to back until first entry is not an IP
//       while ((*Gp)[(*path)[0]].CellNum > 99999){
// 	if((*path)[0] == (*path)[(*path).size()-1])
// 	  (*path).erase((*path).begin());
// 	else{
// 	  (*path).push_back((*path)[0]);
// 	  (*path).erase((*path).begin());
// 	}
//       }
      
//       //make sure last entry is same as first in path
//       if((*path)[0] != (*path)[(*path).size()-1])
// 	(*path).push_back((*path)[0]);
//       if(verbose){
// 	cout << "After rotating,The path contains " ;
// 	for(int i = 0; i < (*path).size(); i++)cout << (*Gp)[(*path)[i]].IsletCellNum << " ";
// 	cout << endl;
//       }
      
//       //Find first IP
//       bool foundIP = 0;
//       while(!foundIP){
//       	while((*Gp)[(*path)[ipn1]].CellNum < 99999)ipn1++; 
//       	//Check if IP is to be left
//       	bool inLeaveIP = 0;
//       	for(int i = 0; i < (*LeaveIP).size();i++){
//       	  if((*path)[ipn1] ==  (*LeaveIP)[i]){
//       	    inLeaveIP = 1;
//       	    if(verbose)cout << "Cell " << (*Gp)[(*path)[ipn1]].IsletCellNum << " is a permanent IP.  Looking for next..." << endl;cout.flush();
//       	    break;
//       	  }
//       	}
// 	if(!inLeaveIP ){// !inIPs){
// 	  IP = (*path)[ipn1];
// 	  int ib = ipn1-1;
// 	  bool foundep1 = 0;
// 	  while(!foundep1){
// 	    if((*Gp)[(*path)[ib]].CellNum > 99999){
// 	      //Check if in LeaveIPs
// 	      for(int i = 0; i < (*LeaveIP).size();i++){
// 		if((*path)[ib]== (*LeaveIP)[i]){
// 		  foundep1 = 1;
// 		  i = (*LeaveIP).size();
// 		}
// 	      }
// 	    }
// 	    else foundep1 = 1;
// 	    if(!foundep1)
// 	      ib--;
//  	  }
// 	  endpt1pos = ib;
// 	  endpt1 = (*path)[ib];
// 	  int ii = ipn1;
// 	  bool foundep2 = 0;
// 	  while(!foundep2){
// 	    if((*Gp)[(*path)[ii]].CellNum > 99999){
// 	      //Check if in LeaveIPs
// 	      for(int i = 0; i < (*LeaveIP).size();i++){
// 		if((*path)[ii]== (*LeaveIP)[i]){
// 		  foundep2 = 1;
// 		  i = (*LeaveIP).size();
// 		}
// 	      }
// 	    }
// 	    else foundep2 = 1;
// 	    if(!foundep2)
// 	      ii++;
// 	  }
// 	  endpt2pos = ii;
// 	  endpt2 = (*path)[ii];
// 	  foundIP = 1;
// 	}
// 	else ipn1++;
//       }
    
//       if(verbose)cout << "The first IP is " << (*Gp)[IP].IsletCellNum << " at position " << ipn1 << " with parents (" <<(*Gp)[(*Gp)[IP].IPparents[0]].IsletCellNum << "," <<(*Gp)[(*Gp)[IP].IPparents[1]].IsletCellNum << ","<<(*Gp)[(*Gp)[IP].IPparents[2]].IsletCellNum << ","<<(*Gp)[(*Gp)[IP].IPparents[3]].IsletCellNum << ")" <<  endl;

//       //if((*Gp)[IP].IsletCellNum == 10193) PrintoutGnuplot(Gp,"temp");
      
//       //remove edges between IP parents
//       bool b = 0;
//       tie(e,b)=edge((*Gp)[IP].IPparents[0],(*Gp)[IP].IPparents[1],Gpcopy2);
//       ClearEdge(&Gpcopy2,e);
//       remove_edge(e,Gpcopy2);
//       tie(e,b)=edge((*Gp)[IP].IPparents[2],(*Gp)[IP].IPparents[3],Gpcopy2);
//       ClearEdge(&Gpcopy2,e);
//       remove_edge(e,Gpcopy2);
//       //PrintoutGnuplot(Gp,"tempbeforeRemoveIPedges");
//       //PrintoutGnuplot(&Gpcopy, "tempRemoveIPedges");
      
//       //See if path exists between IPparents
//       for(int i = 0; i < (*Gp)[IP].IPparents.size(); i++){
// 	for(int j = i+1; j <(*Gp)[IP].IPparents.size(); j++){
// 	  vector<vertex_t> temppath;vectorcount++;if(verbosemem)cout << "Adding temppath: vectorcount = " << vectorcount << endl; cout.flush();

// 	  bool nopath = PathBetweenTwoVertices(&Gpcopy2, (*Gp)[IP].IPparents[i], (*Gp)[IP].IPparents[j],&temppath);
// 	  RemoveVector(&temppath);vectorcount--;if(verbosemem)cout << "Removing temppath: vectorcount = " << vectorcount << endl; cout.flush();
// 	  if(nopath){
// 	    if(verbose)cout << "Path does not exist between " << (*Gp)[(*Gp)[IP].IPparents[i]].IsletCellNum << " and " << (*Gp)[(*Gp)[IP].IPparents[j]].IsletCellNum << endl;
// 	    allconnected = 0;
// 	    i =(*Gp)[IP].IPparents.size() ;
// 	    break;
// 	  }
// 	}
//       }
      
//       if(!allconnected){
// 	(*LeaveIP).push_back(IP);
// 	//add IP to Gpcopy
	
//       }
//     }
    
//     if(allconnected){
//       if(verbose)cout << "Inside allconnected" << endl;cout.flush();
//       if(verbose &&numNonIPs > 1)cout << "The endpts found are " << (*Gp)[endpt1].IsletCellNum << " at pos " << endpt1pos << " and " << (*Gp)[endpt2].IsletCellNum << " at pos "<< endpt2pos << endl;cout.flush();

//       if(numNonIPs > 0 && (*Gp)[endpt1].CellNum > 9999){
//       	//Check if IP's parents match endpt1's parents
//       	if(verbose){
//       	  cout << "Endpt1's IP parents are: ";
//       	  for(int i = 0;i < 4;i++)
//       	    cout << (*Gp)[(*Gp)[endpt1].IPparents[i]].IsletCellNum << " ";
//       	  cout << endl;
//       	}
//       	//check if endpt1 is inline with IP and any IPparents
//       	for(int i = 0; i < (*Gp)[IP].IPparents.size(); i++){
//       	  vertex_t ipp = (*Gp)[IP].IPparents[i];
//       	  double dist = DistanceFromPointToLineSegment((*Gp)[endpt1].x, (*Gp)[endpt1].y,
//       						       (*Gp)[IP].x,(*Gp)[IP].y,
//       						       (*Gp)[ipp].x,(*Gp)[ipp].y);
// 	  if(verbose)cout << "Distance from endpt1 to line segment " << (*Gp)[IP].IsletCellNum << " and " << (*Gp)[ipp].IsletCellNum << " is " << dist << endl; 
//       	  if(fabs(dist) < 1e-4){//Replace IP parent with endpt1
//       	    if(verbose)cout << "updating ip parent" << endl;
//       	    (*Gp)[IP].IPparents[i] =(*Gp)[endpt1].CellNum;
// 	  }
//       	}
//       }

//       //Sort IPparents to be added in to path
//       //Create initial list of IPparents
//       vector<VertexDistance> IPparents;vectorcount++;if(verbosemem)cout << "Adding IPparents: vectorcount = " << vectorcount << endl; cout.flush();
//       int pos = endpt1pos+1;
//       while(pos !=endpt2pos){
//       	if(verbose)cout << "Adding parents from IP at pos " << pos << endl;cout.flush();
//       	for(int p = 0; p <(*Gp)[(*path)[pos]].IPparents.size() ; p++){
//       	  if(verbose) cout << "Checking cell " << (*Gp)[(*Gp)[(*path)[pos]].IPparents[p]].IsletCellNum << endl;cout.flush();
//       	  //Check if parent not in IPparents or an endpt
//       	  bool inParent = 0;
//       	  for(int i = 0; i < IPparents.size();i++){
//       	    if((*Gp)[IPparents[i].vertex].CellNum == (*Gp)[(*Gp)[(*path)[pos]].IPparents[p]].CellNum){
//       	      inParent = 1;
//       	      break;
//       	    }
//       	  }
//       	  if(numNonIPs > 0){
//       	    if((*Gp)[(*Gp)[(*path)[pos]].IPparents[p]].CellNum==(*Gp)[endpt1].CellNum)inParent = 1;
//       	    if((*Gp)[(*Gp)[(*path)[pos]].IPparents[p]].CellNum==(*Gp)[endpt2].CellNum)inParent = 1;
//       	  }
//       	  if(verbose && numNonIPs > 0)cout << "inParent = " << inParent <<" " <<  (*Gp)[(*Gp)[(*path)[pos]].IPparents[p]].IsletCellNum << " " << (*Gp)[endpt1].IsletCellNum << endl;
	  
//       	  if(!inParent){
//       	    if(verbose) cout << "  Adding cell " << (*Gp)[(*Gp)[(*path)[pos]].IPparents[p]].IsletCellNum << endl;
//       	    VertexDistance temp;
//       	    temp.vertex = (*Gp)[(*path)[pos]].IPparents[p];
//       	    if(numNonIPs ==0)temp.distance = 0;
//       	    else temp.distance = Distancexy((*Gp)[endpt1].x, (*Gp)[endpt1].y, (*Gp)[temp.vertex].x, (*Gp)[temp.vertex].y);
//       	    IPparents.push_back(temp);
//       	  }
//       	}
//       	pos++;
//       }
//       if(numNonIPs == 0){
//       	//Find IPparent closest to set and
//       	double mindist = 1000000;
//       	for(int s = 0; s < (*set).size();s++){
//       	  for(int i = 0; i < IPparents.size();i++){
//       	    double dist = Distancexy((*set)[s][0],(*set)[s][1],(*Gp)[IPparents[i].vertex].x, (*Gp)[IPparents[i].vertex].y);
//       	    if(dist < mindist){
//       	      mindist = dist;
//       	      endpt1 = IPparents[i].vertex;
//       	    }
//       	  }
//       	}
//       	//Set IPparent.distance wrt closest
//       	for(int i = 0; i < IPparents.size();i++)
//       	  IPparents[i].distance = Distancexy((*Gp)[endpt1].x, (*Gp)[endpt1].y,(*Gp)[IPparents[i].vertex].x, (*Gp)[IPparents[i].vertex].y);
//       	if(verbose)cout << "Endpt1 is "<<  (*Gp)[endpt1].IsletCellNum  << endl;
//       }
      
//       //sort for closest parent to endpt1
//       sort(IPparents.begin(), IPparents.end(),compareByDistance);
//       if(numNonIPs == 0)IPparents.erase(IPparents.begin());
//       if(verbose) {
//       	cout << "Finished sorting. IPParents contains: " << endl;
//       	for(int i =  0; i < IPparents.size();i++)
//       	  cout << (*Gp)[IPparents[i].vertex].IsletCellNum << "\t" << IPparents[i].distance << endl; 
//       }
      

//       //****************New code for creating temppath************************
//       //if only one nonIP: set closest IPparent to endpt2 and erase
//       if(numNonIPs < 2){
//       	endpt2 = IPparents[0].vertex;
//       	IPparents.erase(IPparents.begin());
//       }
      
//       //Find center of IP parents
//       double xc = 0, yc = 0;
//       for(int i = 0; i < IPparents.size();i++){
//       	xc += (*Gp)[IPparents[i].vertex].x;
//       	yc += (*Gp)[IPparents[i].vertex].y;
//       }
//       xc+= (*Gp)[endpt1].x; yc+= (*Gp)[endpt1].y;
//       xc+= (*Gp)[endpt2].x; yc+= (*Gp)[endpt2].y;
//       xc /= (double)(IPparents.size()+2);
//       yc /= (double)(IPparents.size()+2);
//       if(verbose)cout << "The center of IP parents is at (" << xc << ","<< yc << ")" << endl; cout.flush(); 

//       vector<CellAngle> CAngles;vectorcount++;if(verbosemem)cout << "Adding vector CAngles: vectorcount = " << vectorcount << endl;cout.flush();
//       CellAngle tangle;
//       tangle.Cell = endpt1;
//       tangle.Angle = GetAngleOfLine(xc,yc, (*Gp)[endpt1].x,(*Gp)[endpt1].y);
//       CAngles.push_back(tangle);
//       tangle.Cell = endpt2;
//       tangle.Angle = GetAngleOfLine(xc,yc, (*Gp)[endpt2].x,(*Gp)[endpt2].y);
//       CAngles.push_back(tangle);
//       for(int i = 0; i < IPparents.size();i++){
//       	tangle.Cell = IPparents[i].vertex;
//       	tangle.Angle = GetAngleOfLine(xc,yc, (*Gp)[IPparents[i].vertex].x,(*Gp)[IPparents[i].vertex].y);
//       	CAngles.push_back(tangle);
//       }
//       //sort CAngles
//       sort(CAngles.begin(),CAngles.end(),compareByAngle);
//       if(verbose){
//       	cout << "After sorting CAngles is :" << endl;
//       	    for(int ca =0; ca < CAngles.size();ca++)cout <<  (*Gp)[CAngles[ca].Cell].IsletCellNum << "\t" << CAngles[ca].Angle << endl;
//       	    cout << endl;
//       	  }
//       //Rotate until first entry is endpt1
//       while(CAngles[0].Cell != endpt1){
//       	CAngles.push_back(CAngles[0]);
//       	CAngles.erase(CAngles.begin());
//       }
//       // //if the 2nd entry is endpt2: change directions
//       // if(CAngles[1].Cell==endpt2){
//       // 	sort(CAngles.begin(),CAngles.end(),compareByAngleBiggestFirst);
//       // 	while(CAngles[0].Cell != endpt1){
//       // 	  CAngles.push_back(CAngles[0]);
//       // 	  CAngles.erase(CAngles.begin());
//       // 	}
//       // }
//       //if the last entry isn't endpt2: flip remainder
//       for(int i = 1;i<CAngles.size()-1;i++){
//       	if(CAngles[i].Cell == endpt2){
//       	  if(verbose){
//       	    cout << "Before rotating CAngles is :" ;
//       	    for(int ca =0; ca < CAngles.size();ca++)cout <<  (*Gp)[CAngles[ca].Cell].IsletCellNum << " ";
//       	    cout << endl;
//       	  }
//       	  //copy remainder of list in reverse order
//       	  vector<CellAngle> tempca; vectorcount++;
//       	  for(int j = CAngles.size()-1;j >= i;j--)
//       	    tempca.push_back(CAngles[j]);
//       	  if(verbose){
//       	    cout << "Tempca is :" ;
//       	    for(int ca =0; ca < tempca.size();ca++)cout <<  (*Gp)[tempca[ca].Cell].IsletCellNum << " ";
//       	    cout << endl;
//       	  }
//       	  //remove remainder
//       	  CAngles.erase(CAngles.begin()+i,CAngles.end());
//       	  if(verbose){
//       	    cout << "After removing: CAngles is :" ;
//       	    for(int ca =0; ca < CAngles.size();ca++)cout <<  (*Gp)[CAngles[ca].Cell].IsletCellNum << " ";
//       	    cout << endl;
//       	  }
//       	  //Add in reversed segment
//       	  CAngles.insert(CAngles.begin()+i,tempca.begin(),tempca.end());
//       	  if(verbose){
//       	    cout << "After adding back: CAngles is :" ;
//       	    for(int ca =0; ca < CAngles.size();ca++)cout <<  (*Gp)[CAngles[ca].Cell].IsletCellNum << " ";
//       	    cout << endl;
//       	  }	
//       	  RemoveVector(&tempca);vectorcount--;
//       	}
//       }
//       //Setup temppath
//       vector<Path> Temppath;vectorcount++;if(verbosemem)cout << "Adding Temppath: vectorcount = " << vectorcount << endl; cout.flush();
//       Path temp;pathcount++;if(verbosemem)cout << "Adding temp(Path): pathcount = " << pathcount << endl; cout.flush();
//       for(int i = 0; i < CAngles.size()-1;i++){
//       	temp.cell1 = CAngles[i].Cell;
//       	temp.cell2 = CAngles[i+1].Cell;
//       	Temppath.push_back(temp);
//       }
//       RemoveVector(&CAngles);vectorcount--;if(verbosemem)cout << "Removing vector CAngles: vectorcount = " << vectorcount << endl;cout.flush();


//       //Add path to endpt2
//       //if(numNonIPs > 1){
//       	temp.cell1 = endpt2;
//       	if(endpt2pos < (*path).size()-1)temp.cell2 = (*path)[endpt2pos+1];
//       	else temp.cell2 = (*path)[0];
//       	Temppath.push_back(temp);
//       	//	}
      
      
//       RemoveVector(&IPparents);vectorcount--;if(verbosemem)cout << "Removing IPparents: vectorcount = " << vectorcount << endl; cout.flush();
//       RemovePath(&temp);pathcount--;if(verbosemem)cout << "Removing temp: pathcountcount = " << pathcount << endl; cout.flush();
      
      
//       //New code replacing dijkstra
//       //Find shortest paths between all vertices in temppath and create shortest path vector
//       vector<Path>ShortestPaths; vectorcount++;if(verbosemem)cout << "Adding ShortestPaths: vectorcount = " << vectorcount << endl; cout.flush();
//       //Make sure LeaveIPs are in Gpcopy
//       int lipnum = 0;
//       for(int i = 0; i < (*LeaveIP).size();i++){
//       	if(verbose)cout << "Checking to see if cell " << (*Gp)[(*LeaveIP)[i]].IsletCellNum << " is in Gpcopy" << endl;cout.flush();
//       	bool inGpcopy = 0;
//       	for(vp = vertices(Gpcopy);vp.first!=vp.second;vp.first++){
//       	  if(verbose)cout << "  vp is at cell " << Gpcopy[*vp.first].IsletCellNum << endl; cout.flush(); 
//       	  if((*LeaveIP)[i] == *vp.first){
//       	    cout << "These two are equal" << endl;
//       	    inGpcopy = 1;
//       	    lipnum = i;
//       	    break;
//       	  }
//       	}
//       	if(!inGpcopy){	 
//       	  u = add_vertex(Gpcopy);
//       	  InitializeVertexFromVertex(Gp,&Gpcopy,(*LeaveIP)[i],u,(*Gp)[(*LeaveIP)[i]].CellNum);
//       	  for(int j = 0; j < 4;j++){
//       	    vertex_t ipp = (Gpcopy)[u].IPparents[j];
//       	    bool b;
//       	    tie(e,b)=add_edge(u,ipp,Gpcopy);
//       	    double dist = Distancexy(Gpcopy[u].x, Gpcopy[u].y,Gpcopy[ipp].x,Gpcopy[ipp].y);
//       	    InitializeEdge(&Gpcopy,e,dist,dist,2);
//       	  }
//       	  if(verbose)
//       	    cout << "Added cell " << Gpcopy[u].IsletCellNum << " to Gpcopy" << endl;
	  
//       	}
//       }
//       	for(int i = 0; i < Temppath.size()-1;i++){
//       	  vertex_t cell1 = Temppath[i].cell1,cell2 = Temppath[i].cell2;
//       	if(verbose)cout << "Gpcopy[Temppath].cell1 is " << (Gpcopy)[cell1].IsletCellNum << " and .cell2 is " << (Gpcopy)[cell2].IsletCellNum << endl;
//       	if(verbose)cout << "Gp[Temppath].cell1 is " << (*Gp)[cell1].IsletCellNum << " and .cell2 is " << (*Gp)[cell2].IsletCellNum << endl;
//       	if(Gpcopy[cell1].IsletCellNum != (*Gp)[Temppath[i].cell1].IsletCellNum){
//       	  for(vp = vertices(Gpcopy);vp.first!=vp.second;vp.first++){
//       	    if(Gpcopy[*vp.first].IsletCellNum == (*Gp)[Temppath[i].cell1].IsletCellNum){
//       	      cell1 = *vp.first;
//       	      break;
//       	    }
//       	  }
//       	  if(Gpcopy[cell1].IsletCellNum != (*Gp)[Temppath[i].cell1].IsletCellNum)cout << "Problem: cell1 is not in Gpcopy" << endl;
//       	}
//       	if(Gpcopy[cell2].IsletCellNum != (*Gp)[Temppath[i].cell2].IsletCellNum){
//       	  for(vp = vertices(Gpcopy);vp.first!=vp.second;vp.first++){
//       	    if(Gpcopy[*vp.first].IsletCellNum == (*Gp)[Temppath[i].cell2].IsletCellNum){
//       	      cell2 = *vp.first;
//       	      break;
//       	    }
//       	  }
//       	  if(Gpcopy[cell2].IsletCellNum != (*Gp)[Temppath[i].cell2].IsletCellNum)cout << "Problem: cell2 is not in Gpcopy" << endl;
//       	}
//       	vector<vertex_t> tpath;vectorcount++;if(verbosemem)cout << "Adding tpath: vectorcount = " << vectorcount << endl; cout.flush();
//       	if(verbose)cout << "Finding shortest path between vertices " << (Gpcopy)[cell1].IsletCellNum << " and " << (Gpcopy)[cell2].IsletCellNum << endl;cout.flush();
       
//       	bool nopath = ShortestWeightedPathBetweenTwoVertices(&Gpcopy, cell1, cell2, &tpath);
//       	if(verbose)cout << "Finished ShortestWeightedPath: tpath.size() is" << tpath.size() <<  endl;cout.flush();
//       	double pathdistance = 0;
//       	for(int p = 0; p < tpath.size()-1;p++)
//       	  pathdistance+=Distancexy(Gpcopy[tpath[p]].x,Gpcopy[tpath[p]].y,
//       				   Gpcopy[tpath[p+1]].x,Gpcopy[tpath[p+1]].y);
//       	Temppath[i].distance = pathdistance;
//       	Temppath[i].path.insert(Temppath[i].path.begin(),tpath.begin(),tpath.end());
//       	RemoveVector(&tpath);vectorcount--;if(verbosemem)cout << "Removing tpath: vectorcount = " << vectorcount << endl; cout.flush();
//       }
//       if(verbose){
//       	cout << "After adding all IPparents and determining distances, Temppath contains: " << endl;
//       	for(int i = 0; i < Temppath.size();i++)
//       	  cout << (*Gp)[Temppath[i].cell1].IsletCellNum << "\t" << (*Gp)[Temppath[i].cell2].IsletCellNum << "\t" << Temppath[i].distance << endl; 
//       }
//       for(int i = 0; i < Temppath.size()-1;i++){
//       	// if(verbose){
//       	//   std::pair<out_edge_iter, out_edge_iter> oep;
//       	//   cout << "Determining distance from cell " << Gpcopy[Temppath[i].cell1].IsletCellNum <<" which has degree " << out_degree(Temppath[i].cell1,Gpcopy) <<" and edges " ;
//       	//   for(oep=out_edges(Temppath[i].cell1,Gpcopy);oep.first!=oep.second;oep.first++)cout << Gpcopy[target(*oep.first,Gpcopy)].IsletCellNum << " ";
//       	//   cout << endl;
//       	//   cout << "vectorcount = " << vectorcount << " and pathcount" << pathcount << endl;
//       	// }
	
//       	for(int j = i+2; j < Temppath.size();j++){
//       	  Path temp2; pathcount++;if(verbosemem)cout << "Adding temp2(Path): pathcount = " << pathcount << endl; cout.flush();
//       	  temp2.cell1 = Temppath[i].cell1;
//       	  temp2.cell2 = Temppath[j].cell1;
//       	  vector<vertex_t> tpath2;vectorcount++;if(verbosemem)cout << "Adding tpath2: vectorcount = " << vectorcount << endl; cout.flush();
//       	  vertex_t cell1 = Temppath[i].cell1,cell2 = Temppath[j].cell1;
//       	  if(verbose)cout << "Gpcopy[Temppath].cell1 is " << (Gpcopy)[cell1].IsletCellNum << " and .cell2 is " << (Gpcopy)[cell2].IsletCellNum << endl;
//       	  if(Gpcopy[cell1].IsletCellNum != (*Gp)[Temppath[i].cell1].IsletCellNum){
//       	    for(vp = vertices(Gpcopy);vp.first!=vp.second;vp.first++){
//       	      if(Gpcopy[*vp.first].IsletCellNum == (*Gp)[Temppath[i].cell1].IsletCellNum){
//       		cell1 = *vp.first;
//       		break;
//       	      }
//       	    }
//       	    if(Gpcopy[cell1].IsletCellNum != (*Gp)[Temppath[i].cell1].IsletCellNum)cout << "Problem: cell1 is not in Gpcopy" << endl;
//       	  }
//       	  if(Gpcopy[cell2].IsletCellNum != (*Gp)[Temppath[j].cell1].IsletCellNum ){
//       	    for(vp = vertices(Gpcopy);vp.first!=vp.second;vp.first++){
//       	      if(Gpcopy[*vp.first].IsletCellNum == (*Gp)[Temppath[i].cell2].IsletCellNum){
//       		cell2 = *vp.first;
//       		break;
//       	      }
//       	    }
//       	    if(Gpcopy[cell2].IsletCellNum != (*Gp)[Temppath[j].cell1].IsletCellNum)cout << "Problem: cell2 is not in Gpcopy" << endl;
//       	  }
//       	  bool nopathexists = ShortestWeightedPathBetweenTwoVertices(&Gpcopy, cell1, cell2, &tpath2);
//       	  double pathdistance = 0;
//       	  for(int p = 0; p < tpath2.size()-1;p++)
//       	  pathdistance+=Distancexy(Gpcopy[tpath2[p]].x,Gpcopy[tpath2[p]].y,
//       				   Gpcopy[tpath2[p+1]].x,Gpcopy[tpath2[p+1]].y);
//       	  if(verbose)cout << "the pathdistance between " << Gpcopy[temp2.cell1].IsletCellNum << " and " << Gpcopy[temp2.cell2].IsletCellNum << " is " << pathdistance << endl; 
//       	  temp2.distance=pathdistance;
//       	  //Update total distance if new path added
//       	  temp2.totalDistance = pathdistance;
//       	  for(int ii = 0;ii < i;ii++){
//       	    if(verbose) cout << "ii = " << ii << ": adding " << Temppath[ii].distance << endl; cout.flush();
//       	    temp2.totalDistance += Temppath[ii].distance;
//       	  }
//       	  for(int jj = j; jj < Temppath.size()-1;jj++){
//       	    if(verbose) cout << "jj = " << jj << ": adding " << Temppath[jj].distance << endl; cout.flush();
//       	    temp2.totalDistance += Temppath[jj].distance;
//       	  }
//       	  temp2.path.insert(temp2.path.begin(),tpath2.begin(),tpath2.end());
//       	  ShortestPaths.push_back(temp2);
//       	  //Clear temp2.path
//       	  RemovePath(&temp2);pathcount--;if(verbosemem)cout << "Removing temp2(Path): pathcount = " << pathcount << endl; cout.flush();
//       	  RemoveVector(&tpath2);vectorcount--;if(verbosemem)cout << "Removing  tpath2: vectorcount = " << vectorcount << endl; cout.flush();
//       	}
//       }
      
//       if(verbose){
//       	cout << "After adding shortest paths, Temppath contains: " << endl;
//       	for(int ii = 0; ii < Temppath.size();ii++){
//       	  cout << (*Gp)[Temppath[ii].cell1].IsletCellNum << "\t" << (*Gp)[Temppath[ii].cell2].IsletCellNum << "\t" << Temppath[ii].distance << "\t";
//       	  for(int p = 0;p < Temppath[ii].path.size();p++)
//       	    cout << (*Gp)[Temppath[ii].path[p]].IsletCellNum << ",";
//       	  cout << endl;
//       	}

//       }
      
     
//       while(ShortestPaths.size() > 0){
//       	//sort shortest paths by distance
//       	sort(ShortestPaths.begin(),ShortestPaths.end(),compareByPathtotalDistance);
//       	if(verbose){
//       	  cout << "And sorted shortest path  contains: " << endl;
//       	  for(int ii = 0; ii < ShortestPaths.size();ii++){
//       	    cout << (*Gp)[ShortestPaths[ii].cell1].IsletCellNum << "\t" << (*Gp)[ShortestPaths[ii].cell2].IsletCellNum << "\t" << ShortestPaths[ii].distance << "\t"<< ShortestPaths[ii].totalDistance << "\t";
//       	    for(int p = 0;p < ShortestPaths[ii].path.size();p++)
//       	      cout << (*Gp)[ShortestPaths[ii].path[p]].IsletCellNum << ",";
//       	    cout << endl;
//       	  }
//       	}
//       	// if(((*Gp)[ShortestPaths[0].cell1].CellNum != (*Gp)[endpt1].CellNum ||(*Gp)[ShortestPaths[0].cell2].CellNum != (*Gp)[endpt2].CellNum) &&
//       	//    ((*Gp)[ShortestPaths[0].cell2].CellNum != (*Gp)[endpt1].CellNum ||(*Gp)[ShortestPaths[0].cell1].CellNum != (*Gp)[endpt2].CellNum)){
//       	//Check if ShortestPath distance is less than sum of nearest neighbor distances
//       	double neighdist = 0;
//       	int Cell1pos=0,Cell2pos=0;
//       	for(int c1 = 0; c1 < Temppath.size()-1;c1++){
//       	  neighdist+=Temppath[c1].distance;
//       	  //Find Cell1 and cell2 in Temppath
//       	  if(ShortestPaths[0].cell1==Temppath[c1].cell1)Cell1pos = c1;
//       	  if(ShortestPaths[0].cell2==Temppath[c1].cell2)Cell2pos = c1;
//       	}
	  
	 
//       	//for(int c1 = 0; c1 < Temppath.size();c1++){
//       	// if((*Gp)[ShortestPaths[0].cell1].CellNum == (*Gp)[Temppath[c1].cell1].CellNum){
//       	//     Cell1pos = c1;
//       	//     neighdist += Temppath[c1].distance;
//       	//     break;
//       	//   }
//       	// }
//       	// int pos = Cell1pos;
//       	// while((*Gp)[Temppath[pos].cell2].CellNum!= (*Gp)[ShortestPaths[0].cell2].CellNum){
//       	//   pos++;
//       	//   neighdist += Temppath[pos].distance;
//       	// }
	
//       	if(verbose)cout << "For path between cells " <<(*Gp)[ShortestPaths[0].cell1].IsletCellNum  << " and " <<(*Gp)[ShortestPaths[0].cell2].IsletCellNum << " : neighdist is " << neighdist << endl;
	  
//       	if(ShortestPaths[0].totalDistance < neighdist){
//       	  if(verbose) cout << "shortest path is less than sum of neighbor distances" << endl;
//       	  if(verbose){
//       	    cout << "endpt1 is cell " << (*Gp)[endpt1].IsletCellNum << endl;
//       	    cout << "cell1pos is cell " << (*Gp)[Temppath[Cell1pos].cell1].IsletCellNum << " at position " << Cell1pos << endl;
//       	    cout << "cell2pos is cell " << (*Gp)[Temppath[Cell2pos].cell2].IsletCellNum << " at position " << Cell2pos << endl;
//       	    cout << "endpt2 is cell " << (*Gp)[endpt2].IsletCellNum << endl;
//       	  }
//       	  bool replacePath = 1;
//       	  if(isSetLoop){
//       	    if(verbose)cout << "Checking if new path surrounds set" << endl;
//       	    //Create new possible path
//       	    int tpi = 0;
//       	    vector<vertex_t> PossiblePath;vectorcount++;if(verbosemem)cout << "Adding PossiblePath: vectorcount = " << vectorcount << endl; cout.flush();
//       	    //Add beginning of original path to Possible path
//       	    for(int i = 0; i < endpt1pos;i++)PossiblePath.push_back((*path)[i]);
//       	    if(verbose){
//       	      cout << "After adding beginning:New Possible path is ";
//       	      for(int p = 0; p < PossiblePath.size();p++)cout << (*Gp)[PossiblePath[p]].IsletCellNum << " ";
//       	      cout << endl;
//       	    }
//       	    //before Shortestpath
//       	    while(Temppath[tpi].cell1 !=ShortestPaths[0].cell1){
//       	      PossiblePath.push_back(Temppath[tpi].cell1);
//       	      for(int p = 1; p < Temppath[tpi].path.size()-1;p++)PossiblePath.push_back(Temppath[tpi].path[p]);
//       	      tpi++;
//       	    }
//       	    if(verbose){
//       	      cout << "BeforeShortest:New Possible path is ";
//       	      for(int p = 0; p < PossiblePath.size();p++)cout << (*Gp)[PossiblePath[p]].IsletCellNum << " ";
//       	      cout << endl;
//       	    }
//       	    //Shortestpath
//       	    for(int p = 0;p < ShortestPaths[0].path.size();p++)PossiblePath.push_back(ShortestPaths[0].path[p]);
//       	    if(verbose){
//       	      cout << "AfterShortest:New Possible path is ";
//       	      for(int p = 0; p < PossiblePath.size();p++)cout << (*Gp)[PossiblePath[p]].IsletCellNum << " ";
//       	      cout << endl;
//       	    }
//       	    //rest of Temppath
//       	    while(Temppath[tpi].cell1!=ShortestPaths[0].cell2)tpi++;
//       	    //tpi++;
//       	    if(verbose)cout << "After updating:tpi =" << tpi << endl;
//       	    for(int tp = tpi; tp < Temppath.size();tp++){
//       	      PossiblePath.push_back(Temppath[tp].cell1);
//       	      if(tp < Temppath.size()-1)
//       		for(int p = 1; p < Temppath[tp].path.size()-1;p++)PossiblePath.push_back(Temppath[tp].path[p]);
	      
//       	    }
//       	    //Add rest of original path
//       	    for(int i = endpt2pos+1; i < (*path).size();i++)PossiblePath.push_back((*path)[i]);
//       	    if(PossiblePath[PossiblePath.size()-1]!=PossiblePath[0])PossiblePath.push_back(PossiblePath[0]);
//       	    if(verbose){
//       	      cout << "New Possible path is ";
//       	      for(int p = 0; p < PossiblePath.size();p++)cout << (*Gp)[PossiblePath[p]].IsletCellNum << " ";
//       	      cout << endl;
//       	    }
//       	    //Check if new path contains set
//       	    for(int s = 0; s < (*set).size();s++){
//       	      int wn = WindingNumber(Gp,&PossiblePath,(*set)[s][0],(*set)[s][1]);
//       	      if(verbose)cout << "For point (" << (*set)[s][0] << "," << (*set)[s][1] << ") the wn is " << wn << endl;
//       	      if(wn == 0){
//       		replacePath = 0;
//       		break;
//       	      }
//       	    }
//       	    RemoveVector(&PossiblePath);vectorcount--;if(verbosemem)cout << "Removing PossiblePath: vectorcount = " << vectorcount << endl; cout.flush(); 
//       	  }
//       	  if(replacePath){
//       	    //Determine which cells were removed
//       	    vector<vertex_t> CellsRemoved;vectorcount++;if(verbosemem)cout << "Adding CellsRemoved: vectorcount = " << vectorcount << endl; cout.flush();
//       	    for(int p = Cell1pos+1;p <Cell2pos+1;p++)CellsRemoved.push_back(Temppath[p].cell1);
//       	    if(verbose){
//       	      cout << "Cells being removed are " ;
//       	      for(int p = 0; p < CellsRemoved.size();p++)cout << (*Gp)[CellsRemoved[p]].IsletCellNum << " ";
//       	      cout << endl;
//       	    }
//       	    //Remove from temppath
//       	    for(int p = Cell1pos;p < Cell2pos+1;p++){
//       	      RemovePath(&Temppath[p]);
//       	    }
//       	    Temppath.erase(Temppath.begin()+Cell1pos,Temppath.begin()+Cell2pos+1);
//       	    //Replace with shortest path
//       	    Temppath.insert(Temppath.begin()+Cell1pos,ShortestPaths[0]);
//       	    //Clear out shortest paths containing cells removed
//       	    for(int p = 0; p < CellsRemoved.size();p++){
//       	      for(int p1 = 0;p1 < ShortestPaths.size();p1++){
//       		if((*Gp)[CellsRemoved[p]].CellNum ==(*Gp)[ShortestPaths[p1].cell1].CellNum ||
//       		   (*Gp)[CellsRemoved[p]].CellNum ==(*Gp)[ShortestPaths[p1].cell2].CellNum){
//       		  RemovePath(&ShortestPaths[p1]);
//       		  ShortestPaths.erase(ShortestPaths.begin()+p1);
//       		  p1--;
//       		}
//       	      }
//       	    }
//       	    if(verbose){
//       	      cout << "After removing paths, Temppath contains: " << endl;
//       	      for(int j = 0; j < Temppath.size();j++){
//       		cout << (*Gp)[Temppath[j].cell1].IsletCellNum << "\t" << (*Gp)[Temppath[j].cell2].IsletCellNum << "\t" << Temppath[j].distance << "\t";
//       		for(int p = 0;p < Temppath[j].path.size();p++)
//       		  cout << (*Gp)[Temppath[j].path[p]].IsletCellNum << ",";
//       		cout << endl;
//       	      }
//       	      cout << "And shortest path  contains: " << endl;
//       	      for(int j = 0; j < ShortestPaths.size();j++){
//       		cout << (*Gp)[ShortestPaths[j].cell1].IsletCellNum << "\t" << (*Gp)[ShortestPaths[j].cell2].IsletCellNum << "\t" << ShortestPaths[j].distance << "\t"<< ShortestPaths[j].totalDistance << "\t";
//       		for(int p = 0;p < ShortestPaths[j].path.size();p++)
//       		  cout << (*Gp)[ShortestPaths[j].path[p]].IsletCellNum << ",";
//       		cout << endl;
//       	      }
//       	    }


//       	    //Update Shortest Path Total distances
//       	    for(int sp = 0; sp < ShortestPaths.size();sp++){
//       	      //Update total distance if new path added
//       	      ShortestPaths[sp].totalDistance = ShortestPaths[sp].distance;
//       	      int tp = 0;
//       	      while(Temppath[tp].cell1!=ShortestPaths[sp].cell1){
//       		ShortestPaths[sp].totalDistance +=Temppath[tp].distance;
//       		tp++;
//       	      }
//       	      while(Temppath[tp].cell1!=ShortestPaths[sp].cell2)tp++;
//       	      while(tp < Temppath.size()-1){
//       		ShortestPaths[sp].totalDistance +=Temppath[tp].distance;
//       		tp++;
//       	      }
//       	    }
	    
	    
//       	    RemoveVector(&CellsRemoved);vectorcount--;if(verbosemem)cout << "Removing CellsRemoved: vectorcount = " << vectorcount << endl; cout.flush();
//       	  }
//       	  if(verbose){
//       	    cout << "After removing paths, Temppath contains: " << endl;
//       	    for(int j = 0; j < Temppath.size();j++){
//       	      cout << (*Gp)[Temppath[j].cell1].IsletCellNum << "\t" << (*Gp)[Temppath[j].cell2].IsletCellNum << "\t" << Temppath[j].distance << "\t";
//       	      for(int p = 0;p < Temppath[j].path.size();p++)
//       		cout << (*Gp)[Temppath[j].path[p]].IsletCellNum << ",";
//       	      cout << endl;
//       	    }
//       	    cout << "And shortest path  contains: " << endl;
//       	    for(int j = 0; j < ShortestPaths.size();j++){
//       	      cout << (*Gp)[ShortestPaths[j].cell1].IsletCellNum << "\t" << (*Gp)[ShortestPaths[j].cell2].IsletCellNum << "\t" << ShortestPaths[j].distance << "\t"<< ShortestPaths[j].totalDistance << "\t";
//       	      for(int p = 0;p < ShortestPaths[j].path.size();p++)
//       		cout << (*Gp)[ShortestPaths[j].path[p]].IsletCellNum << ",";
//       	      cout << endl;
//       	    }
//       	  }
//       	}
//       	//Remove the first shortest path
//       	RemovePath(&ShortestPaths[0]);
//       	ShortestPaths.erase(ShortestPaths.begin());
	
//       	if(verbose){
//       	  cout << "After path shortest path  contains: " << endl;
//       	  for(int i = 0; i < ShortestPaths.size();i++){
//       	    cout << (*Gp)[ShortestPaths[i].cell1].IsletCellNum << "\t" << (*Gp)[ShortestPaths[i].cell2].IsletCellNum << "\t" << ShortestPaths[i].distance << "\t" << ShortestPaths[i].totalDistance << "\t";
//       	    for(int p = 0;p < ShortestPaths[i].path.size();p++)
//       	      cout << (*Gp)[ShortestPaths[i].path[p]].IsletCellNum << ",";
//       	    cout << endl;
//       	  }
//       	}
//       }
//       RemoveVector(&ShortestPaths);vectorcount--;if(verbosemem)cout << "Removing ShortestPaths: vectorcount = " << vectorcount << endl; cout.flush();
        
//       //Remove IPs
//       (*path).erase((*path).begin()+endpt1pos+1,(*path).begin()+endpt2pos);
//       //insert new path
//       vector<vertex_t> newpath;vectorcount++;if(verbosemem)cout << "Adding newpath: vectorcount = " << vectorcount << endl; cout.flush();
//       if(numNonIPs ==0)newpath.push_back(Temppath[0].cell1);
//       for(int p = 0; p < Temppath.size()-2;p++)
//       	newpath.insert(newpath.end(),Temppath[p].path.begin()+1,Temppath[p].path.end());
//       if(Temppath[Temppath.size()-2].path.size() > 2)
//       	newpath.insert(newpath.end(),Temppath[Temppath.size()-2].path.begin()+1,Temppath[Temppath.size()-2].path.end()-1);
//       if(numNonIPs <2)newpath.push_back(Temppath[Temppath.size()-1].cell1);
//       if(numNonIPs == 0)newpath.push_back(Temppath[0].cell1);
//       (*path).insert((*path).begin()+endpt1pos+1,newpath.begin(),newpath.end());
//       //Remove newpath
//       RemoveVector(&newpath);vectorcount--;if(verbosemem)cout << "Removing newpath: vectorcount = " << vectorcount << endl; cout.flush();
//       if(verbose){
//       	cout << "Path now contains: " ;
//       	for(int i = 0; i < (*path).size(); i++)cout << (*Gp)[(*path)[i]].IsletCellNum << " ";
//       	cout << endl;
//       }
      
//       //Clear out Temppath
//       RemoveVector(&Temppath);vectorcount--;if(verbosemem)cout << "Removing Temppath: vectorcount = " << vectorcount << endl; cout.flush();
      
//     }
    
//     //check if still IPs
//     ContainsIPs = 0;
//     for(int i = 0; i < (*path).size();i++){
//       if((*Gp)[(*path)[i]].CellNum > 99999){
//     	bool leave = 0;
//     	for(int j = 0; j < (*LeaveIP).size();j++){
//     	  if((*path)[i] == (*LeaveIP)[j]){
//     	    leave = 1;
//     	    break;
//     	  }
//     	}
//     	if(!leave){
//     	  ContainsIPs = 1; 
//     	  break;
//     	}
//       }
//     }
//     RemoveGraph(&Gpcopy2); graphcount--;if(verbosemem)cout << "Removing Gpcopy2 (graph): graphcount = " << graphcount << endl; cout.flush();
//     // iter++;//**********************Make sure to remove******************//
//     // if(iter==12)ContainsIPs = 0;//**********************Make sure to remove******************//
//   }
//   //Remove Graph
//   RemoveGraph(&Gpcopy);graphcount--;if(verbosemem)cout << "Removing Gpcopy (graph): graphcount = " << graphcount << endl; cout.flush();
//  if(verbose){
//    cout << "After CleanIPs: Path now contains: " ;
//    for(int i = 0; i < (*path).size(); i++)cout << (*Gp)[(*path)[i]].IsletCellNum << " ";
//    cout << endl;
//  }

//   if(vectorcount != 0 || pathcount != 0 || graphcount!=0) cout << " ***Problem with vectors and/or paths and/or graphs! vectorcount = " << vectorcount << " and pathcount = " << pathcount << " and graphcount = " << graphcount << " in CleanIntersectionPoints" << endl;
//  if(verbose)
//    cout << "**************End of CleanIntersectionPointsInPaths (graphcount = " << graphcount << " vectorcount = " << vectorcount << " and pathcount = " << pathcount << " **********************" << endl;cout.flush();

//   return 0;
// }

// //*********************************************************************************************//
// //*********************************************************************************************//
// //*********************************************************************************************//
// void CreateMaskOfIslet(Graph* G,vector<vector<bool> > *Imask){
//   bool verbose = 0;
  
//   //Find max,min x and max,min y
//   double xmax = -100000, xmin = 100000, ymax = -100000, ymin = 100000;
//   for(int i = 0; i < num_vertices(*G);i++){
//     if((*G)[i].x > xmax)xmax = (*G)[i].x;
//     if((*G)[i].x < xmin)xmin = (*G)[i].x;
//     if((*G)[i].y > ymax)ymax = (*G)[i].y;
//     if((*G)[i].y < ymin)ymin = (*G)[i].y;
//   }
  
//   vector<vector<bool> > Mask((int)(xmax-xmin+1),vector<bool>((int)(ymax-ymin+1),0));
  
//   for(int i = 0; i < num_vertices(*G);i++){
//     for(int j= -10; j < 11;j++){
//       if((int)(*G)[i].x -(int)xmin -j > 0 && (int)(*G)[i].x -(int)xmin -j < Mask.size()){
// 	for(int k = -11; k < 10;k++){
// 	  if((int)(*G)[i].y -(int)ymin -k > 0 && (int)(*G)[i].y -(int)ymin -k < Mask[j].size() )
// 	    Mask[(int)(*G)[i].x -(int)xmin -j][(int)(*G)[i].y -(int)ymin -k] = 1;
// 	}
//       }
//     }
//   }
  
//   Imask = & Mask;
//   return;
// }
// //*********************************************************************************************//
// //*********************************************************************************************//
// //*********************************************************************************************//

// void FindHullPath(IsletGraph * I){
//   int vectorcount = 0;
//   bool verbose=0,verbosemem=0;
//   std::pair <vertex_iter,vertex_iter> vp,vp2;
//   Graph * G = &((*I).allGraph);
//   if(verbose) cout << "inside FindHullPath" << endl; cout.flush();
//   //Find furthest distance between points
//   double maxdist = 0;
//   vertex_t jmax , kmax;
//   for(vp = vertices(*G);vp.first!=vp.second;vp.first++){
//     double x1 = (*G)[*vp.first].x, y1 = (*G)[*vp.first].y;
//     vp2=vertices(*G);
//     for(vp2.first = vp.first+1; vp2.first!=vp2.second;vp2.first++){
//       double x2 = (*G)[*vp2.first].x, y2 = (*G)[*vp2.first].y;
//       double radius = Distancexy(x1,y1,x2,y2);
//       if(radius > maxdist){
// 	maxdist = radius;
// 	jmax = *vp.first;
// 	kmax = *vp2.first;
//       }
//     }
//   }
//   if(verbose)cout << "The furthest distance is between cells " << (*G)[jmax].IsletCellNum << " and "<< (*G)[kmax].IsletCellNum << endl; cout.flush();
//   //find furthest point from line 
//   double x1 = (*G)[jmax].x, y1 = (*G)[jmax].y;
//   double x2 = (*G)[kmax].x, y2 = (*G)[kmax].y;
//   maxdist = 0;
//   vertex_t pmax;
//   for(vp = vertices(*G);vp.first!=vp.second;vp.first++){
//     //if(verbose)cout << "checking cell " << (*G) 
//     double x0 = (*G)[*vp.first].x, y0 = (*G)[*vp.first].y;
//     double distance = DistanceFromPointToLineSegment(x0,y0,x1,y1,x2,y2);
//     if(distance > maxdist){
//       maxdist = distance;
//       pmax = *vp.first;
//       //if(verbose)cout << "updating pmax to cell " << (*G)[pmax].IsletCellNum << " with distance " << maxdist <<endl;cout.flush();
//     }
//   }
  
//   //Create original path
//   (*I).hullPath.push_back(jmax);
//   (*I).hullPath.push_back(pmax);
//   (*I).hullPath.push_back(kmax);
//   (*I).hullPath.push_back(jmax);
//   int fignum = 0;
//   if(verbose)PrintoutGnuplotHull3x4(I,fignum,"temp");
//   //Create vector of points not within loop
//   vector<ParentEdge> PEdges;vectorcount++;
//   for(vp = vertices(*G);vp.first!=vp.second;vp.first++){
//     //Check if vertex is in loop
//     bool inloop = 0;
//     for(int i = 0; i < (*I).hullPath.size();i++){
//       if(*vp.first == (*I).hullPath[i]){
// 	inloop = 1;
// 	break;
//       }
//     }
//     if(!inloop){
//       //Calculate winding number for each point to see if located inside the loop
//       int wn = WindingNumber(&((*I).allGraph),&((*I).hullPath),(*G)[*vp.first].x, (*G)[*vp.first].y);
//       if(wn == 0){//point is outside of the loop
// 	ParentEdge PEdge;
// 	PEdge.v = *vp.first;
// 	PEdges.push_back(PEdge);
//       }
//     }
//   }
//   if(verbose)cout << "There are " << PEdges.size() << " vertices outside the loop" << endl;
//   //For each vertex outside loop,Find shortest distance to an edge 
//   maxdist = 0;
//   int Pedgemax = 0;
//   for(int i = 0; i < PEdges.size();i++){
//     double mindist = 100000;
//     vertex_t ev1, ev2;
//     for(int j = 0; j < (*I).hullPath.size()-1; j++){
//       double distance = DistanceFromPointToLineSegment((*I).allGraph[PEdges[i].v].x, (*I).allGraph[PEdges[i].v].y, (*I).allGraph[(*I).hullPath[j]].x,(*I).allGraph[(*I).hullPath[j]].y,(*I).allGraph[(*I).hullPath[j+1]].x,(*I).allGraph[(*I).hullPath[j+1]].y);
//       // if(verbose)cout << "distance for cell " << (*G)[PEdges[i].v].IsletCellNum << " is " << distance << endl; 
//       if(distance < mindist){
// 	//if(verbose)cout << "distance is less than mindist  =" << mindist << endl;
// 	mindist = distance;
// 	PEdges[i].distance =mindist;
// 	PEdges[i].ev1 = (*I).hullPath[j];
// 	PEdges[i].ev2 = (*I).hullPath[j+1];
// 	if(mindist > maxdist && mindist != 100000){
// 	  maxdist = mindist;
// 	  Pedgemax = i;
// 	}
//       }
//     }
//   }
//   if(verbose)cout << "pedgemax is cell " << (*G)[PEdges[Pedgemax].v].IsletCellNum << " with closest edge " << (*G)[PEdges[Pedgemax].ev1].IsletCellNum << " and " << (*G)[PEdges[Pedgemax].ev2].IsletCellNum << endl;
//   bool allinloop = 0;
//   if(PEdges.size()== 0)allinloop= 1;
//   while (!allinloop){
//     if(verbose)cout << "inside while loop" << endl;cout.flush();
//     fignum++;
//     //add maximal Pedge to path
//     int iv1 = 0;
//     while((*I).hullPath[iv1] != PEdges[Pedgemax].ev1)iv1++;
//     if(verbose){
//       cout << "Before updating path (iv1 = " << iv1 << "): The new path is " << endl;
//       for (int p = 0; p < (*I).hullPath.size();p++)
// 	cout << (*G)[(*I).hullPath[p]].IsletCellNum << " " ;
//     }
//     (*I).hullPath.insert((*I).hullPath.begin()+iv1+1,PEdges[Pedgemax].v);
//     if(verbose){
//       cout << "Updated path: The new path is " << endl;
//       for (int p = 0; p < (*I).hullPath.size();p++)
// 	cout << (*G)[(*I).hullPath[p]].IsletCellNum << " " ;
//     }
//     double maxdist = 0;
//     vertex_t Pedgemaxvold = PEdges[Pedgemax].v, Pedgemaxv = 0;
//     if(verbose)cout << "pedgemaxvold is cell " << (*G)[Pedgemaxvold].IsletCellNum  << endl;
//     //check vertices that were closest to old edge
//     vertex_t ev1max = PEdges[Pedgemax].ev1;
//     for(int i = PEdges.size()-1; i > -1;i--){
//       cout << "Cell " << (*G)[PEdges[i].v].IsletCellNum << " is closest to line " <<  (*G)[PEdges[i].ev1].IsletCellNum << " and " << (*G)[PEdges[i].ev2].IsletCellNum << endl;
//       if(PEdges[i].ev1 == ev1max  && i !=Pedgemax){
// 	if(verbose)cout << "inside update" << endl; 
// 	//check winding number
// 	int wn = WindingNumber(&((*I).allGraph),&((*I).hullPath),(*G)[PEdges[i].v].x, (*G)[PEdges[i].v].y);
// 	if (wn != 0){ //point is inside loop
// 	  //Remove entry from list
// 	  if(verbose)cout << "Removing cell " << (*G)[PEdges[i].v].IsletCellNum << endl;
// 	  PEdges.erase(PEdges.begin()+i);
	  
// 	}
// 	else{//still outside need to update PEdge
// 	  double Dist_e1 = DistanceFromPointToLineSegment((*I).allGraph[PEdges[i].v].x, (*I).allGraph[PEdges[i].v].y, (*I).allGraph[(*I).hullPath[iv1]].x,(*I).allGraph[(*I).hullPath[iv1]].y,(*I).allGraph[(*I).hullPath[iv1+1]].x,(*I).allGraph[(*I).hullPath[iv1+1]].y);
// 	  double Dist_e2 = DistanceFromPointToLineSegment((*I).allGraph[PEdges[i].v].x, (*I).allGraph[PEdges[i].v].y, (*I).allGraph[(*I).hullPath[iv1+1]].x,(*I).allGraph[(*I).hullPath[iv1+1]].y,(*I).allGraph[(*I).hullPath[iv1+2]].x,(*I).allGraph[(*I).hullPath[iv1+2]].y);
// 	  if(Dist_e1 < Dist_e2){
// 	    PEdges[i].ev2 = (*I).hullPath[iv1+1];
// 	    PEdges[i].distance = Dist_e1;
// 	  }
// 	  else{
// 	    PEdges[i].ev1 = (*I).hullPath[iv1+1];
// 	    PEdges[i].distance = Dist_e2;
// 	  }
// 	  cout << "Now, Cell " << (*G)[PEdges[i].v].IsletCellNum << " is closest to line " <<  (*G)[PEdges[i].ev1].IsletCellNum << " and " << (*G)[PEdges[i].ev2].IsletCellNum << endl;
// 	} 
//       }
//       //update Pedgemaxnew
//       if(PEdges[i].v != Pedgemaxvold && i < PEdges.size()){
// 	if(PEdges[i].distance > maxdist){
// 	  maxdist = PEdges[i].distance;
// 	  Pedgemaxv = PEdges[i].v;
// 	}
//       }
//     }
    
//     //Remove Pedgemaxold
//     int pv1 = 0; 
//     while (PEdges[pv1].v != Pedgemaxvold)pv1++;
//     if(verbose)cout << "Erasing pedgemax at Cell " << (*G)[PEdges[pv1].v].IsletCellNum << endl; cout.flush();
//     PEdges.erase(PEdges.begin()+pv1);
    
//     if(fignum < 12) PrintoutGnuplotHull3x4(I,fignum,"temp");
//     if(PEdges.size() == 0)allinloop = 1;
//     else{
//       int pv1 = 0; 
//       while (PEdges[pv1].v != Pedgemaxv)pv1++;
//       Pedgemax = pv1;
//       if(verbose)cout << "Updating pedgemax to " << (*G)[Pedgemax].IsletCellNum << endl; cout.flush();
//     }
//   }
//   cout << "Finished finding hull: fignum is " << fignum << endl; 
//   return;
// }   



// //*********************************************************************************************//
// //*********************************************************************************************//
// //*********************************************************************************************//


// //*********************************************************************************************//
// //*********************************************************************************************//
// //*********************************************************************************************//

// void FindNextCellInPath(IsletGraph* Ip,Graph*Gp, vector<vertex_t>*Path, double* MidAngle, int pathnum,vertex_t* nextv){
//   bool verbose=0,verbosemem=0;
//   int vectorcount = 0;
//   std::pair<out_edge_iter, out_edge_iter> oep;
//   vertex_t u = (*Path)[(*Path).size()-1];
//   if(verbose)cout << " For vertex " << (*Gp)[u].IsletCellNum << " with degree " << out_degree(u,*Gp) << endl;cout.flush();
//   MakeGraphLocallyPlanar(Ip,Gp,u);
//   if(verbose)cout << "Finished making " << (*Gp)[u].IsletCellNum << " locally planar" <<endl;cout.flush();
//   //find angle of each edge of u and sort
//   vector <CellAngle> CellAngles;vectorcount++;if(verbosemem)cout << "Adding CellAngles: vectorcount = " << vectorcount << endl; cout.flush();
//   for(oep = out_edges(u,(*Gp));oep.first!=oep.second; oep.first++){
//     CellAngle temp;
//     temp.Cell = target(*oep.first,(*Gp));
//     temp.Angle = GetAngleOfLine((*Gp)[u].x,(*Gp)[u].y, (*Gp)[target(*oep.first,(*Gp))].x,(*Gp)[target(*oep.first,(*Gp))].y);
//     CellAngles.push_back(temp);
//   }
//   sort(CellAngles.begin(),CellAngles.end(),compareByAngle);
//   if(verbose){
//     cout << "The sorted angles are " << endl;
//     for(int ca = 0; ca < CellAngles.size();ca++)
//       cout << (*Gp)[CellAngles[ca].Cell].IsletCellNum << "\t" << CellAngles[ca].Angle << endl;
//   }
  
//   //find previous leg of path in list and add next angle to path (if pathsize = 1 use 3pi/2)
//   int CApos = 0;
//   if(verbose)cout << "The path has " << (*Path).size() << " vertices and pathnum = " << pathnum << endl;
//   if((*Path).size() == 1){
//     if(pathnum!= 1){
//       for(int i = 0; i < CellAngles.size();i++){
// 	if(verbose)cout << "CellAngles[" << i << "].Cell is " <<(*Gp)[CellAngles[i].Cell].IsletCellNum <<" and nextv  = " <<  (*Gp)[*nextv].IsletCellNum << endl;
// 	if((*Gp)[CellAngles[i].Cell].IsletCellNum == (*Gp)[*nextv].IsletCellNum){
// 	  if(verbose)cout << "They are equal" << endl;
// 	  CApos = (i+1)%CellAngles.size();
// 	  break;
// 	}
//       }
//     }
//     else{
//       if(verbose)cout << "pathnum = 1 and the last Cell angle is " << CellAngles[CellAngles.size()-1].Angle << endl; 
//       if(CellAngles[CellAngles.size()-1].Angle > 3.0*PI/2.0){
// 	while(CellAngles[CApos].Angle <3.0*PI/2.0)CApos++;
//       }
//     }
//     (*Path).push_back(CellAngles[CApos].Cell);
//     (*nextv) = CellAngles[CApos].Cell;
//     if(verbose)cout << "Adding cell " << (*Gp)[(*Path)[(*Path).size()-1]].IsletCellNum << " to path " << endl;cout.flush();
//   }
//   else{
//     int newpos = 0;
//     if(CellAngles[CellAngles.size()-1].Cell != (*Path)[(*Path).size()-2]){
//       for(int ca = 0; ca < CellAngles.size();ca++){
// 	if((*Path)[(*Path).size()-2] == CellAngles[ca].Cell){
// 	  //   cout << "Previous path is from " << CellAngles[ca].Cell << " and next angle is at " << CellAngles[ca+1].Cell << endl;
// 	  newpos = ca+1;
	  
// 	  break;
// 	}
//       }
//     }
//     if(CellAngles.size() == 1)*MidAngle = mod(CellAngles[0].Angle+PI,2*PI);
//     else if(newpos==0)*MidAngle = mod((CellAngles[0].Angle+(2*PI - CellAngles[CellAngles.size()-1].Angle))/2.0 +CellAngles[CellAngles.size()-1].Angle,2*PI);
//     else *MidAngle = (CellAngles[newpos].Angle+CellAngles[newpos-1].Angle)/2.0;
//     if(verbose)cout << "MidAngle = " << *MidAngle <<  endl;
//     (*Path).push_back(CellAngles[newpos].Cell);
//   }
//   RemoveVector(&CellAngles);vectorcount--;if(verbosemem)cout << "Removing CellAngles: vectorcount = " << vectorcount << endl; cout.flush();

//   return;
// }
  


// //*********************************************************************************************//
// //*********************************************************************************************//
// //*********************************************************************************************//





// //*********************************************************************************************//
// //*********************************************************************************************//
// //*********************************************************************************************//


// void RemoveVerticesAwayFromSets(Graph * G,string settype, double Threshold){
//   bool verbose=0,verbosemem=0;
//   int vectorcount = 0;
//   vector<int> VerticesToRemove;vectorcount++;if(verbosemem)cout << "Adding vector VerticesToRemove: vectorcount = " << vectorcount << endl;cout.flush();

//   for(int i = 0; i < num_vertices(*G); i++){
//     if((settype == "ad" && (*G)[i].type =="b") ||
//        (settype == "b" && (*G)[i].type == "a") ||
//        (settype == "b" && (*G)[i].type == "d")){
//       bool WithinThreshold = 0;
//       for(int j = 0; j < num_vertices(*G);j++){
// 	if((settype == "b" &&(*G)[j].type == "b") ||
// 	   (settype == "ad" && (*G)[j].type == "a") ||
// 	   (settype == "ad" && (*G)[j].type == "d") ){
// 	  double dist = Distancexy((*G)[i].x, (*G)[i].y, (*G)[j].x, (*G)[j].y);
// 	  if(dist < Threshold){
// 	    WithinThreshold = 1;
// 	    j = num_vertices(*G);
// 	  }
// 	}
//       }
//       if(!WithinThreshold)VerticesToRemove.push_back(i);
//     }
//   }
//   sort(VerticesToRemove.begin(),VerticesToRemove.end());
//   //RemoveVertices
//   for(int i = VerticesToRemove.size()-1; i > -1; i--){
//     ClearVertex(G,VerticesToRemove[i]);
//     clear_vertex(VerticesToRemove[i],(*G));
//     remove_vertex(VerticesToRemove[i],(*G));
//   }
//   if(verbose)cout << "Removed " << VerticesToRemove.size() << " vertices" << endl;
//   RemoveVector(&VerticesToRemove);vectorcount--;if(verbosemem)cout << "Removing vector VerticesToRemove: vectorcount = " << vectorcount << endl;cout.flush();
//   if(vectorcount != 0) cout << "Problem with vectors! vectorcount = " << vectorcount << " in RemoveVerticesAwayFromSets" << endl;
//   return;
// }
// //*********************************************************************************************//
// //*********************************************************************************************//
// //*********************************************************************************************//


// void DeterminePlanarGraph(IsletGraph * I, Graph * Gp, Graph *PGp, string fprefix, string sthresh){

//   bool verbose=0;
//   if(verbose)cout << "Inside DeterminePlanarGraph..." << endl;cout.flush();
//   if(num_vertices(*Gp)>0){
//     //Names should be fprefix.Islet#.type.gr.verts
//     string sIslet =  static_cast<ostringstream*>( &(ostringstream() << (*I).IsletNum) )->str();
//     string type;
//     if((*Gp)[0].type == "b")type  = "b";
//     else type = "ad";
//     string vertfile = "./PlanarGraphs/"+fprefix+".Islet"+sIslet+"."+type+"."+sthresh+".IPvertices"; 
//     string edgefile = "./PlanarGraphs/"+fprefix+".Islet"+sIslet+"."+type+"."+sthresh+".edges";
//     bool vertfileexists = file_exists(vertfile.c_str());
//     bool edgefileexists = file_exists(edgefile.c_str());
    
    
//     //if(verbose)cout << "Finished copying graph vertices.  There are " << num_vertices(*PGp) << "vertices" <<  endl; cout.flush();
//     if(vertfileexists && edgefileexists){
//       if(verbose)cout << "Verts and edges files exist! " << endl;cout.flush();
//       CopyGraphVertices(Gp,PGp);
      
//       //Read in IP verts file
//       AddIPVertices(PGp,vertfile);
//       if(verbose)cout << "Finished adding IP vertices from " << vertfile<< endl; cout.flush();
//       //Read in edge file
//       AddIPEdges(PGp,edgefile);
//       if(verbose)cout << "Finished adding IP edges from " <<edgefile << endl; cout.flush();
//       //PrintoutGnuplotWithIPs(PGp,"ReadinPlanar."+type);
//     }
//     else{
//       CopyGraph(Gp,PGp);
//       if(verbose)cout << "Finished copying graph edges" << endl; cout.flush();
//       MakeGraphPlanar(I,PGp);
//       if(verbose)cout << "Finished making graph planar" << endl;cout.flush();
//       //PrintoutGnuplotWithIPs(PGp,"calculatedPlanar."+type);
//       if(num_vertices(*Gp) > 500){
// 	//Write out file
// 	PrintoutIPVertices(PGp,vertfile);
// 	if(verbose)cout << "Finished printing out vertices to " << vertfile <<  endl;cout.flush();
// 	PrintoutIPEdges(PGp,edgefile);
// 	if(verbose)cout << "Finished printing out edges to " << edgefile <<  endl;cout.flush();
//       }
//     }
//   }
//   return;
// }

// //*********************************************************************************************//
// //*********************************************************************************************//
// //*********************************************************************************************//

// void MakeGraphLocallyPlanar(IsletGraph * I,Graph * Gp, vertex_t u1){
//   std::pair<out_edge_iter, out_edge_iter> oep1,oep2;
//   std::pair<vertex_iter,vertex_iter>vp;
//   // std::pair<edge_iter, edge_iter> ep2;
//   vertex_t v1, u2, v2, u;
//   edge_t e;
//   //PrintoutGnuplot(Gp,"temp");
//   bool verbose=0,verbosemem=0;
//   if((*Gp)[u1].IsletCellNum == 43)verbose=0;
//   int num = 100000;
//   sort((*I).IPs.begin(),(*I).IPs.end(),compareByCellNum);
//   for(int i = 0; i < (*I).IPs.size();i++){
//     if((*I).IPs[i].CellNum == num)num++;
//     else break;
//   }
 
//   vector<vertex_t> *Cp;
//   Graph *RG;
//   if(verbose)  cout << "The compnum for u (Cell " << (*Gp)[u1].IsletCellNum << ") is " << (*Gp)[u1].CompNum << " with type " << (*Gp)[u1].type << endl;
//   if((*Gp)[u1].type == "b"){
//     Cp = &((*I).bComponents[(*Gp)[u1].CompNum]);
//     RG = &(*I).bGraph;
//   }
//   else{
//     Cp = &((*I).adComponents[(*Gp)[u1].CompNum]);
//     RG = &(*I).adGraph;
//   }
//   if(verbose){
//     cout << "Inside MakeGraphLocallyPlanar" << endl;
//     cout << "The component contains" << endl;
//     for(int i = 0;i < (*Cp).size();i++)
//       cout << (*RG)[(*Cp)[i]].IsletCellNum << " ";cout.flush();
//     cout << endl;
//   }
//   //PrintoutGnuplotWithIPs(Gp,"temp");
//   bool intersections = 0;
//   double Ix = 0,Iy = 0;
//   edge_t ei, ej;
//   for(oep1 = out_edges(u1,*Gp);oep1.first!=oep1.second;oep1.first++){
//     v1 = target(*oep1.first, (*Gp));
//     if(verbose)cout << " Finding 1st intersection:v1 is " << (*Gp)[v1].IsletCellNum << " at (" << (*Gp)[v1].x << "," << (*Gp)[v1].y << ")" <<  endl; cout.flush();
//     //Search for intersections within the component
//     for(int c = 0; c < (*Cp).size();c++){
//       vertex_t u2r = (*Cp)[c];
//       for(vp = vertices(*Gp);vp.first!=vp.second;vp.first++){
//   	if((*RG)[u2r].IsletCellNum == (*Gp)[*vp.first].IsletCellNum){
//   	  u2 = *vp.first;
//   	  break;
//   	}
//        }
//       if(u2!=u1 && u2!=v1){
//   	if(verbose)cout << "  u2 is " << (*Gp)[u2].IsletCellNum << " at (" << (*Gp)[u2].x << "," << (*Gp)[u2].y << ")" <<  endl; cout.flush();
//   	for(oep2 = out_edges(u2,*Gp);oep2.first!=oep2.second;oep2.first++){
//   	  v2 = target(*oep2.first, (*Gp));
//   	  if(verbose)cout << "   v2 is " << (*Gp)[v2].IsletCellNum << " at (" << (*Gp)[v2].x << "," << (*Gp)[v2].y << ")" <<  endl; cout.flush();
//   	  if(u1 !=v2 && v1!=v2 && (*Gp)[u2].CellNum < (*Gp)[v2].CellNum){
//   	    bool intersect = 0;
//   	    if((*Gp)[u1].IsletCellNum == 43)verbose=0;
//   	    else verbose=0;
//   	    intersect = IntersectionPointBetweenLines((*Gp)[u1].x, (*Gp)[u1].y,(*Gp)[v1].x, (*Gp)[v1].y,
//   						      (*Gp)[u2].x, (*Gp)[u2].y,(*Gp)[v2].x, (*Gp)[v2].y,
//   						      &Ix, &Iy,verbose);
//   	    if(verbose)cout << "intersect = " << intersect <<  endl;
//   	    if((*Gp)[u1].IsletCellNum == 43)verbose=0;;
//   	    if(intersect){
//   	      if((*Gp)[u1].IsletCellNum == 43)verbose=0;;
//   	      if(verbose)cout << "There is an intersection between (" << (*Gp)[u1].IsletCellNum << "," << (*Gp)[v1].IsletCellNum << ") and (" << (*Gp)[u2].IsletCellNum << "," << (*Gp)[v2].IsletCellNum << ") with IP = (" << Ix << "," << Iy << endl;
//   	      intersections = 1;
//   	      ei = *oep1.first; ej = *oep2.first;
//   	      while(oep1.first!= oep1.second)oep1.first++;
//   	      oep1.first--;
//   	      c = (*Cp).size();
//   	      if((*Gp)[u1].IsletCellNum == 43)verbose=0;;
//   	      break;
//   	    }
//   	  }
//   	}
//       }
//     }
//      if(!intersections){
//       //Check intersection points
//       vp = vertices(*Gp);
//       while(vp.first!=vp.second &&(*Gp)[*vp.first].CellNum < 100000 )vp.first++;
//       if(vp.first!=vp.second){
//     	for(;vp.first !=vp.second;vp.first++){
//     	  vertex_t u2 = *vp.first;
//     	  if(verbose)cout << "  Checking IPs:u2 is " << (*Gp)[u2].IsletCellNum << " at (" << (*Gp)[u2].x << "," << (*Gp)[u2].y << ")" <<  endl; cout.flush();
//     	  for(oep2 = out_edges(u2,*Gp);oep2.first!=oep2.second;oep2.first++){
//     	    v2 = target(*oep2.first, (*Gp));
//     	    if(verbose)cout << "   v2 is " << (*Gp)[v2].IsletCellNum << " at (" << (*Gp)[v2].x << "," << (*Gp)[v2].y << ")" <<  endl; cout.flush();
//     	    if(u1 !=v2 && v1!=v2 && (*Gp)[u2].CellNum < (*Gp)[v2].CellNum ){
//     	      bool intersect = 0;
//     	      verbose=0;
//     	      intersect = IntersectionPointBetweenLines((*Gp)[u1].x, (*Gp)[u1].y,(*Gp)[v1].x, (*Gp)[v1].y,
//     	    						(*Gp)[u2].x, (*Gp)[u2].y,(*Gp)[v2].x, (*Gp)[v2].y,
//     	    						&Ix, &Iy,verbose);
//     	      if(verbose)cout << "intersect = " << intersect <<  endl;
//     	      if((*Gp)[u1].IsletCellNum == 43)verbose=0;;
	     
//     	      if(intersect){
//     	      	if((*Gp)[u1].IsletCellNum == 43)verbose=0;;
//     	      	if(verbose)cout << "There is an intersection between (" << (*Gp)[u1].IsletCellNum << "," << (*Gp)[v1].IsletCellNum << ") and (" << (*Gp)[u2].IsletCellNum << "," << (*Gp)[v2].IsletCellNum << ") with IP = (" << Ix << "," << Iy << endl;
//     	      	intersections = 1;
//     	      	ei = *oep1.first; ej = *oep2.first;
//     	      	while(oep1.first!= oep1.second)oep1.first++;
//     	      	oep1.first--;
//     	      	while(vp.first!= vp.second)vp.first++;
//     	      	vp.first--;
//     	      	if((*Gp)[u1].IsletCellNum == 43)verbose=0;;
//     	      	break;
//     	      }
//     	    }
//     	  }
//     	}
//       }
//     }
//    }
//   if(verbose && !intersections)cout << "There are no intersections" << endl;
//   while(intersections){
//     if(verbose)cout << "Adding point at (" << Ix << "," << Iy << ")" << endl;cout.flush(); 
//     AddIPtoGraph(I, Gp, Ix, Iy, ei, ej);
//     intersections=0;
//     //Search to see if any other intersections
//     for(oep1 = out_edges(u1,*Gp);oep1.first!=oep1.second;oep1.first++){
//       v1 = target(*oep1.first, (*Gp));
//       if(verbose)cout << " Checking further:v1 is " << (*Gp)[v1].IsletCellNum << " at (" << (*Gp)[v1].x << "," << (*Gp)[v1].y << ")" <<  endl; cout.flush();
//       //Search for intersections within the component
//       for(int c = 0; c < (*Cp).size();c++){
//   	vertex_t u2r = (*Cp)[c];
//   	for(vp = vertices(*Gp);vp.first!=vp.second;vp.first++){
//   	  if((*RG)[u2r].IsletCellNum == (*Gp)[*vp.first].IsletCellNum){
//   	    u2 = *vp.first;
//   	    if(verbose)cout << "  u2 is " << (*Gp)[u2].IsletCellNum << " at (" << (*Gp)[u2].x << "," << (*Gp)[u2].y << ")" <<  endl; cout.flush();
//   	    break;
//   	  }
//   	}
//   	if(u2!=u1 && u2!=v1){ 
//   	  for(oep2 = out_edges(u2,*Gp);oep2.first!=oep2.second;oep2.first++){
//   	    if(verbose)cout << "   v2 is " << (*Gp)[v2].IsletCellNum << " at (" << (*Gp)[v2].x << "," << (*Gp)[v2].y << ")" <<  endl; cout.flush();
//   	    v2 = target(*oep2.first, (*Gp));
//   	    if(u1 !=v2 && v1!=v2){
//   	      bool intersect = 0;
//   	      verbose=0;
//   	      intersect = IntersectionPointBetweenLines((*Gp)[u1].x, (*Gp)[u1].y,(*Gp)[v1].x, (*Gp)[v1].y,
//   							(*Gp)[u2].x, (*Gp)[u2].y,(*Gp)[v2].x, (*Gp)[v2].y,
//   							&Ix, &Iy,verbose);
//   	      if(verbose)cout << "intersect = " << intersect <<  endl;
//   	      if((*Gp)[u1].IsletCellNum == 43)verbose=0;;
//   	      if(intersect){
//   		if((*Gp)[u1].IsletCellNum == 43)verbose=0;;
//   		if(verbose)cout << "There is an intersection between (" << (*Gp)[u1].IsletCellNum << "," << (*Gp)[v1].IsletCellNum << ") and (" << (*Gp)[u2].IsletCellNum << "," << (*Gp)[v2].IsletCellNum << ") with IP = (" << Ix << "," << Iy << endl;
//   		intersections = 1;
//   		ei = *oep1.first; ej = *oep2.first;
//   		while(oep1.first!= oep1.second)oep1.first++;
//   		oep1.first--;
//   		c = (*Cp).size();
//   		if((*Gp)[u1].IsletCellNum == 43)verbose=0;;
//   		break;
//   	      }
//   	    }
//   	  }
//   	}
//       }
//       if(!intersections){
//   	if(verbose)cout << "no intersections between regular cells, checking IPs" << endl; cout.flush();
//   	if(verbose)cout << "Checking further IPs:u1 is " << (*Gp)[u1].IsletCellNum << " at (" << (*Gp)[u1].x << "," << (*Gp)[u1].y << ")" <<  endl; cout.flush();
//   	if(verbose)cout << " v1 is " << (*Gp)[v1].IsletCellNum << " at (" << (*Gp)[v1].x << "," << (*Gp)[v1].y << ")" <<  endl; cout.flush();
//   	//Check intersection points
//   	vp = vertices(*Gp);
//   	while((*Gp)[*vp.first].CellNum < 100000 && vp.first!=vp.second)vp.first++;
//   	if(vp.first!=vp.second){
//   	  for(;vp.first !=vp.second;vp.first++){
//   	    vertex_t u2 = *vp.first;
//   	    if(u1 != u2 && v1!=u2){
//   	      if(verbose) cout << "  u2 is " << (*Gp)[u2].IsletCellNum << " at (" << (*Gp)[u2].x << "," << (*Gp)[u2].y << ")" <<  endl; cout.flush();
//   	      for(oep2 = out_edges(u2,*Gp);oep2.first!=oep2.second;oep2.first++){
//   		v2 = target(*oep2.first, (*Gp));
//   		if(u1 !=v2 && v1!=v2){
//   		  if(verbose) cout << "   v2 is " << (*Gp)[v2].IsletCellNum << " at (" << (*Gp)[v2].x << "," << (*Gp)[v2].y << ")" <<  endl; cout.flush();
//   		  bool intersect = 0;
//   		  verbose=0;
//   		  intersect = IntersectionPointBetweenLines((*Gp)[u1].x, (*Gp)[u1].y,(*Gp)[v1].x, (*Gp)[v1].y,
//   							    (*Gp)[u2].x, (*Gp)[u2].y,(*Gp)[v2].x, (*Gp)[v2].y,
//   							    &Ix, &Iy,verbose);
//   		  if(verbose)cout << "intersect = " << intersect <<  endl;
		 
//   		  if((*Gp)[u1].IsletCellNum == 43)verbose=0;;
//   		  if(intersect){
//   		    if((*Gp)[u1].IsletCellNum == 43)verbose=0;;
//   		    if(verbose)cout << "There is an intersection between (" << (*Gp)[u1].IsletCellNum << "," << (*Gp)[v1].IsletCellNum << ") and (" << (*Gp)[u2].IsletCellNum << "," << (*Gp)[v2].IsletCellNum << ") with IP = (" << Ix << "," << Iy << endl; cout.flush();
//   		    intersections = 1;
//   		    ei = *oep1.first; ej = *oep2.first;
//   		    while(oep1.first!= oep1.second)oep1.first++;
//   		    oep1.first--;
//   		    while(vp.first!= vp.second)vp.first++;
//   		    vp.first--;
//   		    if((*Gp)[u1].IsletCellNum == 43)verbose=0;
//   		    break;
//   		  }
//   		}
//   	      }
//   	    }
//   	  }
//   	}
//       }
//     }
//   }
//   if(verbose)cout << "Finished MakeGraphPlanar " << endl;
//   cout.flush();
//   return;
// }

// //*********************************************************************************************//
// //*********************************************************************************************//
// //*********************************************************************************************//

// void AddIPtoGraph(IsletGraph * I, Graph * Gp,double Ix, double Iy,edge_t ei,edge_t ej){
//   vertex_t u1,u2,v1,v2,u;
//   edge_t e;
//   bool verbose=0,verbosemem=0;
//   if(verbose)cout << "Fixing intersection (" << Ix << "," << Iy << ")" << endl; 
//   //Check if u is already in IPs
//   bool inIPs = 0;
//   for(int i = 0;i < (*I).IPs.size();i++){
//     if(Ix == (*I).IPs[i].x && Iy == (*I).IPs[i].y){ //intersection point is already an IP of the islet
//       if(verbose)cout << " Adding an old IP, num =  " << (*I).IPs[i].CellNum <<  endl;
//       inIPs = 1;
//       u =add_vertex(*Gp);
//       (*Gp)[u].x = Ix;
//       (*Gp)[u].y = Iy;
//       (*Gp)[u].CellNum = (*I).IPs[i].CellNum;
//       (*Gp)[u].IsletCellNum =(*I).IPs[i].IsletCellNum;
//       (*Gp)[u].CompNum =(*Gp)[source(ei,*Gp)].CompNum;
//       if((*Gp)[source(ei,(*Gp))].type == "b")(*Gp)[u].type = "b";
//       else (*Gp)[u].type = "ad";
//       // for(int p = 0;p<(*I).IPs[i].IPparents.size();p++){
//       // 	cout << "Adding IPparent " << (*I).IPs[i].IPparents[p] << endl;
//       // 	(*Gp)[u].IPparents.push_back((*I).IPs[i].IPparents[p]);
//       // }
//       // //u1 = (*Gp)[u].IPparents[0];v1 = (*Gp)[u].IPparents[1];
//       // //u2 = (*Gp)[u].IPparents[2];v2 = (*Gp)[u].IPparents[3];
//       if(verbose)cout << "finished adding" << endl;
//       break;
//     }
//   }
//   if(verbose)cout <<"outside loop" << endl;
//   if(!inIPs){
//     int num=100000;
//     //update num
//     sort((*I).IPs.begin(),(*I).IPs.end(),compareByCellNum);
//     int ips = 0;
//     // num++;
//     if((*I).IPs.size()>0){
//       while((*I).IPs[ips].CellNum < num)ips++;
//       for(int i = ips; i < (*I).IPs.size();i++){
// 	if((*I).IPs[i].CellNum == num)num++;
// 	else break;
//       }
//     }

//     std::pair<vertex_iter,vertex_iter> vp;
//     int max = 99999;
//     for(vp = vertices(*Gp);vp.first!=vp.second; vp.first++){
//       if((*Gp)[*vp.first].IsletCellNum > max){
// 	max = (*Gp)[*vp.first].IsletCellNum;
// 	if(verbose)cout << "The maximal Islet cell num is " << max << endl;
//       }
//     }
//     //Check IPs also
//     if(verbose)cout << "For IPs:" << endl;
//     for(int ip = 0; ip < (*I).IPs.size();ip++){
//       if(max+1 == (*I).IPs[ip].CellNum)max++;
//       if(verbose)cout << "The maximal Islet cell num is " << max << " and (*I).IPs[ip].CellNum is " << (*I).IPs[ip].CellNum << endl;
//     }
//     // vp = vertices(*Gp);
//     // while((*Gp)[*vp.first].CellNum < num && vp.first!=vp.second)vp.first++;
//     // if(verbose)cout <<"num = " << num << " and vp is at " << (*Gp)[*vp.first].IsletCellNum << endl;
//     // while( (*Gp)[*vp.first].CellNum == num ){
//     //   if(verbose)cout <<"Inside while: num = " << num << " and vp is at " << (*Gp)[*vp.first].IsletCellNum << endl;
//     //   vp.first++;
//     //   num++;
//     //   if(vp.first==vp.second)break;
//     // }
//     num = max+1;
//     if(verbose)cout <<"After while: num = " << num << " and vp is at " << (*Gp)[*vp.first].IsletCellNum << endl;
//     //if(num == 100013)PrintoutGnuplotWithIPs(Gp,"temp.IPs_100010");
//     //Insert point at intersection
//     u =add_vertex(*Gp);
//     (*Gp)[u].x = Ix;
//     (*Gp)[u].y = Iy;
//     (*Gp)[u].CellNum = num;
//     (*Gp)[u].IsletCellNum = num;
//     (*Gp)[u].vother = num;
//     (*Gp)[u].CompNum = (*Gp)[source(ei,(*Gp))].CompNum;
//     if((*Gp)[source(ei,(*Gp))].type == "b")(*Gp)[u].type = "b";
//     else (*Gp)[u].type = "ad";
//     if(verbose) cout <<  "Adding a new IP: num = " << num << " at point (" << (*Gp)[u].x << "," << (*Gp)[u].y << ")" <<  endl;
//   }
//   //Add IPparents
//   u1 = source(ei, (*Gp));
//   v1 = target(ei, (*Gp));
//   u2 = source(ej, (*Gp));
//   v2 = target(ej, (*Gp));
  
//   if(verbose)cout << " For cell " << (*Gp)[u].IsletCellNum << ": u1 =  "<< (*Gp)[u1].IsletCellNum  << " v1 =  "<< (*Gp)[v1].IsletCellNum  <<" u2 =  "<< (*Gp)[u2].IsletCellNum  << " v2 =  "<< (*Gp)[v2].IsletCellNum  << endl; 
//   if((*Gp)[u1].CellNum > 99999 && (*Gp)[v1].CellNum > 99999){
//     if(verbose)cout << "u1 and v1 are IPs" << endl;
//     //check parents for matching set
//     if(verbose){
//       cout << "\t Parents of u1 are " ;
//       for (int i = 0; i < (*Gp)[u1].IPparents.size(); i++)cout << (*Gp)[(*Gp)[u1].IPparents[i]].IsletCellNum << " ";
//       cout << endl;
//       cout << "\t Parents of v1 are " ;
//       for (int i = 0; i < (*Gp)[v1].IPparents.size(); i++)cout << (*Gp)[(*Gp)[v1].IPparents[i]].IsletCellNum << " ";
//       cout << endl;
//     }
//     if((*Gp)[u1].IPparents[0] == (*Gp)[v1].IPparents[0] &&
//        (*Gp)[u1].IPparents[1] == (*Gp)[v1].IPparents[1]){
//       (*Gp)[u].IPparents.push_back((*Gp)[u1].IPparents[0]);(*Gp)[u].IPparents.push_back((*Gp)[u1].IPparents[1]);
//     }
//     else if((*Gp)[u1].IPparents[0] == (*Gp)[v1].IPparents[1] &&
// 	    (*Gp)[u1].IPparents[1] == (*Gp)[v1].IPparents[0]){
//       (*Gp)[u].IPparents.push_back((*Gp)[u1].IPparents[0]);(*Gp)[u].IPparents.push_back((*Gp)[u1].IPparents[1]);
//     }
//     else if((*Gp)[u1].IPparents[0] == (*Gp)[v1].IPparents[2] &&
// 	    (*Gp)[u1].IPparents[1] == (*Gp)[v1].IPparents[3]){
//       (*Gp)[u].IPparents.push_back((*Gp)[u1].IPparents[0]);(*Gp)[u].IPparents.push_back((*Gp)[u1].IPparents[1]);
//     }
//     else if((*Gp)[u1].IPparents[0] == (*Gp)[v1].IPparents[3] &&
// 	    (*Gp)[u1].IPparents[1] == (*Gp)[v1].IPparents[2]){
//       (*Gp)[u].IPparents.push_back((*Gp)[u1].IPparents[0]);(*Gp)[u].IPparents.push_back((*Gp)[u1].IPparents[1]);
//     }
    
//     else if((*Gp)[u1].IPparents[2] == (*Gp)[v1].IPparents[0] &&
// 	    (*Gp)[u1].IPparents[3] == (*Gp)[v1].IPparents[1]){
//       (*Gp)[u].IPparents.push_back((*Gp)[u1].IPparents[2]);(*Gp)[u].IPparents.push_back((*Gp)[u1].IPparents[3]);
//     }
    
//     else if((*Gp)[u1].IPparents[2] == (*Gp)[v1].IPparents[1] &&
// 	    (*Gp)[u1].IPparents[3] == (*Gp)[v1].IPparents[0]){
//       (*Gp)[u].IPparents.push_back((*Gp)[u1].IPparents[2]);(*Gp)[u].IPparents.push_back((*Gp)[u1].IPparents[3]);
//     }
//     else if((*Gp)[u1].IPparents[2] ==(*Gp)[v1].IPparents[2] &&
// 	    (*Gp)[u1].IPparents[3] ==(*Gp)[v1].IPparents[3]){
//       (*Gp)[u].IPparents.push_back((*Gp)[u1].IPparents[2]);(*Gp)[u].IPparents.push_back((*Gp)[u1].IPparents[3]);
//     }
    
//     else if((*Gp)[u1].IPparents[2] ==(*Gp)[v1].IPparents[3] &&
// 	    (*Gp)[u1].IPparents[3] ==(*Gp)[v1].IPparents[2]){
//       (*Gp)[u].IPparents.push_back((*Gp)[u1].IPparents[2]);(*Gp)[u].IPparents.push_back((*Gp)[u1].IPparents[3]);
//     }
//   }
//   else{
//     if((*Gp)[u1].CellNum < 99999) (*Gp)[u].IPparents.push_back(source(ei, (*Gp)));//  not IP
//     else{
//       if(verbose){
// 	cout << "\t Parents of u1 are " ;
// 	for (int i = 0; i < (*Gp)[u1].IPparents.size(); i++)cout << (*Gp)[(*Gp)[u1].IPparents[i]].IsletCellNum << " ";
// 	cout << endl;
//       }
//       //Check the parents of u1 for v1
//       if((*Gp)[u1].IPparents[0] == v1){
// 	if(verbose) cout << "\t\tAdding " << (*Gp)[(*Gp)[u1].IPparents[1]].CellNum << endl;
// 	(*Gp)[u].IPparents.push_back((*Gp)[u1].IPparents[1]);
//       }
//       else if((*Gp)[u1].IPparents[1] == v1){
// 	if(verbose) cout << "\t\tAdding " << (*Gp)[(*Gp)[u1].IPparents[0]].CellNum << endl;
// 	(*Gp)[u].IPparents.push_back((*Gp)[u1].IPparents[0]);
//       }
//       else if((*Gp)[u1].IPparents[2] == v1){
// 	if(verbose) cout << "\t\tAdding " << (*Gp)[(*Gp)[u1].IPparents[3]].CellNum << endl;
// 	(*Gp)[u].IPparents.push_back((*Gp)[u1].IPparents[3]);
//       }
//       else if((*Gp)[u1].IPparents[3] == v1){
// 	if(verbose) cout << "\t\tAdding " << (*Gp)[(*Gp)[u1].IPparents[2]].CellNum << endl;
// 	(*Gp)[u].IPparents.push_back((*Gp)[u1].IPparents[2]);
//       }
//       else if(verbose) cout << " For cell " << (*Gp)[u].IsletCellNum << " not adding u1 parent for edge " <<  (*Gp)[u1].IsletCellNum << "," << (*Gp)[v1].IsletCellNum << endl;
//     }
//     if((*Gp)[v1].CellNum < 99999) (*Gp)[u].IPparents.push_back(target(ei, (*Gp)));//  not IP
//     else{
//       //Check the parents of v1 for u1
//       if((*Gp)[v1].IPparents[0] == u1)
// 	(*Gp)[u].IPparents.push_back((*Gp)[v1].IPparents[1]);
//       else if((*Gp)[v1].IPparents[1] == u1)
// 	(*Gp)[u].IPparents.push_back((*Gp)[v1].IPparents[0]);
//       else if((*Gp)[v1].IPparents[2] == u1)
// 	(*Gp)[u].IPparents.push_back((*Gp)[v1].IPparents[3]);
//       else if((*Gp)[v1].IPparents[3] == u1)
// 	(*Gp)[u].IPparents.push_back((*Gp)[v1].IPparents[2]);
//       else if(verbose) cout << " For cell " << (*Gp)[u].IsletCellNum << " not adding v1 parent for edge " <<  (*Gp)[u1].IsletCellNum << "," << (*Gp)[v1].IsletCellNum << endl;
//     }
//   }
//   if((*Gp)[u2].CellNum > 99999 && (*Gp)[v2].CellNum > 99999){
//     if(verbose) cout << "u2 and v2 are IPs" << endl;
//     //check parents for matching set
//     if(verbose){
//       cout << "\t Parents of u2 are " ;
//       for (int i = 0; i < (*Gp)[u2].IPparents.size(); i++)cout << (*Gp)[(*Gp)[u2].IPparents[i]].IsletCellNum << " ";
//       cout << endl;
//       cout << "\t Parents of v2 are " ;
//       for (int i = 0; i < (*Gp)[v2].IPparents.size(); i++)cout << (*Gp)[(*Gp)[v2].IPparents[i]].IsletCellNum << " ";
//       cout << endl;
//     }
//     if((*Gp)[u2].IPparents[0] ==(*Gp)[v2].IPparents[0] &&
//        (*Gp)[u2].IPparents[1] ==(*Gp)[v2].IPparents[1]){
//       (*Gp)[u].IPparents.push_back((*Gp)[u2].IPparents[0]);(*Gp)[u].IPparents.push_back((*Gp)[u2].IPparents[1]);
//     }
//     else if((*Gp)[u2].IPparents[0] ==(*Gp)[v2].IPparents[1] &&
// 	    (*Gp)[u2].IPparents[1] ==(*Gp)[v2].IPparents[0]){
//       (*Gp)[u].IPparents.push_back((*Gp)[u2].IPparents[0]);(*Gp)[u].IPparents.push_back((*Gp)[u2].IPparents[1]);
//     }
//     else if((*Gp)[u2].IPparents[0] ==(*Gp)[v2].IPparents[2] &&
// 	    (*Gp)[u2].IPparents[1] ==(*Gp)[v2].IPparents[3]){
//       (*Gp)[u].IPparents.push_back((*Gp)[u2].IPparents[0]);(*Gp)[u].IPparents.push_back((*Gp)[u2].IPparents[1]);
//     }
//     else if((*Gp)[u2].IPparents[0] ==(*Gp)[v2].IPparents[3] &&
// 	    (*Gp)[u2].IPparents[1] ==(*Gp)[v2].IPparents[2]){
//       (*Gp)[u].IPparents.push_back((*Gp)[u2].IPparents[0]);(*Gp)[u].IPparents.push_back((*Gp)[u2].IPparents[1]);
//     }
    
//     else if((*Gp)[u2].IPparents[2] ==(*Gp)[v2].IPparents[0] &&
// 	    (*Gp)[u2].IPparents[3] ==(*Gp)[v2].IPparents[1]){
//       (*Gp)[u].IPparents.push_back((*Gp)[u2].IPparents[2]);(*Gp)[u].IPparents.push_back((*Gp)[u2].IPparents[3]);
//     }
    
//     else if((*Gp)[u2].IPparents[2] ==(*Gp)[v2].IPparents[1] &&
// 	    (*Gp)[u2].IPparents[3] ==(*Gp)[v2].IPparents[0]){
//       (*Gp)[u].IPparents.push_back((*Gp)[u2].IPparents[2]);(*Gp)[u].IPparents.push_back((*Gp)[u2].IPparents[3]);
//     }
//     else if((*Gp)[u2].IPparents[2] ==(*Gp)[v2].IPparents[2] &&
// 	    (*Gp)[u2].IPparents[3] ==(*Gp)[v2].IPparents[3]){
//       (*Gp)[u].IPparents.push_back((*Gp)[u2].IPparents[2]);(*Gp)[u].IPparents.push_back((*Gp)[u2].IPparents[3]);
//     }
    
//     else if((*Gp)[u2].IPparents[2] ==(*Gp)[v2].IPparents[3] &&
// 	    (*Gp)[u2].IPparents[3] ==(*Gp)[v2].IPparents[2]){
//       (*Gp)[u].IPparents.push_back((*Gp)[u2].IPparents[2]);(*Gp)[u].IPparents.push_back((*Gp)[u2].IPparents[3]);
//     }
//   }
//   else{
//     if((*Gp)[u2].CellNum < 99999) (*Gp)[u].IPparents.push_back(source(ej, (*Gp)));//  not IP
//     else{
//       //Check the parents of u2 for v2
//       if((*Gp)[u2].IPparents[0] == v2)
// 	(*Gp)[u].IPparents.push_back((*Gp)[u2].IPparents[1]);
//       else if((*Gp)[u2].IPparents[1] == v2)
// 	(*Gp)[u].IPparents.push_back((*Gp)[u2].IPparents[0]);
//       else if((*Gp)[u2].IPparents[2] == v2)
// 	(*Gp)[u].IPparents.push_back((*Gp)[u2].IPparents[3]);
//       else if((*Gp)[u2].IPparents[3] == v2)
// 	(*Gp)[u].IPparents.push_back((*Gp)[u2].IPparents[2]);
//       else  if(verbose)cout << " For cell " << (*Gp)[u].IsletCellNum << " not adding u2 parent for edge " <<  (*Gp)[u2].IsletCellNum << "," << (*Gp)[v2].IsletCellNum << endl;
//     }
//     if((*Gp)[v2].CellNum < 99999) (*Gp)[u].IPparents.push_back(target(ej, (*Gp)));//  not IP
//     else{
//       if(verbose){
// 	cout << "\t Parents of v2 are " ;
// 	for (int i = 0; i < (*Gp)[v2].IPparents.size(); i++)cout << (*Gp)[(*Gp)[v2].IPparents[i]].CellNum << " ";
// 	cout << endl;
//       }
//       //Check the parents of v2 for u2
//       if((*Gp)[v2].IPparents[0] == u2){
// 	if(verbose) cout << "\t\tAdding " << (*Gp)[(*Gp)[v2].IPparents[1]].CellNum << endl;
// 	(*Gp)[u].IPparents.push_back((*Gp)[v2].IPparents[1]);
//       }
//       else if((*Gp)[v2].IPparents[1] == u2){
// 	if(verbose) cout << "\t\tAdding " << (*Gp)[(*Gp)[v2].IPparents[0]].CellNum << endl;
// 	(*Gp)[u].IPparents.push_back((*Gp)[v2].IPparents[0]);
//       }
//       else if((*Gp)[v2].IPparents[2] == u2){
// 	if(verbose) cout << "\t\tAdding " << (*Gp)[(*Gp)[v2].IPparents[3]].CellNum << endl;
// 	(*Gp)[u].IPparents.push_back((*Gp)[v2].IPparents[3]);
//       }
//       else if((*Gp)[v2].IPparents[3] == u2){
// 	if(verbose) cout << "\t\tAdding " << (*Gp)[(*Gp)[v2].IPparents[2]].CellNum << endl;
// 	(*Gp)[u].IPparents.push_back((*Gp)[v2].IPparents[2]);
//       }
//       else  if(verbose)cout << " For cell " << (*Gp)[u].CellNum << " not adding v2 parent for edge " <<  (*Gp)[u2].CellNum << "," << (*Gp)[v2].CellNum << endl;
//     }
//   }
  
//   if((*Gp)[u].IPparents.size() !=4){
//     cout << "*****Problem with IPparents for cell " <<  (*Gp)[u].IsletCellNum << "*********" << endl;
//     cout << " only " << (*Gp)[u].IPparents.size() <<" parents added " << endl;
//   }
//   if(verbose){
//     cout << "Added Parents " ;
//     for (int i = 0; i < (*Gp)[u].IPparents.size(); i++)cout << (*Gp)[(*Gp)[u].IPparents[i]].IsletCellNum << " ";
//     cout << endl << endl;
//   }
  
  
//   //Add edges between points and intersection
//   bool b = 0;
//   verbose=0;
//   if(verbosemem)cout << "Adding edges between (" << (*Gp)[u1].IsletCellNum << "," << (*Gp)[u].IsletCellNum << ") ("<< (*Gp)[u2].IsletCellNum << "," << (*Gp)[u].IsletCellNum << ") ("<< (*Gp)[v1].IsletCellNum << "," << (*Gp)[u].IsletCellNum << ") ("<< (*Gp)[v2].IsletCellNum << "," << (*Gp)[u].IsletCellNum << ") " << endl;cout.flush();
//   tie(e,b) = add_edge(u1, u, (*Gp));
//   double dist = Distancexy((*Gp)[u1].x, (*Gp)[u1].y, (*Gp)[u].x, (*Gp)[u].y);
//   InitializeEdge(Gp,e, dist,dist,2); 
//   tie(e,b) = add_edge(u2, u, (*Gp));
//   dist = Distancexy((*Gp)[u2].x, (*Gp)[u2].y, (*Gp)[u].x, (*Gp)[u].y);
//   InitializeEdge(Gp,e, dist,dist,2); 
//   tie(e,b) = add_edge(v1, u, (*Gp));
//   dist = Distancexy((*Gp)[v1].x, (*Gp)[v1].y, (*Gp)[u].x, (*Gp)[u].y);
//   InitializeEdge(Gp,e, dist,dist,2); 
//   tie(e,b) = add_edge(v2, u, (*Gp));
//   dist = Distancexy((*Gp)[v2].x, (*Gp)[v2].y, (*Gp)[u].x, (*Gp)[u].y);
//   InitializeEdge(Gp,e, dist,dist,2); 
  
//   //Remove edges ei and ej
//   ClearEdge(Gp,ei);
//   ClearEdge(Gp,ej);
//   remove_edge(ei,(*Gp));
//   remove_edge(ej,(*Gp));
  
//   return;
// }

// //*********************************************************************************************//
// //*********************************************************************************************//
// //*********************************************************************************************//

// void MakeGraphPlanar(IsletGraph * I,Graph * Gp){
//   std::pair<edge_iter, edge_iter> ep1,ep2;
//   vertex_t u1, v1, u2, v2, u;
//   edge_t e;
  
//     bool verbose=0,verbosemem=0;

//   int num = 100000;
//   sort((*I).IPs.begin(),(*I).IPs.end(),compareByCellNum);
//   for(int i = 0; i < (*I).IPs.size();i++){
//     if((*I).IPs[i].CellNum == num)num++;
//     else break;
//   }

//   if(verbose)cout << "Inside MakeGraphPlanar" << endl;

//   bool intersections = 0;
//   double Ix = 0,Iy = 0;
//   edge_t ei, ej;
//   for(ep1 = edges(*Gp);ep1.first!=ep1.second;ep1.first++){
//     u1 = source(*ep1.first, (*Gp));
//     v1 = target(*ep1.first, (*Gp));
//     ep2 = edges(*Gp);
//     while(ep2.first!=ep1.first)ep2.first++;
//     ep2.first++;
//     for(; ep2.first!=ep2.second;ep2.first++){
//       u2 = source(*ep2.first, (*Gp));
//       v2 = target(*ep2.first, (*Gp));
//       bool intersect = 0;
//       verbose=0;
//       intersect = IntersectionPointBetweenLines((*Gp)[u1].x, (*Gp)[u1].y,(*Gp)[v1].x, (*Gp)[v1].y,
// 						(*Gp)[u2].x, (*Gp)[u2].y,(*Gp)[v2].x, (*Gp)[v2].y,
// 						&Ix, &Iy,verbose);
      
//       if(verbose)cout << "intersect = " << intersect <<  endl;
//       verbose=0;
//       if(intersect){
// 	if(verbose)cout << "There is an intersection between (" << (*Gp)[u1].IsletCellNum << "," << (*Gp)[v1].IsletCellNum << ") and (" << (*Gp)[u2].IsletCellNum << "," << (*Gp)[v2].IsletCellNum << ") with IP = (" << Ix << "," << Iy << endl;
// 	intersections = 1;
// 	ei = *ep1.first; ej = *ep2.first;
// 	while(ep1.first!= ep1.second)ep1.first++;
// 	ep1.first--;
// 	break;
//       }
//     }
//   }
//   while(intersections){
//     verbose=0;
//     if(verbose)cout << "Fixing intersection (" << Ix << "," << Iy << ")" << endl; 
//     //Check if u is already in IPs
//     bool inIPs = 0;
//     for(int i = 0;i < (*I).IPs.size();i++){
//       if(Ix == (*I).IPs[i].x && Iy == (*I).IPs[i].y){ //intersection point is already an IP of the islet
// 	if(verbose)cout << " Adding an old IP, num =  " << (*I).IPs[i].CellNum <<  endl;
// 	inIPs = 1;
// 	u =add_vertex(*Gp);
// 	(*Gp)[u].x = Ix;
// 	(*Gp)[u].y = Iy;
// 	(*Gp)[u].CellNum = (*I).IPs[i].CellNum;
// 	(*Gp)[u].IsletCellNum =(*I).IPs[i].IsletCellNum;
// 	(*Gp)[u].CompNum =(*Gp)[(*I).IPs[i].IPparents[0]].CompNum;
// 	if((*Gp)[source(ei,(*Gp))].type == "b")(*Gp)[u].type = "b";
// 	else (*Gp)[u].type = "ad";
// 	if(verbose)cout << "finished adding" << endl;
// 	break;
//       }
//     }
//     if(verbose)cout <<"outside loop" << endl;
//     if(!inIPs){
//       //Insert point at intersection
//       u =add_vertex(*Gp);
//       (*Gp)[u].x = Ix;
//       (*Gp)[u].y = Iy;
//       (*Gp)[u].CellNum = num;
//       (*Gp)[u].IsletCellNum = num;
//       (*Gp)[u].vother = num;
//       (*Gp)[u].CompNum = (*Gp)[source(ei,(*Gp))].CompNum;
//       if((*Gp)[source(ei,(*Gp))].type == "b")(*Gp)[u].type = "b";
//       else (*Gp)[u].type = "ad";
//       if(verbose) cout <<  "Adding a new IP: num = " << num << endl;
//       //Add IPparents
//       u1 = source(ei, (*Gp));
//       v1 = target(ei, (*Gp));
//       u2 = source(ej, (*Gp));
//       v2 = target(ej, (*Gp));
     
//       //update num
//       sort((*I).IPs.begin(),(*I).IPs.end(),compareByCellNum);
//       int ips = 0;
//       num++;
//       if((*I).IPs.size()>0){
// 	while((*I).IPs[ips].CellNum < num)ips++;
// 	for(int i = ips; i < (*I).IPs.size();i++){
// 	  if((*I).IPs[i].CellNum == num)num++;
// 	  else break;
// 	}
//       }
//     }
//     if(verbose)cout << " For cell " << (*Gp)[u].IsletCellNum << ": u1 =  "<< (*Gp)[u1].IsletCellNum  << " v1 =  "<< (*Gp)[v1].IsletCellNum  <<" u2 =  "<< (*Gp)[u2].IsletCellNum  << " v2 =  "<< (*Gp)[v2].IsletCellNum  << endl; 
//     if((*Gp)[u1].CellNum > 99999 && (*Gp)[v1].CellNum > 99999){
//       if(verbose)cout << "u1 and v1 are IPs" << endl;
//       //check parents for matching set
//       if(verbose){
// 	cout << "\t Parents of u1 are " ;
// 	for (int i = 0; i < (*Gp)[u1].IPparents.size(); i++)cout << (*Gp)[(*Gp)[u1].IPparents[i]].IsletCellNum << " ";
// 	cout << endl;
// 	cout << "\t Parents of v1 are " ;
// 	for (int i = 0; i < (*Gp)[v1].IPparents.size(); i++)cout << (*Gp)[(*Gp)[v1].IPparents[i]].IsletCellNum << " ";
// 	cout << endl;
//       }
//       if((*Gp)[u1].IPparents[0] == (*Gp)[v1].IPparents[0] &&
// 	 (*Gp)[u1].IPparents[1] == (*Gp)[v1].IPparents[1]){
// 	(*Gp)[u].IPparents.push_back((*Gp)[u1].IPparents[0]);(*Gp)[u].IPparents.push_back((*Gp)[u1].IPparents[1]);
//       }
//       else if((*Gp)[u1].IPparents[0] == (*Gp)[v1].IPparents[1] &&
// 	      (*Gp)[u1].IPparents[1] == (*Gp)[v1].IPparents[0]){
// 	(*Gp)[u].IPparents.push_back((*Gp)[u1].IPparents[0]);(*Gp)[u].IPparents.push_back((*Gp)[u1].IPparents[1]);
//       }
//       else if((*Gp)[u1].IPparents[0] == (*Gp)[v1].IPparents[2] &&
// 	      (*Gp)[u1].IPparents[1] == (*Gp)[v1].IPparents[3]){
// 	(*Gp)[u].IPparents.push_back((*Gp)[u1].IPparents[0]);(*Gp)[u].IPparents.push_back((*Gp)[u1].IPparents[1]);
//       }
//       else if((*Gp)[u1].IPparents[0] == (*Gp)[v1].IPparents[3] &&
// 	      (*Gp)[u1].IPparents[1] == (*Gp)[v1].IPparents[2]){
// 	(*Gp)[u].IPparents.push_back((*Gp)[u1].IPparents[0]);(*Gp)[u].IPparents.push_back((*Gp)[u1].IPparents[1]);
//       }
      
//       else if((*Gp)[u1].IPparents[2] == (*Gp)[v1].IPparents[0] &&
// 	      (*Gp)[u1].IPparents[3] == (*Gp)[v1].IPparents[1]){
// 	(*Gp)[u].IPparents.push_back((*Gp)[u1].IPparents[2]);(*Gp)[u].IPparents.push_back((*Gp)[u1].IPparents[3]);
//       }
      
//       else if((*Gp)[u1].IPparents[2] == (*Gp)[v1].IPparents[1] &&
// 	      (*Gp)[u1].IPparents[3] == (*Gp)[v1].IPparents[0]){
// 	(*Gp)[u].IPparents.push_back((*Gp)[u1].IPparents[2]);(*Gp)[u].IPparents.push_back((*Gp)[u1].IPparents[3]);
//       }
//       else if((*Gp)[u1].IPparents[2] ==(*Gp)[v1].IPparents[2] &&
// 	      (*Gp)[u1].IPparents[3] ==(*Gp)[v1].IPparents[3]){
// 	(*Gp)[u].IPparents.push_back((*Gp)[u1].IPparents[2]);(*Gp)[u].IPparents.push_back((*Gp)[u1].IPparents[3]);
//       }
      
//       else if((*Gp)[u1].IPparents[2] ==(*Gp)[v1].IPparents[3] &&
// 	      (*Gp)[u1].IPparents[3] ==(*Gp)[v1].IPparents[2]){
// 	(*Gp)[u].IPparents.push_back((*Gp)[u1].IPparents[2]);(*Gp)[u].IPparents.push_back((*Gp)[u1].IPparents[3]);
//       }
//     }
//     else{
//       if((*Gp)[u1].CellNum < 99999) (*Gp)[u].IPparents.push_back(source(ei, (*Gp)));//  not IP
//       else{
// 	if(verbose){
// 	  cout << "\t Parents of u1 are " ;
// 	  for (int i = 0; i < (*Gp)[u1].IPparents.size(); i++)cout << (*Gp)[(*Gp)[u1].IPparents[i]].IsletCellNum << " ";
// 	  cout << endl;
// 	}
// 	//Check the parents of u1 for v1
// 	if((*Gp)[u1].IPparents[0] == v1){
// 	  if(verbose) cout << "\t\tAdding " << (*Gp)[(*Gp)[u1].IPparents[1]].CellNum << endl;
// 	  (*Gp)[u].IPparents.push_back((*Gp)[u1].IPparents[1]);
// 	}
// 	else if((*Gp)[u1].IPparents[1] == v1){
// 	  if(verbose) cout << "\t\tAdding " << (*Gp)[(*Gp)[u1].IPparents[0]].CellNum << endl;
// 	  (*Gp)[u].IPparents.push_back((*Gp)[u1].IPparents[0]);
// 	}
// 	else if((*Gp)[u1].IPparents[2] == v1){
// 	  if(verbose) cout << "\t\tAdding " << (*Gp)[(*Gp)[u1].IPparents[3]].CellNum << endl;
// 	  (*Gp)[u].IPparents.push_back((*Gp)[u1].IPparents[3]);
// 	}
// 	else if((*Gp)[u1].IPparents[3] == v1){
// 	  if(verbose) cout << "\t\tAdding " << (*Gp)[(*Gp)[u1].IPparents[2]].CellNum << endl;
// 	  (*Gp)[u].IPparents.push_back((*Gp)[u1].IPparents[2]);
// 	}
// 	else if(verbose) cout << " For cell " << (*Gp)[u].IsletCellNum << " not adding u1 parent for edge " <<  (*Gp)[u1].IsletCellNum << "," << (*Gp)[v1].IsletCellNum << endl;
//       }
//       if((*Gp)[v1].CellNum < 99999) (*Gp)[u].IPparents.push_back(target(ei, (*Gp)));//  not IP
//       else{
// 	//Check the parents of v1 for u1
// 	if((*Gp)[v1].IPparents[0] == u1)
// 	  (*Gp)[u].IPparents.push_back((*Gp)[v1].IPparents[1]);
// 	else if((*Gp)[v1].IPparents[1] == u1)
// 	  (*Gp)[u].IPparents.push_back((*Gp)[v1].IPparents[0]);
// 	else if((*Gp)[v1].IPparents[2] == u1)
// 	  (*Gp)[u].IPparents.push_back((*Gp)[v1].IPparents[3]);
// 	else if((*Gp)[v1].IPparents[3] == u1)
// 	  (*Gp)[u].IPparents.push_back((*Gp)[v1].IPparents[2]);
// 	else if(verbose) cout << " For cell " << (*Gp)[u].IsletCellNum << " not adding v1 parent for edge " <<  (*Gp)[u1].IsletCellNum << "," << (*Gp)[v1].IsletCellNum << endl;
//       }
//     }
//     if((*Gp)[u2].CellNum > 99999 && (*Gp)[v2].CellNum > 99999){
//       if(verbose) cout << "u2 and v2 are IPs" << endl;
//       //check parents for matching set
//       if(verbose){
// 	cout << "\t Parents of u2 are " ;
// 	for (int i = 0; i < (*Gp)[u2].IPparents.size(); i++)cout << (*Gp)[(*Gp)[u2].IPparents[i]].IsletCellNum << " ";
// 	cout << endl;
// 	cout << "\t Parents of v2 are " ;
// 	for (int i = 0; i < (*Gp)[v2].IPparents.size(); i++)cout << (*Gp)[(*Gp)[v2].IPparents[i]].IsletCellNum << " ";
// 	cout << endl;
//       }
//       if((*Gp)[u2].IPparents[0] ==(*Gp)[v2].IPparents[0] &&
// 	 (*Gp)[u2].IPparents[1] ==(*Gp)[v2].IPparents[1]){
// 	(*Gp)[u].IPparents.push_back((*Gp)[u2].IPparents[0]);(*Gp)[u].IPparents.push_back((*Gp)[u2].IPparents[1]);
//       }
//       else if((*Gp)[u2].IPparents[0] ==(*Gp)[v2].IPparents[1] &&
// 	      (*Gp)[u2].IPparents[1] ==(*Gp)[v2].IPparents[0]){
// 	(*Gp)[u].IPparents.push_back((*Gp)[u2].IPparents[0]);(*Gp)[u].IPparents.push_back((*Gp)[u2].IPparents[1]);
//       }
//       else if((*Gp)[u2].IPparents[0] ==(*Gp)[v2].IPparents[2] &&
// 	      (*Gp)[u2].IPparents[1] ==(*Gp)[v2].IPparents[3]){
// 	(*Gp)[u].IPparents.push_back((*Gp)[u2].IPparents[0]);(*Gp)[u].IPparents.push_back((*Gp)[u2].IPparents[1]);
//       }
//       else if((*Gp)[u2].IPparents[0] ==(*Gp)[v2].IPparents[3] &&
// 	      (*Gp)[u2].IPparents[1] ==(*Gp)[v2].IPparents[2]){
// 	(*Gp)[u].IPparents.push_back((*Gp)[u2].IPparents[0]);(*Gp)[u].IPparents.push_back((*Gp)[u2].IPparents[1]);
//       }
      
//       else if((*Gp)[u2].IPparents[2] ==(*Gp)[v2].IPparents[0] &&
// 	      (*Gp)[u2].IPparents[3] ==(*Gp)[v2].IPparents[1]){
// 	(*Gp)[u].IPparents.push_back((*Gp)[u2].IPparents[2]);(*Gp)[u].IPparents.push_back((*Gp)[u2].IPparents[3]);
//       }
      
//       else if((*Gp)[u2].IPparents[2] ==(*Gp)[v2].IPparents[1] &&
// 	      (*Gp)[u2].IPparents[3] ==(*Gp)[v2].IPparents[0]){
// 	(*Gp)[u].IPparents.push_back((*Gp)[u2].IPparents[2]);(*Gp)[u].IPparents.push_back((*Gp)[u2].IPparents[3]);
//       }
//       else if((*Gp)[u2].IPparents[2] ==(*Gp)[v2].IPparents[2] &&
// 	      (*Gp)[u2].IPparents[3] ==(*Gp)[v2].IPparents[3]){
// 	(*Gp)[u].IPparents.push_back((*Gp)[u2].IPparents[2]);(*Gp)[u].IPparents.push_back((*Gp)[u2].IPparents[3]);
//       }
      
//       else if((*Gp)[u2].IPparents[2] ==(*Gp)[v2].IPparents[3] &&
// 	      (*Gp)[u2].IPparents[3] ==(*Gp)[v2].IPparents[2]){
// 	(*Gp)[u].IPparents.push_back((*Gp)[u2].IPparents[2]);(*Gp)[u].IPparents.push_back((*Gp)[u2].IPparents[3]);
//       }
//     }
//     else{
//       if((*Gp)[u2].CellNum < 99999) (*Gp)[u].IPparents.push_back(source(ej, (*Gp)));//  not IP
//       else{
// 	//Check the parents of u2 for v2
// 	if((*Gp)[u2].IPparents[0] == v2)
// 	  (*Gp)[u].IPparents.push_back((*Gp)[u2].IPparents[1]);
// 	else if((*Gp)[u2].IPparents[1] == v2)
// 	  (*Gp)[u].IPparents.push_back((*Gp)[u2].IPparents[0]);
// 	else if((*Gp)[u2].IPparents[2] == v2)
// 	  (*Gp)[u].IPparents.push_back((*Gp)[u2].IPparents[3]);
// 	else if((*Gp)[u2].IPparents[3] == v2)
// 	  (*Gp)[u].IPparents.push_back((*Gp)[u2].IPparents[2]);
// 	else  if(verbose)cout << " For cell " << (*Gp)[u].IsletCellNum << " not adding u2 parent for edge " <<  (*Gp)[u2].IsletCellNum << "," << (*Gp)[v2].IsletCellNum << endl;
//       }
//       if((*Gp)[v2].CellNum < 99999) (*Gp)[u].IPparents.push_back(target(ej, (*Gp)));//  not IP
//       else{
// 	if(verbose){
// 	  cout << "\t Parents of v2 are " ;
// 	  for (int i = 0; i < (*Gp)[v2].IPparents.size(); i++)cout << (*Gp)[(*Gp)[v2].IPparents[i]].CellNum << " ";
// 	  cout << endl;
// 	}
// 	//Check the parents of v2 for u2
// 	if((*Gp)[v2].IPparents[0] == u2){
// 	  if(verbose) cout << "\t\tAdding " << (*Gp)[(*Gp)[v2].IPparents[1]].CellNum << endl;
// 	  (*Gp)[u].IPparents.push_back((*Gp)[v2].IPparents[1]);
// 	}
// 	else if((*Gp)[v2].IPparents[1] == u2){
// 	  if(verbose) cout << "\t\tAdding " << (*Gp)[(*Gp)[v2].IPparents[0]].CellNum << endl;
// 	  (*Gp)[u].IPparents.push_back((*Gp)[v2].IPparents[0]);
// 	}
// 	else if((*Gp)[v2].IPparents[2] == u2){
// 	  if(verbose) cout << "\t\tAdding " << (*Gp)[(*Gp)[v2].IPparents[3]].CellNum << endl;
// 	  (*Gp)[u].IPparents.push_back((*Gp)[v2].IPparents[3]);
// 	}
// 	else if((*Gp)[v2].IPparents[3] == u2){
// 	  if(verbose) cout << "\t\tAdding " << (*Gp)[(*Gp)[v2].IPparents[2]].CellNum << endl;
// 	  (*Gp)[u].IPparents.push_back((*Gp)[v2].IPparents[2]);
// 	}
// 	else  if(verbose)cout << " For cell " << (*Gp)[u].CellNum << " not adding v2 parent for edge " <<  (*Gp)[u2].CellNum << "," << (*Gp)[v2].CellNum << endl;
//       }
//     }
//     if((*Gp)[u].IPparents.size() !=4){
// 	cout << "*****Problem with IPparents for cell " <<  (*Gp)[u].IsletCellNum << "*********" << endl;
// 	cout << " only " << (*Gp)[u].IPparents.size() <<" parents added " << endl;
//     }
//     if(verbose){
//       cout << "Added Parents " ;
//       for (int i = 0; i < (*Gp)[u].IPparents.size(); i++)cout << (*Gp)[(*Gp)[u].IPparents[i]].IsletCellNum << " ";
//       cout << endl << endl;
//     }
    
   
//     //Add edges between points and intersection
//     bool b = 0;
//     verbose=0;
//     if(verbosemem)cout << "Adding edges between (" << (*Gp)[u1].IsletCellNum << "," << (*Gp)[u].IsletCellNum << ") ("<< (*Gp)[u2].IsletCellNum << "," << (*Gp)[u].IsletCellNum << ") ("<< (*Gp)[v1].IsletCellNum << "," << (*Gp)[u].IsletCellNum << ") ("<< (*Gp)[v2].IsletCellNum << "," << (*Gp)[u].IsletCellNum << ") " << endl;cout.flush();
//     tie(e,b) = add_edge(u1, u, (*Gp));
//     double dist = Distancexy((*Gp)[u1].x, (*Gp)[u1].y, (*Gp)[u].x, (*Gp)[u].y);
//     InitializeEdge(Gp,e, dist,dist,2); 
//     tie(e,b) = add_edge(u2, u, (*Gp));
//     dist = Distancexy((*Gp)[u2].x, (*Gp)[u2].y, (*Gp)[u].x, (*Gp)[u].y);
//     InitializeEdge(Gp,e, dist,dist,2); 
//     tie(e,b) = add_edge(v1, u, (*Gp));
//     dist = Distancexy((*Gp)[v1].x, (*Gp)[v1].y, (*Gp)[u].x, (*Gp)[u].y);
//     InitializeEdge(Gp,e, dist,dist,2); 
//     tie(e,b) = add_edge(v2, u, (*Gp));
//     dist = Distancexy((*Gp)[v2].x, (*Gp)[v2].y, (*Gp)[u].x, (*Gp)[u].y);
//     InitializeEdge(Gp,e, dist,dist,2); 
   
//     //Remove edges ei and ej
//     ClearEdge(Gp,ei);
//     ClearEdge(Gp,ej);
//     remove_edge(ei,(*Gp));
//     remove_edge(ej,(*Gp));
    
//     //Check if any other intersections
//     for(ep1 = edges(*Gp);ep1.first!=ep1.second;ep1.first++){
//       u1 = source(*ep1.first, (*Gp));
//       v1 = target(*ep1.first, (*Gp));
//       ep2 = edges(*Gp);
//       while(ep2.first!=ep1.first)ep2.first++;
//       ep2.first++;
//       for(; ep2.first!=ep2.second;ep2.first++){
// 	u2 = source(*ep2.first, (*Gp));
// 	v2 = target(*ep2.first, (*Gp));
// 	bool intersect = 0;
// 	intersections = 0;
// 	verbose=0;
// 	intersect = IntersectionPointBetweenLines((*Gp)[u1].x, (*Gp)[u1].y,(*Gp)[v1].x, (*Gp)[v1].y,
// 						  (*Gp)[u2].x, (*Gp)[u2].y,(*Gp)[v2].x, (*Gp)[v2].y,
// 						  &Ix, &Iy,verbose);
// 	verbose=0;
// 	if(intersect){
// 	  if(verbose)cout << "There is an intersection between (" << (*Gp)[u1].IsletCellNum << "," << (*Gp)[v1].IsletCellNum << ") and (" << (*Gp)[u2].IsletCellNum << "," << (*Gp)[v2].IsletCellNum << ") with IP = (" << Ix << "," << Iy << ")" << endl;
// 	  intersections = 1;
// 	  ei = *ep1.first; ej = *ep2.first;
// 	  while(ep1.first!= ep1.second)ep1.first++;
// 	  ep1.first--;
// 	  break;
// 	}
//       }
//     }
//   }
//   if(verbose)cout << "Finished MakeGraphPlanar " << endl;
//   cout.flush();
//   return;
// }

// //*********************************************************************************************//
// //*********************************************************************************************//
// //*********************************************************************************************//



// double MinBetween15and35(vector<vector<double> > * vec){
//    double minx = 0;
//    double yval = 100000;
//    int xpos = 0;
//    while((*vec)[xpos][0] < 15)xpos++;
//    while((*vec)[xpos][0] <= 35){
//      if((*vec)[xpos][1] < yval){
//        minx = (*vec)[xpos][0];
//        yval = (*vec)[xpos][1];
//      }
//      xpos++;
//    }
//    return minx;
// }
   

// //*********************************************************************************************//
// //*********************************************************************************************//
// //*********************************************************************************************//

// double MinBetweenPeaks2and3(vector<vector<double> > * vec, vector<vector<double> > * Peaks){
//   double minx = 0;
//   if((*Peaks).size() > 2){
//     int peak2xi = 0; 
//     while((*vec)[peak2xi][0] < (*Peaks)[1][0])peak2xi++;
//     int peak3xi = peak2xi;
//     while ((*vec)[peak3xi][0] < (*Peaks)[2][0])peak3xi++;
//     double miny = 1000000;
//     for(int i = peak2xi; i < peak3xi;i++){
//       if((*vec)[i][1] < miny){
// 	miny = (*vec)[i][1];
// 	minx = (*vec)[i][0];
//       }
//     }
//   }
//   else if((*Peaks).size() > 1){
//     int peak2xi = 0; 
//     while((*vec)[peak2xi][0] < (*Peaks)[1][0])peak2xi++;
//     int min3xi = peak2xi;
//     while ((*vec)[min3xi][1] > 1.0)min3xi++;
//     minx = min3xi;
//   }
//   else if((*Peaks).size() > 0){
//     int peak1xi = 0; 
//     while((*vec)[peak1xi][0] < (*Peaks)[0][0])peak1xi++;
//     int min3xi = peak1xi;
//     while ((*vec)[min3xi][1] > 1.0)min3xi++;
//     minx = min3xi;
//   }
//   else minx = 20;
//   return minx;
// }

// //*********************************************************************************************//
// //*********************************************************************************************//
// //*********************************************************************************************//


    


