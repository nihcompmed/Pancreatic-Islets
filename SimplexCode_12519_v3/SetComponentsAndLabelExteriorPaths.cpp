#include "./GraphNN_VersionForPaper.h"

#define PI 3.14159

using namespace boost;


void SetComponentsAndLabelExteriorPaths(IsletGraph* IGp,string type, string fprefix, string sthresh){

  std::pair<vertex_iter, vertex_iter> vp1,vp;
  std::pair<out_edge_iter, out_edge_iter> oep;
  vertex_t u;
  edge_t e;
  std::pair<edge_iter, edge_iter> ep;
  int graphcount=0;
  
  int vectorcount = 0;
  bool verbose=0,verbosemem=0;
  
  cout << "Setting " << type << "," << sthresh << " components and labels..." ;
  if(verbose)cout << endl;
  
  //*****************************************************************************//
  //                       Set Graph and component pointers
  //*****************************************************************************//
  
  Graph * Gp;
  vector<vector<vertex_t> > * Cp,* Pp,* Mp;
  vector<int> * Mtp, *Mtncomp;
  vector<vector<double> > * MtPerc;
  if(type == "all"){
    Gp = &(*IGp).allGraph;
    Cp = &(*IGp).allComponents;
    Pp = &(*IGp).allExteriorPaths;
  }
  if(type == "b"){
    Gp = &(*IGp).bGraph;
    Cp = &(*IGp).bComponents;
    Pp = &(*IGp).bExteriorPaths;
    Mp = &(*IGp).bMantle;
    Mtp = &(*IGp).bMantle_type;
    MtPerc = &(*IGp).bMantlePercent;
    //Mtncomp = &(*IGp).bMantle_ncomp;
  }
  if(type == "ad"){
    Gp = &(*IGp).adGraph;
    Cp = &(*IGp).adComponents;
    Pp = &(*IGp).adExteriorPaths;
    Mp = &(*IGp).adMantle;
    Mtp= &(*IGp).adMantle_type;
    MtPerc = &(*IGp).adMantlePercent;
    //Mtncomp = &(*IGp).adMantle_ncomp;
  }
  if(verbose) cout << "finished setting graphs and components \n" <<endl ;
  //Clear out previous components
  for(int c = (*Cp).size()-1; c > -1;c--)
    (*Cp)[c].clear();
  (*Cp).clear();
  //Clear out previous paths
  for(int p = (*Pp).size()-1; p > -1;p--)
    (*Pp)[p].clear();
  (*Pp).clear();

  //Compute components of graphs 
  std::vector<int> component(num_vertices(*Gp), 0);vectorcount++;if(verbosemem)cout << "Adding vector component: vectorcount = " << vectorcount << endl;cout.flush();
  int num = 0;
  num = connected_components_DS(Gp, &component);
  vector<int>NumCellsinComp(num,0), CellNuminComp(num,0); vectorcount+=2;if(verbosemem)cout << "Adding vectors NumCellsinComp and CellNuminComp: vectorcount = " << vectorcount << endl;cout.flush();
  if(verbose)cout << "There are " << num << " total " << type << " components" << endl;

  //Add component vectors to *Gp[Gno]
  vector<vertex_t> temp; vectorcount++;if(verbosemem)cout << "Adding vector temp: vectorcount = " << vectorcount << endl;cout.flush();
  vector<int> tempint; vectorcount++;if(verbosemem)cout << "Adding vector tempint: vectorcount = " << vectorcount << endl;cout.flush();
  vector<double> tempdouble; vectorcount++;if(verbosemem)cout << "Adding vector tempdouble: vectorcount = " << vectorcount << endl;cout.flush();
  for(int j = 0; j < num;j++){
    (*Cp).push_back(temp);
    (*Mp).push_back(temp);
    (*Mtp).push_back(0);
    //(*Mtncomp).push_back(0);
    (*MtPerc).push_back(tempdouble);
  }

  RemoveVector(&temp);vectorcount--;if(verbosemem)cout << "Removing vector temp: vectorcount = " << vectorcount << endl;cout.flush();
  RemoveVector(&tempint); vectorcount--;if(verbosemem)cout << "Removing vector tempint: vectorcount = " << vectorcount << endl;cout.flush();
  RemoveVector(&tempdouble); vectorcount--;if(verbosemem)cout << "Removing vector tempdouble: vectorcount = " << vectorcount << endl;cout.flush();
  int cn = 0;
  for(int j = 0; j < num_vertices(*Gp);j++)
    NumCellsinComp[component[j]]++;

  for(vp1=vertices(*Gp);vp1.first!=vp1.second;vp1.first++){
    (*Gp)[cn].CompNum = component[cn];
    //(*Gp)[cn].ncellcomp = NumCellsinComp[component[cn]];
    CellNuminComp[component[cn]]++;
    //Add vertices Components
    (*Cp)[component[cn]].push_back(*vp1.first);
    //if(verbose)cout << "Adding " << (*Gp)[*vp1.first].IsletCellNum << " to component " << component[cn] << endl;
    cn++;
  }
  
  // //Add in existing IPs from IsletGraph
  // if(verbose)cout << "There are " << (*IGp).IPs.size() << " listed IPs" << endl;
  // for(int i = 0; i < (*IGp).IPs.size();i++){
  //   if((*IGp).IPs[i].type == type){
  //     (*IGp).IPs[i].CompNum = (*Gp)[(*IGp).IPs[i].IPparents[0]].CompNum;
  //     if(verbose)cout << "Setting cell " << (*IGp).IPs[i].IsletCellNum << " comp num to " << (*IGp).IPs[i].CompNum  << endl;
  //     if(verbose){
  // 	cout << "parents are " ;
  // 	for(int p = 0;p<(*IGp).IPs[i].IPparents.size();p++)cout <<(*Gp)[(*IGp).IPs[i].IPparents[p]].IsletCellNum << " " ; 
  // 	cout << endl;
  //     }
  //   }
  // }
  if(verbose && (*Cp).size() > 0){
    cout << "The components are " << endl;
    for(int c = 0; c < (*Cp).size();c++){
      cout << c << ": " ;
      for(int cc = 0; cc < (*Cp)[c].size();cc++)cout << (*Gp)[(*Cp)[c][cc]].IsletCellNum << " " ;
      cout << endl;cout.flush();
    }
  }
  //PrintoutGnuplotWithIPs(Gp,"temp.comps."+type+"."+sthresh);   
  
  //Clear vectors
  RemoveVector(&component);vectorcount--;if(verbosemem)cout << "Removing vector component: vectorcount = " << vectorcount << endl;cout.flush();
  RemoveVector(&NumCellsinComp);vectorcount--;if(verbosemem)cout << "Removing vector NumCellsinComp: vectorcount = " << vectorcount << endl;cout.flush();
  RemoveVector(&CellNuminComp);vectorcount--;if(verbosemem)cout << "Removing vector CellNuminComp: vectorcount = " << vectorcount << endl;cout.flush();


  //*************************************************************************//
  //                          Finding Exterior Path                          //
  //*************************************************************************//
  
  ////Label cells as either isolated or exterior (default; interior cells will be updated later)
  ////If isolated or only 2 cells in component add exterior path
  for(int c = 0; c < (*Cp).size(); c++){
    verbose=0;

    //cout << "For component " << c << " which has " << (*Cp)[c].size() << " cells " << endl; cout.flush();
    
    //PrintoutGnuplotWithIPs(&PGp,"temp");
    //Label isolated cells
    if((*Cp)[c].size() == 1){
      //  (*Gp)[(*Cp)[c][0]].componentlabel = 0;
      vector<vertex_t> newPath; vectorcount++;if(verbosemem)cout << "Adding vector newPath: vectorcount = " << vectorcount << endl;cout.flush();
      newPath.push_back((*Cp)[c][0]);
      newPath.push_back((*Cp)[c][0]);
      (*Pp).push_back(newPath);
      RemoveVector(&newPath);vectorcount--; if(verbosemem)cout << "Removing vector newPath: vectorcount = " << vectorcount << endl;cout.flush();  
    }
    else if ((*Cp)[c].size() ==2){
      vector<vertex_t> newPath; vectorcount++;if(verbosemem)cout << "Adding vector newPath: vectorcount = " << vectorcount << endl;cout.flush();
      for(int i = 0; i < (*Cp)[c].size();i++){
  	      //(*Gp)[(*Cp)[c][i]].componentlabel = 2;
  	      newPath.push_back((*Cp)[c][i]);
      }
      newPath.push_back((*Cp)[c][0]);
      (*Pp).push_back(newPath);
      RemoveVector(&newPath);vectorcount--;if(verbosemem)cout << "Removing vector newPath: vectorcount = " << vectorcount << endl;cout.flush();
    }
    //If more than 2 cells in component, determine exterior path
    else{
    
      //Find path and check
      int pathnum = 1;
      vector<vertex_t> Pathold;vectorcount++;if(verbosemem)cout << "Adding vector Pathold: vectorcount = " << vectorcount << endl;cout.flush();
      //Find starting point
      double miny = 100000;
      vertex_t startv,nextv;

      for(int i = 0; i < (*Cp)[c].size();i++){
	        if((*Gp)[(*Cp)[c][i]].y <miny){
            miny = (*Gp)[i].y;
            startv = (*Cp)[c][i];
          }
      }

      if(verbose)cout << "The starting point is " << ((*Gp))[startv].IsletCellNum << endl;
      FindExteriorPath(IGp,&(*Gp),startv, &Pathold, &nextv, pathnum);
      //Add path to Exterior path list
      (*Pp).push_back(Pathold);

      verbose=0;
      if(verbose){
	        cout << "For component " << c << " containing " ;
	        for(int p = 0; p < (*Cp)[c].size();p++)
	          cout << (*Gp)[(*Cp)[c][p]].IsletCellNum << " ";
	        cout << endl;
              	cout << "Path contains ";
              	for(int p = 0; p < Pathold.size();p++){
	          if((*Gp)[Pathold[p]].CellNum > 99999) cout <<  ((*Gp))[Pathold[p]].IsletCellNum  << " ";
	          else cout << ((*Gp))[Pathold[p]].IsletCellNum  << " ";
	          //cout << endl;
	        }
	        cout << endl;
	        cout.flush();
      }

      verbose=0;

      //Remove path vectors
      RemoveVector(&Pathold);vectorcount--;if(verbosemem)cout << "Removing vector Pathold: vectorcount = " << vectorcount << endl;cout.flush();
    }
  }
   
  //RemoveGraph(&PGp);graphcount--;if(verbosemem)cout << "Removing graph PGp: graphcount = " << vectorcount << endl;cout.flush();
  if(verbose){
    cout << "The components are " << endl;
    for(int i = 0; i <(*Cp).size();i++){
	cout << "Component " << i << ":  " ;
	for(int j = 0; j < (*Cp)[i].size();j++){
	  cout << (*Gp)[(*Cp)[i][j]].IsletCellNum << " ";
	}
	cout << endl;
    }
    cout << endl;
    cout << "The paths are " << endl;
    for(int i = 0; i <(*Pp).size();i++){
      cout << "Path " << i << ":  " ;
      for(int j = 0; j < (*Pp)[i].size();j++){
	if((*Pp)[i][j] < 100000)
	  cout << (*Gp)[(*Pp)[i][j]].IsletCellNum << " ";
	else cout << (*Pp)[i][j] << " ";
      }
      cout << endl;
    }
    cout << endl;
    PrintoutGnuplotExteriorPaths(IGp,"temp", "./");
  }
  
  if(vectorcount != 0 || graphcount != 0) cout << "Problem with vectors and or graphs! vectorcount = " << vectorcount << " and graphcount = " << graphcount << " in SetComponents" << endl;
  if(verbosemem)cout << "*****Finished SetComponents" << endl;cout.flush();
  return;
}
  
//*********************************************************************************************//
//*********************************************************************************************//
//*********************************************************************************************//
void FindExteriorPath(IsletGraph* Ip,Graph* Gp ,vertex_t startv, vector<vertex_t>*Path , vertex_t *nextv , int pathnum){
  
  vertex_t u;
  bool verbose=0,verbosemem=0;
  int vectorcount = 0;
  
  (*Path).push_back(startv);

  double MidAngle = 0,MidAngleBegin = 0,MidAngleEnd = 99999;
  bool foundpath = 0;

  while(!foundpath){
    while((*Path).size() == 1 ||(*Path)[(*Path).size()-1]!=startv){// ||((*Path)[(*Path).size()-1]==startv && PathMidAngle[0] !=PathMidAngle[PathMidAngle.size()-1]) ){
        FindNextCellInPath(Ip,Gp, Path, &MidAngle,pathnum,nextv);

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
  }
  //PrintoutGnuplotWithIPs(Gp,"temp.ExtPath");
  if(vectorcount!=0)cout << "Problem in FindExteriorPath: vectorcount = " << vectorcount << endl;
  return;
}

//*********************************************************************************************//
//*********************************************************************************************//
//*********************************************************************************************//

void MakeGraphLocallyPlanar(IsletGraph * I,Graph * Gp, vertex_t u1){

  std::pair<out_edge_iter, out_edge_iter> oep1,oep2;
  std::pair<vertex_iter,vertex_iter>vp;
  // std::pair<edge_iter, edge_iter> ep2;
  vertex_t v1, u2, v2, u;
  edge_t e;
  //PrintoutGnuplot(Gp,"temp");
  bool verbose=0,verbosemem=0;

  //if((*Gp)[u1].IsletCellNum == 43)verbose=0;

  int num = 100000;
  // sort((*I).IPs.begin(),(*I).IPs.end(),compareByCellNum);
  // for(int i = 0; i < (*I).IPs.size();i++){
  //   if((*I).IPs[i].CellNum == num)num++;
  //   else break;
  // }
 
  vector<vertex_t> *Cp;
  //Graph *RG;
  
  if(verbose)  cout << "The compnum for u (Cell " << (*Gp)[u1].IsletCellNum << ") is " << (*Gp)[u1].CompNum << " with type " << (*Gp)[u1].type << endl;
  if((*Gp)[u1].type == "b"){
    Cp = &((*I).bComponents[(*Gp)[u1].CompNum]);
    //RG = &(*I).bGraph;
  }
  else{
    Cp = &((*I).adComponents[(*Gp)[u1].CompNum]);
    //RG = &(*I).adGraph;
  }
  if(verbose){
    cout << "Inside MakeGraphLocallyPlanar" << endl;
    cout << "The component contains" << endl;
    for(int i = 0;i < (*Cp).size();i++)
      cout << (*Gp)[(*Cp)[i]].IsletCellNum << " ";cout.flush();
    cout << endl;
  }
  //PrintoutGnuplotWithIPs(Gp,"temp");
  bool intersections;
  double Ix = 0,Iy = 0;
  edge_t ei, ej;

  do{

      intersections = 0;

      for(oep1 = out_edges(u1,*Gp);oep1.first!=oep1.second;oep1.first++){

          v1 = target(*oep1.first, (*Gp));

          //Search for intersections within the component
          for(int c = 0; c < (*Cp).size();c++){

              vertex_t u2r = (*Cp)[c];
              for(vp = vertices(*Gp);vp.first!=vp.second;vp.first++){
              	if((*Gp)[u2r].IsletCellNum == (*Gp)[*vp.first].IsletCellNum){
              	  u2 = *vp.first;
              	  break;
              	}
              }

              if(u2!=u1 && u2!=v1){

  	              for(oep2 = out_edges(u2,*Gp);oep2.first!=oep2.second;oep2.first++){

  	                  v2 = target(*oep2.first, (*Gp));

  	                  if(u1 !=v2 && v1!=v2 && (*Gp)[u2].CellNum < (*Gp)[v2].CellNum){

                            
  	                      intersections = IntersectionPointBetweenLines((*Gp)[u1].x, (*Gp)[u1].y,(*Gp)[v1].x, (*Gp)[v1].y,
  	      					                  (*Gp)[u2].x, (*Gp)[u2].y,(*Gp)[v2].x, (*Gp)[v2].y,
  	      					                  &Ix, &Iy,verbose);

                          if (intersections){
  	                          ei = *oep1.first; ej = *oep2.first;

                              //cout << "Adding point at (" << Ix << "," << Iy << ")" << endl;cout.flush(); 
                              AddIPtoGraph(I, Gp, Ix, Iy, ei, ej);

                              break;
                          }

                      }

                  }

              }

              if (intersections) break;

          }

          if (intersections) break;
          else{//Search for intersection with IP points

  	          vp = vertices(*Gp);
  	          while((*Gp)[*vp.first].CellNum < 100000 && vp.first!=vp.second)vp.first++;

  	          if(vp.first!=vp.second){
  	              for(;vp.first !=vp.second;vp.first++){
  	                  vertex_t u2 = *vp.first;
  	                  if(u1 != u2 && v1!=u2){
  	                    for(oep2 = out_edges(u2,*Gp);oep2.first!=oep2.second;oep2.first++){
  	      	                v2 = target(*oep2.first, (*Gp));
  	      	                if(u1 !=v2 && v1!=v2){

  	      	                  intersections = IntersectionPointBetweenLines((*Gp)[u1].x, (*Gp)[u1].y,(*Gp)[v1].x, (*Gp)[v1].y,
  	      	                					    (*Gp)[u2].x, (*Gp)[u2].y,(*Gp)[v2].x, (*Gp)[v2].y,
  	      	                					    &Ix, &Iy,verbose);
		       
  	      	                  if(intersections){

  	      	                    ei = *oep1.first; ej = *oep2.first;

                                //cout << "Adding point at (" << Ix << "," << Iy << ")" << endl;cout.flush(); 
                                AddIPtoGraph(I, Gp, Ix, Iy, ei, ej);
  	      	                    break;

  	      	                  }

  	      	                }

  	                    }

  	                  }

                      if (intersections) break;

                  }

              }

          }

      }


  }while(intersections);


  if(verbose)cout << "Finished MakeGraphPlanar " << endl;
  cout.flush();
  return;


}

//*********************************************************************************************//
//*********************************************************************************************//
//*********************************************************************************************//

void FindNextCellInPath(IsletGraph* Ip,Graph*Gp, vector<vertex_t>*Path, double* MidAngle, int pathnum,vertex_t* nextv){

  bool verbose=0,verbosemem=0;
  int vectorcount = 0;
  std::pair<out_edge_iter, out_edge_iter> oep;

  vertex_t u = (*Path)[(*Path).size()-1];

  if(verbose)cout << " For vertex " << (*Gp)[u].IsletCellNum << " with degree " << out_degree(u,*Gp) << endl;cout.flush();

  MakeGraphLocallyPlanar(Ip,Gp,u);

  if(verbose)cout << "Finished making " << (*Gp)[u].IsletCellNum << " locally planar" <<endl;cout.flush();

  //find angle of each edge of u and sort
  vector <CellAngle> CellAngles;vectorcount++;if(verbosemem)cout << "Adding CellAngles: vectorcount = " << vectorcount << endl; cout.flush();
  for(oep = out_edges(u,(*Gp));oep.first!=oep.second; oep.first++){
      CellAngle temp;
      temp.Cell = target(*oep.first,(*Gp));
      temp.Angle = GetAngleOfLine((*Gp)[u].x,(*Gp)[u].y, (*Gp)[target(*oep.first,(*Gp))].x,(*Gp)[target(*oep.first,(*Gp))].y);
      if (verbose) cout << "Edge with " << (*Gp)[temp.Cell].IsletCellNum << " angle is " << temp.Angle << endl;
      CellAngles.push_back(temp);
  }

  sort(CellAngles.begin(),CellAngles.end(),compareByAngle);

  if(verbose){
    cout << "The sorted angles are " << endl;
    for(int ca = 0; ca < CellAngles.size();ca++)
      cout << (*Gp)[CellAngles[ca].Cell].IsletCellNum << "\t" << CellAngles[ca].Angle << endl;
  }
  
  //find previous leg of path in list and add next angle to path (if pathsize = 1 use 3pi/2)
  int CApos = 0;
  if(verbose)cout << "The path has " << (*Path).size() << " vertices and pathnum = " << pathnum << endl;

  if((*Path).size() == 1){
      if(pathnum!= 1){
        for(int i = 0; i < CellAngles.size();i++){
	          if(verbose)
            {
                cout << "CellAngles[" << i << "].Cell is "\
                     <<(*Gp)[CellAngles[i].Cell].IsletCellNum\
                     <<" and nextv  = " <<  (*Gp)[*nextv].IsletCellNum << endl;
            }
	          if((*Gp)[CellAngles[i].Cell].IsletCellNum == (*Gp)[*nextv].IsletCellNum){
	            if(verbose)cout << "They are equal" << endl;
	            CApos = (i+1)%CellAngles.size();
	            break;
	          }
        }
      }
      else{
        if(verbose)cout << "pathnum = 1 and the last Cell angle is " << CellAngles[CellAngles.size()-1].Angle << endl; 
        if(CellAngles[CellAngles.size()-1].Angle > 3.0*PI/2.0){
	          while(CellAngles[CApos].Angle <3.0*PI/2.0)CApos++;
        }
      }
      (*Path).push_back(CellAngles[CApos].Cell);
      (*nextv) = CellAngles[CApos].Cell;
      if(verbose)cout << "Adding cell " << (*Gp)[(*Path)[(*Path).size()-1]].IsletCellNum << " to path " << endl;cout.flush();
  }
  else{

      int newpos = 0;
      if(CellAngles[CellAngles.size()-1].Cell != (*Path)[(*Path).size()-2]){
        for(int ca = 0; ca < CellAngles.size();ca++){
	          if((*Path)[(*Path).size()-2] == CellAngles[ca].Cell){
	            //   cout << "Previous path is from " << CellAngles[ca].Cell << " and next angle is at " << CellAngles[ca+1].Cell << endl;
	            newpos = ca+1;
	            
	            break;
	          }
        }
      }
      if(CellAngles.size() == 1){
          *MidAngle = mod(CellAngles[0].Angle+PI,2*PI);   
      }
      else if(newpos==0)
      {
          *MidAngle =\
              mod((CellAngles[0].Angle+(2*PI - CellAngles[CellAngles.size()-1].Angle))/2.0 +CellAngles[CellAngles.size()-1].Angle,2*PI);
      }
      else{
          *MidAngle = (CellAngles[newpos].Angle+CellAngles[newpos-1].Angle)/2.0;
      }
      if(verbose)cout << "MidAngle = " << *MidAngle <<  endl;
      (*Path).push_back(CellAngles[newpos].Cell);

  }

  RemoveVector(&CellAngles);vectorcount--;if(verbosemem)cout << "Removing CellAngles: vectorcount = " << vectorcount << endl; cout.flush();

  return;
}


  //*********************************************************************************************//
//*********************************************************************************************//
//*********************************************************************************************//

void AddIPtoGraph(IsletGraph * I, Graph * Gp,double Ix, double Iy,edge_t ei,edge_t ej){

  //cout << "Adding point at (" << Ix << "," << Iy << ")" << endl;cout.flush(); 

  vertex_t u1,u2,v1,v2,u;
  edge_t e;
  bool verbose=0,verbosemem=0;
  if(verbose)cout << "Fixing intersection (" << Ix << "," << Iy << ")" << endl; 

  //Find Cell num
  int maxcellnum= 0,maxIsletCellNum = 0;
  for(int i = 0; i < num_vertices(*Gp);i++){
    if((*Gp)[i].CellNum > maxcellnum)maxcellnum = (*Gp)[i].CellNum;
    if((*Gp)[i].IsletCellNum > maxIsletCellNum)maxIsletCellNum = (*Gp)[i].IsletCellNum;
  }

  if(maxIsletCellNum < 100000)maxIsletCellNum = 99999;

  //Add vertex
  u =add_vertex(*Gp);
  (*Gp)[u].x = Ix;
  (*Gp)[u].y = Iy;
  (*Gp)[u].CellNum = maxcellnum+1;
  (*Gp)[u].IsletCellNum =maxIsletCellNum+1;
  (*Gp)[u].CompNum =(*Gp)[source(ei,*Gp)].CompNum;
  if((*Gp)[source(ei,(*Gp))].type == "b"){
    (*Gp)[u].type = "b";
    (*I).bComponents[(*Gp)[source(ei,*Gp)].CompNum].push_back(u);
  }
  else{
    (*Gp)[u].type = "ad";
    (*I).adComponents[(*Gp)[source(ei,*Gp)].CompNum].push_back(u);
  }

  
  //Add edges between points and intersection
  bool b = 0;
  verbose=0;
  u1 = source(ei,*Gp);u2 = target(ei,*Gp);
  v1 = source(ej,*Gp);v2 = target(ej,*Gp);
  if(verbosemem)cout << "Adding edges between (" << (*Gp)[u1].IsletCellNum << "," << (*Gp)[u].IsletCellNum << ") ("<< (*Gp)[u2].IsletCellNum << "," << (*Gp)[u].IsletCellNum << ") ("<< (*Gp)[v1].IsletCellNum << "," << (*Gp)[u].IsletCellNum << ") ("<< (*Gp)[v2].IsletCellNum << "," << (*Gp)[u].IsletCellNum << ") " << endl;cout.flush();
  tie(e,b) = add_edge(u1, u, (*Gp));
  double dist = Distancexy((*Gp)[u1].x, (*Gp)[u1].y, (*Gp)[u].x, (*Gp)[u].y);
  InitializeEdge(Gp,e, dist); 
  tie(e,b) = add_edge(u2, u, (*Gp));
  dist = Distancexy((*Gp)[u2].x, (*Gp)[u2].y, (*Gp)[u].x, (*Gp)[u].y);
  InitializeEdge(Gp,e, dist); 
  tie(e,b) = add_edge(v1, u, (*Gp));
  dist = Distancexy((*Gp)[v1].x, (*Gp)[v1].y, (*Gp)[u].x, (*Gp)[u].y);
  InitializeEdge(Gp,e, dist); 
  tie(e,b) = add_edge(v2, u, (*Gp));
  dist = Distancexy((*Gp)[v2].x, (*Gp)[v2].y, (*Gp)[u].x, (*Gp)[u].y);
  InitializeEdge(Gp,e, dist); 
  
  //Remove edges ei and ej
  ClearEdge(Gp,ei);
  ClearEdge(Gp,ej);
  remove_edge(ei,(*Gp));
  remove_edge(ej,(*Gp));
  
  return;
}
 //*********************************************************************************************//
//*********************************************************************************************//
//*********************************************************************************************//

int connected_components_DS(Graph *Gp, vector<int> *components){
  int vectorcount = 0;
  bool verbose=0;
  //if(verbose)PrintoutGnuplot(Gp,"connectedcomponentsDS");

  std::pair<vertex_iter, vertex_iter> vp;
  vector<vertexListWeight> Assigned;vectorcount++; //contains vertex_t vertex, int Parent (which will be component number); and Weight (which will be distance);
  int Lastcompnum = -1;
  for(vp = vertices(*Gp);vp.first!=vp.second;vp.first++){
    vertex_t u = *vp.first;
    if(verbose)cout << "Finding component for " << (*Gp)[u].IsletCellNum <<endl;
    //Calculate distance between u and Assigned vertices
    for(int i = 0; i < Assigned.size();i++)
      Assigned[i].Weight = Distancexy((*Gp)[u].x, (*Gp)[u].y, (*Gp)[Assigned[i].vertex].x,(*Gp)[Assigned[i].vertex].y);
    //Sort Assigned by distance
    sort(Assigned.begin(),Assigned.end(),compareByWeight);
    if(verbose){
      cout << "The sorted Assigned vertices are:" << endl;
      for(int i = 0;i < Assigned.size();i++)
    	cout << (*Gp)[Assigned[i].vertex].IsletCellNum << "\t" << Assigned[i].Parent << "\t" << Assigned[i].Weight << endl;
    }
    //see if path exists between u and Assigned
    bool NoPathExists = 1;
    for(int i = 0; i < Assigned.size();i++){
      vector<vertex_t>path;vectorcount++;
      NoPathExists = PathBetweenTwoVertices(Gp, u, Assigned[i].vertex, &path);
      if(!NoPathExists){//Assign u to i's component
	if(verbose)cout << "A path exists between " << (*Gp)[u].IsletCellNum << " and " << (*Gp)[Assigned[i].vertex].IsletCellNum << endl;
	vertexListWeight vlw;
	vlw.vertex = u;
	vlw.Parent = Assigned[i].Parent;
	Assigned.push_back(vlw);
	RemoveVector(&path);vectorcount--;
	break;
      }
      RemoveVector(&path);vectorcount--;
    }
    if(NoPathExists){
      vertexListWeight vlw;
      vlw.vertex = u;
      Lastcompnum++;
      vlw.Parent = Lastcompnum;
      Assigned.push_back(vlw);
      if(verbose)cout << "Assigning " << (*Gp)[u].IsletCellNum << " to component " << Lastcompnum  << endl;
    }
  }
  for(int i = 0; i < Assigned.size();i++)
    (*components)[(*Gp)[Assigned[i].vertex].CellNum] = Assigned[i].Parent;
  if(verbose){
    sort(Assigned.begin(),Assigned.end(),compareByParent);
    cout << "The components are "<< endl;
    for(int i = 0; i < Assigned.size();i++)
      cout << (*Gp)[Assigned[i].vertex].IsletCellNum << "\t" << Assigned[i].Parent << endl;
  }
  RemoveVector(&Assigned);vectorcount--;
  cout << "Finished finding components" <<endl;
  if(vectorcount>0)cout << "Problem with connected_components_DS! vectorcount = " << vectorcount << endl; 
  return Lastcompnum+1;
}

//*********************************************************************************************//
//*********************************************************************************************//
//*********************************************************************************************//

bool PathBetweenTwoVertices(Graph * Gp, vertex_t u, vertex_t v, vector<vertex_t> * path){
  //Determines if path between two vertices exists
  //quicker than Shortest weighted path
  //returns a path between the two vertices in *path

  bool verbose=0,verbosemem=0;
  int vectorcount = 0;
if(verbosemem)cout << "*****Inside PathBetweenTwoVertices" << endl;
  vector<vertexListWeight> que1, que2;vectorcount+=2;if(verbosemem)cout << "Adding que1,que2: vectorcount = " << vectorcount << endl;cout.flush();

  
  if(verbose)cout << "Inside PathBetweenTwoVertices...u = " << (*Gp)[u].IsletCellNum << " and " << (*Gp)[v].IsletCellNum <<  endl;cout.flush();
  //string su = static_cast<ostringstream*>( &(ostringstream() <<(*Gp)[u].IsletCellNum ) )->str();
  //string sv = static_cast<ostringstream*>( &(ostringstream() <<(*Gp)[v].IsletCellNum ) )->str();
 
  //if(verbose)PrintoutGnuplot(Gp, "temp.Path"+su+"_"+sv);
  //Setup que1
  vertexListWeight temp;
  temp.vertex = u; temp.Parent = -1;
  que1.push_back(temp);
  std::pair<out_edge_iter, out_edge_iter> oep;
  bool foundPathIn1 = 0, foundPathIn2 = 0;
  for(oep=out_edges(u,*Gp); oep.first!=oep.second; oep.first++){
    vertex_t u1 = target(*oep.first,*Gp);
    temp.vertex = u1; temp.Parent = 0;
    que1.push_back(temp);
    if(u1 == v){
      if(verbose)cout << "Found in path 1" << endl;
      foundPathIn1 = 1;
      break;
    }
  }
  if(!foundPathIn1){
    temp.vertex = v; temp.Parent = -1;
    que2.push_back(temp);
    std::pair<out_edge_iter, out_edge_iter> oep;
    for(oep=out_edges(v,*Gp); oep.first!=oep.second; oep.first++){
      vertex_t v1 = target(*oep.first,*Gp);
      temp.vertex = v1; temp.Parent = 0;
      que2.push_back(temp);
      if(v1 == u){
	foundPathIn2 = 1;
	if(verbose)cout << "Found in path 2" << endl;
	break;
      }
    }
  }
  int i1 = 1, i2 = 1;
  while(!foundPathIn1 && !foundPathIn2 && i1 < que1.size() && i2 < que2.size()){
    if(verbose)cout << "*** i1 = " << i1 << " with vertex " << (*Gp)[que1[i1].vertex].IsletCellNum << " and i2 = " << i2 << " with vertex " << (*Gp)[que2[i2].vertex].IsletCellNum <<" ***" << endl;
    for(oep=out_edges(que1[i1].vertex,*Gp); oep.first!=oep.second; oep.first++){
      vertex_t u1 = target(*oep.first,*Gp);
      //Check if u1 is in que1
      bool inque1 = 0;
      for(int i = 0; i < que1.size();i++){
	if(que1[i].vertex == u1){
	  inque1 = 1;
	  break;
	}
      }
      if(!inque1){
	temp.vertex = u1; temp.Parent = i1;
	que1.push_back(temp);
      }
      if(u1 == v){
	foundPathIn1 = 1;
	if(verbose)cout << "Found in path 1" << endl;
	break;
      }
    }
    if(!foundPathIn1){
      temp.vertex = v; temp.Parent = 0;
      que2.push_back(temp);
      std::pair<out_edge_iter, out_edge_iter> oep;
      for(oep=out_edges(que2[i2].vertex,*Gp); oep.first!=oep.second; oep.first++){
	vertex_t v1 = target(*oep.first,*Gp);
	//Check if v1 is in que2
	bool inque2 = 0;
	for(int i = 0; i < que2.size();i++){
	  if(que2[i].vertex == v1){
	  inque2 = 1;
	  break;
	  }
	}
	if(!inque2){
	  temp.vertex = v1; temp.Parent = i2;
	  que2.push_back(temp);
	}
	if(v1 == u){
	  foundPathIn2 = 1;
	  if(verbose)cout << "Found in path 2" << endl;
	  break;
	}
      }
    }
    if(!foundPathIn1 && !foundPathIn2){
      i1++;i2++;
    }
  }
 
  if(que1.size() == i1 || que2.size() == i2){
    if(verbose && que1.size() == i1)cout << "que1.size = i1: no path" << endl;
    else if(verbose && que2.size() == i2)cout << "que2.size = i2: no path" << endl;
    RemoveVector(&que1);vectorcount--;if(verbosemem)cout << "Removing que1: vectorcount = " << vectorcount << endl;cout.flush();
    RemoveVector(&que2);vectorcount--;if(verbosemem)cout << "Removingque2: vectorcount = " << vectorcount << endl;cout.flush();
    if(vectorcount!=0)cout << "Problem with vectorcount! vectorcount = " << vectorcount << endl;
    return 1;
  }
  else if(foundPathIn1){
    if(verbose){
      cout << "que1 consists of " << endl;
      for(int i = 0; i < que1.size(); i++)
	cout << (*Gp)[que1[i].vertex].IsletCellNum << "\t" << que1[i].Parent << endl;
      cout.flush();
    }
    vector<vertex_t> temp2;vectorcount++;if(verbosemem)cout << "Adding temp2: vectorcount = " << vectorcount << endl;cout.flush();
    int i = que1.size()-1;
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
    int i = que2.size()-1;
    while(i !=-1){
      (*path).push_back(que2[i].vertex);
      i = que2[i].Parent;
    }
  }
  if(verbose){
    cout << "The path from vertex " << (*Gp)[u].IsletCellNum << " to " << (*Gp)[v].IsletCellNum << " is " << endl;cout.flush();
    for(int i = 0; i < (*path).size();i++)
      cout << (*Gp)[(*path)[i]].IsletCellNum << " ";cout.flush();
    cout << endl;cout.flush();
  }
  RemoveVector(&que1);vectorcount--;if(verbosemem)cout << "Removing que1: vectorcount = " << vectorcount << endl;cout.flush();
  RemoveVector(&que2);vectorcount--;if(verbosemem)cout << "Removing que2: vectorcount = " << vectorcount << endl;cout.flush();
  if(vectorcount!=0)cout << "Problem in PathBetweenTwoVertices with vectorcount! vectorcount = " << vectorcount << endl;cout.flush();

  return 0;
}

