#include "./GraphNN_VersionForPaper.h"

#define PI 3.14159

using namespace boost;

int insert_in_set(int v1, int v2){
    return ((v1 * (v1-1))/2) + v2;
}

// Default sthresh = 'gr' and type = 'all'
void FindAD_BetaSets(IsletGraph* I,string fprefix, string sthresh, string type){
  
  std::pair<edge_iter, edge_iter> ep;
  edge_t e;
  vertex_t u,v;
  bool verbose=0;

  //cout << "Finding AD and Beta Sets and Loops..." << endl;
  if(verbose)cout << "There are " << (*I).nbeta << " betas, " << (*I).nalpha << " alphas and " << (*I).ndelta << " deltas " << endl;

  // Do only if # beta cells > 5 and ad cells > 5
  if((*I).nbeta > 5 && (*I).nalpha+(*I).ndelta > 5){
    //Find beta sets
    if(type  == "b" ||type =="all")FindSetsAndLoops(I,"b",fprefix, sthresh);
    if(type  == "ad" ||type =="all")FindSetsAndLoops(I,"ad",fprefix, sthresh);
    
    //PrintoutforFindBetasetLoops(&(*Gpb)[i],&(*Gpad)[i],i,&betasets_all,&betasetloops_all,&adsets_all,&adsetloops_all,prefix);
  }
  return;
}

//*********************************************************************************************//
//*********************************************************************************************//
//*********************************************************************************************//



void FindCellSets(Graph* Tree,Graph* setG,Graph *loopG, vector<vector<vertex_t> > * sets,vector<int>*loopcomp, double thresh){

  std::pair<vertex_iter, vertex_iter> lvp,vp ;//vp1,vp2, vtemp,vtemp2;
  std::pair<out_edge_iter, out_edge_iter> oep;
  edge_t e;
  vertex_t u,v;
  int vectorcount = 0;
  bool verbose=0,verbosemem=0;

  if(verbose)cout << "*****************Inside FindCellSets***************" << endl;cout.flush();
  if(verbose){
    cout << "setG is of type " ;
    if((*setG)[0].type == "b")cout << "b" <<endl;
    else cout << "ad" << endl;
    PrintoutGnuplot(loopG,"loopG");
  }
 

  //find all non-tree edges and store in F
  vector<edge_t > F; vectorcount++;if(verbosemem)cout << "Adding vector F: vectorcount =" << vectorcount << endl;cout.flush();//all non-tree edges  
  for(lvp = vertices(*loopG);lvp.first!=lvp.second;lvp.first++){
      u = *lvp.first;
      if(verbose)cout << "Checking cell " << (*loopG)[u].IsletCellNum << " with degree " << degree((*loopG)[u].CellNum,(*loopG)) <<  endl;
      for(oep = out_edges(u,*loopG);oep.first != oep.second; oep.first++){
          v = target(*oep.first,*loopG);
          if(verbose)cout << " with cell " << (*loopG)[v].IsletCellNum << endl;
          //Check if edge is in treenode
          bool b;
          tie(e,b) = edge(u,(*Tree)[v].CellNum,(*Tree));
          if(!b){
 	            if(verbose)cout << "  which does not have an edge in the tree" << endl;
 	            tie(e,b) = edge(u,v,(*loopG));
	            //Check if edge is already in F
	            bool inF = 0;
	            for(int j = 0; j < F.size();j++){
	              if(source(F[j],(*loopG)) == source(e,(*loopG)) && target(F[j],(*loopG)) == target(e,(*loopG))) inF = 1;
	              else if(source(F[j],(*loopG)) == target(e,(*loopG)) && target(F[j],(*loopG)) == source(e,(*loopG))) inF = 1;
	            }
	            if(!inF){

                  ///////// CHECK MY MANU FOR LONG EDGES ////////////////
                  //double test_dist = Distancexy((*loopG)[u].x,(*loopG)[u].y,(*loopG)[v].x,(*loopG)[v].y);
                  //if (test_dist > thresh){
                  //    cout << "Length of edge in F is " << test_dist << endl;
                  //}
                  ///////////////////////////////////////////////////////
                  F.push_back(e);
              }
          }
      }
  }

  //cout << "Press key to continue finding loops" << endl;
  //getchar();

  //cout << "Thresh is " << thresh << endl;
  //getchar();

  //Find all loops 
  vector<vector<vertex_t> > Tloops; vectorcount++;if(verbosemem)cout << "Adding vector Tloops: vectorcount =" << vectorcount << endl;cout.flush();
  //For each edge in F: find minimal path between vertices on tree
  for(int i = 0; i < F.size();i++){
    u = source(F[i],(*Tree));
    v = target(F[i],(*Tree));
    //if(verbose)
      //cout << "For edge ("  << (*Tree)[u].IsletCellNum << "," << (*Tree)[v].IsletCellNum << ")" << endl;
      //double test_dist = Distancexy((*Tree)[u].x,(*loopG)[u].y,(*Tree)[v].x,(*loopG)[v].y);
      //cout << "Length of edge in F is " << test_dist << endl;

    //New code replacing dijkstra
    vector<vertex_t> path;vectorcount++;if(verbosemem)cout << "Adding vector path: vectorcount =" << vectorcount << endl;cout.flush();
    bool nopath = ShortestWeightedPathBetweenTwoVertices(Tree, u, v,&path);
    //cout << "no path flag " << nopath << endl;
    
    // EDITED MY MANU TO ADD PATHS ONLY WHEN FOUND
    if (!nopath){

      // ADDED BY MANU
      // Check that each edge in path is within threshold
      int flag_add_path = 1;
      path.push_back(u);

      // CHECK PATH EDGE LENGTHS
      for(int k = 0; k < path.size()-1;k++){
          //cout << (*Tree)[path[k]].IsletCellNum << " ";
          double test_dist = Distancexy((*Tree)[path[k]].x,(*Tree)[path[k]].y,(*Tree)[path[k+1]].x,(*Tree)[path[k+1]].y);
          if (test_dist > thresh){
            flag_add_path = 0;
            break;
            //cout << "edge " << (*Tree)[path[k]].IsletCellNum << ", " << (*Tree)[path[k+1]].IsletCellNum << endl;
            //cout << "Length of edge in loop path is " << test_dist << endl;
            //getchar();
          }
      }
      if (flag_add_path){
          (Tloops).push_back(path);
      }
      //else{
      //  //cout << "Ignoring path because of long edges:" ;
      //  for(int k = 0; k < path.size();k++)cout << (*Tree)[path[k]].IsletCellNum << ",";
      //  //cout << endl;
      //}

    }
    ////////////////////////////////////////
    
    //Printout path
    if(verbose){
      cout << "The path is " ;
      for(int k = 0; k < path.size();k++)cout << (*Tree)[path[k]].IsletCellNum << " ";
      cout << "Press key to continue" << endl;
      getchar();
    }
    RemoveVector(&path); vectorcount--;if(verbosemem)cout << "Removing vector path: vectorcount =" << vectorcount << endl;cout.flush();
  }
  RemoveVector(&F);vectorcount--;if(verbosemem)cout << "Removing vector F: vectorcount =" << vectorcount << endl;cout.flush();

  //cout << "Press key to continue" << endl;
  //getchar();

  
  //Find cells with loops surrounding them
  vector<vertex_t> cellWithLoop; vectorcount++;
  if(verbosemem)cout << "Adding vector cellWithLoop: vectorcount =" << vectorcount << endl;cout.flush(); //cells that have at least 1 loop around them
  vector<vector<int> > loopsAroundCell; vectorcount++;
  if(verbosemem)cout << "Adding vector loopsAroundCell: vectorcount =" << vectorcount << endl;cout.flush(); //corresponding loop#(s) (from Tloops);

  //bool morethan1loop = 0;

  ///////////////////////////////////
  //// Printing Loop info, Manu
  ///////////////////////////////////
  //cout << "Printing all loops" << endl;
  //for(int i = 0; i < (Tloops).size();i++){

  //  cout << "Loop " << i << ":";
  //  for (int j = 0; j < Tloops[i].size(); j++){
  //        cout << (*Tree)[Tloops[i][j]].IsletCellNum << "," ;
  //  }

  //  cout << endl;
  //    
  //}
  ////exit(1);
  ///////////////////////////////////


  for(vp = vertices(*setG); vp.first != vp.second; vp.first++){

    //cout << "Cell:" << (*setG)[*vp.first].IsletCellNum << " is in loops:";

    double cellx = (*setG)[*vp.first].x;
    double celly = (*setG)[*vp.first].y;
    //Check loops
    for(int i = 0; i < (Tloops).size();i++){
      //if(verbose)cout << "i = " << i << endl;
      int wn1 = WindingNumber(Tree,&(Tloops)[i],cellx,celly);
      if(wn1 !=0){ //beta cell is inside loop
	        if(verbose)cout << " Cell " << (*setG)[*vp.first].IsletCellNum << " (" << (*setG)[*vp.first].x << "," << (*setG)[*vp.first].y << ") is inside loop " << endl;

          //cout << i << ",";

	        //check if b is already in list
	        if(cellWithLoop.size()==0 || cellWithLoop[cellWithLoop.size()-1] != *vp.first){
	          vector<int> temp;vectorcount++;
	          temp.push_back(i);
	          loopsAroundCell.push_back(temp);
	          RemoveVector(&temp);vectorcount--;
	          cellWithLoop.push_back(*vp.first);
	        }
	        else
	          loopsAroundCell[loopsAroundCell.size()-1].push_back(i);
      }
    }
    //cout << endl;

    //if(loopsAroundCell.size() > 0 && loopsAroundCell[loopsAroundCell.size()-1].size()>1)
    //  morethan1loop = 1;
    
  }



  RemoveVector2d(&Tloops);vectorcount--;if(verbosemem)cout << "Removing vector Tloops: vectorcount =" << vectorcount << endl;cout.flush();

  if(verbose){
    cout << "b's  with loops are " << endl;
    for(int m = 0; m < cellWithLoop.size();m++){
      cout << (*setG)[cellWithLoop[m]].IsletCellNum << ": " ;
      for(int n = 0; n < loopsAroundCell[m].size();n++)
	      cout << loopsAroundCell[m][n] << " ";
      cout << endl; 
    }
    cout << endl;
  }

  //count how many betas have the same loops
  vector<vector<vertex_t> > sameloops;vectorcount++;

  if(verbosemem)cout << "Adding vector sameloops: vectorcount =" << vectorcount << endl;cout.flush();

  int nsameloops = 0;
  //if(verbose){
  //  cout << "After creating sameloops, they are " << endl;
  //  for(int r = 0; r < sameloops.size();r++){
  //    for(int rr = 0; r < sameloops[r].size();rr++){
	//        cout << (*setG)[sameloops[r][rr]].IsletCellNum << " ";
  //    }
  //    cout << endl;
  //  }
  //  cout << endl;
  //}

  for(int m1 = 0;m1 < cellWithLoop.size();m1++){
    if(verbose)cout << "m1 = Cell " <<(*setG)[cellWithLoop[m1]].IsletCellNum << endl;
    bool insameloop = 0;
    for(int r = 0; r < sameloops.size();r++){
      for(int rr = 0; rr < sameloops[r].size();rr++){
	        if(cellWithLoop[m1] == sameloops[r][rr]){
	          if(verbose)cout << " m1 is in sameloops " << r << endl;
	          insameloop = 1;
	          r= sameloops.size();
	          break;
	        }
      }
    }	
    if(!insameloop){
      vector<vertex_t> temp; vectorcount++;
      temp.push_back(cellWithLoop[m1]);
      sameloops.push_back(temp);
      nsameloops++;
      RemoveVector(&temp);vectorcount--;
      for(int m2 =m1+1 ;m2 < cellWithLoop.size();m2++){
	        insameloop = 0;
	        for(int r = 0; r < sameloops.size();r++){
	          for(int rr = 0; rr < sameloops[r].size();rr++){
	            if(cellWithLoop[m2] == sameloops[r][rr]){
	              // cout << " m2 is in sameloops " << r << endl;
	              insameloop = 1;
	              r= sameloops.size();
	              break;
	            }
	          }
	        }
	        if(!insameloop){
	          if(verbose)cout << " m2 = Cell " << cellWithLoop[m2] << endl;
	          if(loopsAroundCell[m1].size()==loopsAroundCell[m2].size()){
	            if(verbose)cout << " m1 and m2 are same size" << endl;
	            bool samelp = 1;
	            for(int n = 0;n < loopsAroundCell[m1].size();n++){
	              if(loopsAroundCell[m1][n]!=loopsAroundCell[m2][n]){
	        	if(verbose)cout << "m1 and m2 are different" << endl;
	        	samelp =0;
	        	break;
	              }
	            }
	            if(samelp){
	              if(verbose)cout << "m1 and m2 are the same " << endl;
	              sameloops[sameloops.size()-1].push_back(cellWithLoop[m2]);
	              nsameloops++;
	            }
	            if(nsameloops == cellWithLoop.size()){
	              m1 = cellWithLoop.size();
	              m2 = cellWithLoop.size();
	            }  
	          }
	        }
      }
    }
  }

  //Printout same loops
  if(verbose){
    cout << "The same loops are " << endl;
    for(int r = 0; r < sameloops.size();r++){
      for(int rr = 0; rr < sameloops[r].size();rr++){
	        cout << sameloops[r][rr] << " ";
      }
      cout << endl;
    }
    cout << endl;
  }

  //cout << "Copying sameloops to beta sets" << endl;
  
  //Copy sameloops to betasets
  for(int r = 0; r < sameloops.size();r++){

    //cout << "set:" << r << " has b cells:";

    vector<vertex_t> temploop;vectorcount++;

    for(int rr = 0; rr < sameloops[r].size();rr++){
      temploop.push_back(sameloops[r][rr]);
      //cout << (*setG)[sameloops[r][rr]].IsletCellNum << ",";
    }
    //cout << endl;

    (*sets).push_back(temploop);

    // THIS LOOKS SUSPICIOUS -- MANU
    // WHY IS IT PUSHING ALWAYS COMP 0???
    //cout << "LOOP COMP NUM:" << (*loopG)[0].CompNum << endl;

    (*loopcomp).push_back((*loopG)[0].CompNum);

    RemoveVector(&temploop);vectorcount--;
    
  }

  //exit(1);

  if(verbose){
    for(int i =0;i < (*sets).size();i++){
      cout << "Entries of betaset " << i << " are:" << endl;
      for(int j = 0; j < (*sets)[i].size();j++){
	        cout << (*setG)[(*sets)[i][j]].IsletCellNum << ",";
      }
      cout << endl;
    }
    cout << endl;
  }
  //
  
  //exit(1);

  //remove cellWithLoop, loopsAroundCell, and sameloop
  RemoveVector(&cellWithLoop); vectorcount--;if(verbosemem)cout << "Removing vector cellWithLoop: vectorcount =" << vectorcount << endl;cout.flush();
  RemoveVector2d(&loopsAroundCell);vectorcount--;if(verbosemem)cout << "Removing vector loopsAroundCell: vectorcount =" << vectorcount << endl;cout.flush();
  RemoveVector2d(&sameloops);vectorcount--;if(verbosemem)cout << "Removing vector sameloops: vectorcount =" << vectorcount << endl;cout.flush();

  if(vectorcount>0)cout << "Problem with vectors! vectorcount = " << vectorcount << " in FindPrelimCellSet" << endl;
  return;
}	  

// //*********************************************************************************************//
// //*********************************************************************************************//
// //*********************************************************************************************//
int FindSetLoops(IsletGraph*I, Graph* setG,Graph *loopG, vector<vector<vertex_t> > *sets,vector<vector<vertex_t> >*setloops,vector<int>*loopcomp,string fprefix,string sthresh){
  std::pair<vertex_iter, vertex_iter> vp, vp1;
  edge_t e; 
  bool verbose=0,verbosemem=0;
  int graphcount = 0, vectorcount = 0;
  
  //cout << "Finding set loops..." << endl;
 
  //Make copy of loopG which will have points added to it from MakeGraphPlanar
  //Graph GpCopy; graphcount++;
  //CopyGraph(loopG,&GpCopy);
  if(verbose)cout << "Finished making copy" << endl;cout.flush();
  if(verbose)PrintoutGnuplotIslet(I,"tempprob","./", "labels");
  if(verbose)PrintoutGnuplot(setG,"tempprobset");
  if(verbose)PrintoutGnuplot(loopG,"tempprobloop");
    
  if(verbose){
    cout << "Inside FindSetLoops with " << (*sets).size()<< " sets consisting of sets " <<endl;
    for(int i = 0; i < (*sets).size();i++){
      cout << "set " << i << ": " ;
      for(int j = 0; j < (*sets)[i].size();j++)cout << (*setG)[(*sets)[i][j]].IsletCellNum << " " ;
      cout << endl;
    }
  }
  
  //exit(1);
  
  //Find closest set-loopcomp cell pair  
  for(int i = 0; i < (*sets).size();i++){

    vertex_t setmin = 0;
    vertex_t startv;
    if(verbose){
        cout << "**Checking set " << i << " containing:";
        for(int j = 0; j < (*sets)[i].size();j++) cout << (*setG)[(*sets)[i][j]].IsletCellNum << "," ;
        cout << endl;
    }

    //getchar();

    bool testset = 0;

    vector <vector<vertex_t> > PairsNotToUse;vectorcount++;
    
    //unordered_set <int> PairsNotToUse ;

    int pathnum = 1;
    vertex_t nextv;
    vector<vertex_t> Pathold;vectorcount++;
    if(verbosemem)cout << "Adding Pathold: vectorcount = " << vectorcount << endl; cout.flush();

    while (testset == 0)
    {
      double mindist = 100000;
      for (int j =0; j < (*sets)[i].size();j++)
      {

	        if(verbose)cout<< (*setG)[(*sets)[i][j]].IsletCellNum << " ";cout.flush();

	        for (vp = vertices((*loopG)); vp.first!=vp.second;vp.first++)
          {

            // THIS LOOKS SUSPICIOUS -- MANU
	            if((*loopG)[*vp.first].CompNum == (*loopcomp)[i])
              {
	                if(verbose) cout << "Checking cell " << ((*loopG))[*vp.first].IsletCellNum << endl;

	                if(out_degree(*vp.first,*loopG) > 0 || ((*loopG))[*vp.first].IsletCellNum > 10000)
                  {// wards off scragglers that might be located inside the ring
	                      double dist = Distancexy(((*setG))[(*sets)[i][j]].x,((*setG))[(*sets)[i][j]].y,(*loopG)[*vp.first].x,(*loopG)[*vp.first].y);

	                      if(dist < mindist)
                        {
		                        //Check to make sure pair is not in PairsNotToUse
		                        bool inPairs = 0;
                            //cout << "Size of PairsNotToUse " << PairsNotToUse.size() << endl;

		                        for(int p = 0; p < PairsNotToUse.size();p++){
		                            if(verbose)cout << "Checking pairs not to use " << (*setG)[PairsNotToUse[p][0]].IsletCellNum << " and " <<(*loopG)[PairsNotToUse[p][1]].IsletCellNum << endl; 
		                            if(PairsNotToUse[p][0] == (*sets)[i][j] && PairsNotToUse[p][1] == *vp.first){
		                                inPairs = 1;
		                                break;
		                            }
		                        }

                            //int minn = std::min(((*setG))[(*sets)[i][j]].IsletCellNum, (*loopG)[*vp.first].IsletCellNum);
                            //int maxx = std::max(((*setG))[(*sets)[i][j]].IsletCellNum, (*loopG)[*vp.first].IsletCellNum);

                            //cout << "Size of PairsNotToUse is " << PairsNotToUse.size() << endl;
                            //if (PairsNotToUse.find(insert_in_set(minn, maxx)) != PairsNotToUse.end()) {
                            ////std::cout << "element found." << std::endl;
                            //    inPairs = 1;
                            //}


		                        if(!inPairs){
		                            mindist = dist;
		                            setmin = (*sets)[i][j];
		                            startv=*vp.first;

		                            if(verbose) cout << "Making " << (*loopG)[*vp.first].IsletCellNum << " locally planar" << endl;cout.flush();

		                            MakeGraphLocallyPlanar(I,&(*loopG),*vp.first);

		                        }

	                      }

	                }

	            }

	        }

      }

      //if(verbose) 
      
      //cout << "   The closest pair are cell " <<((*setG))[setmin].IsletCellNum << " with loopcell " << (*loopG)[startv].IsletCellNum << endl;

      //getchar();
      
      //Find initial path
      double cellx = ((*setG))[setmin].x,celly = ((*setG))[setmin].y;
      if(verbose)cout << "Before finding initial path" << endl;

      FindInteriorPath(I,&(*loopG),startv, &Pathold, &nextv, pathnum,cellx,celly);

      if(verbose){
	        cout << "The initial path for set " << i << " is:";
	        for (int i = 0; i < Pathold.size();i++)
	          cout << ((*loopG))[Pathold[i]].IsletCellNum << ",";
	        cout << endl;
          cout << "Press key to continue" << endl;
          getchar();
      }
      
      //Check if initial path contains set
      bool containsSet = 1;
      for(int j = 0; j <(*sets)[i].size();j++){
	        //Calculate winding number
	        int wn1 = WindingNumber(&(*loopG),&Pathold,(*setG)[(*sets)[i][j]].x,(*setG)[(*sets)[i][j]].y);

	        if(wn1 ==0){ //set vertex is outside of loop
	          if(verbose) cout<< "Cell " << (*setG)[(*sets)[i][j]].IsletCellNum << " is not inside the loop" << endl;
	          containsSet = 0;
	          vector<vertex_t> temp;vectorcount++;
	          temp.push_back(setmin);temp.push_back(startv);
	          if(verbose) cout << "Adding " << (*setG)[setmin].IsletCellNum << "," << (*loopG)[startv].IsletCellNum << " to PairsNotToUse "<< endl;
            //getchar();

            //int minn = std::min((*setG)[setmin].IsletCellNum, (*loopG)[startv].IsletCellNum);
            //int maxx = std::max((*setG)[setmin].IsletCellNum, (*loopG)[startv].IsletCellNum);
            //PairsNotToUse.insert(insert_in_set(minn, maxx));

	          PairsNotToUse.push_back(temp);
            
	          RemoveVector(&temp);vectorcount--;
	          while(Pathold.size() > 0)Pathold.erase(Pathold.begin());
	          break;
	        }

      }

      if(containsSet)testset = 1;
    }

    RemoveVector2d(&PairsNotToUse);vectorcount--;
    //PairsNotToUse.clear();
   
    if(verbose){
      cout << "The path found before removing IPs is " ;
      for(int i = 0; i < Pathold.size();i++)cout << (*loopG)[Pathold[i]].IsletCellNum << " ";
      cout << endl;
    }
    
   
    //Add set loop
    (*setloops).push_back(Pathold);
    RemoveVector(&Pathold);vectorcount--;if(verbosemem)cout << "Removing Pathold: vectorcount = " << vectorcount << endl; cout.flush();
  }
  
  
  
  
  if(vectorcount != 0 || graphcount != 0) cout << "Problem in FindSetLoops with vectors and/or graphs! graphcount = " << graphcount << " and vectorcount = " << vectorcount << endl;
  cout << "Finished finding set loops" << endl;
   return 0;
}

//*********************************************************************************************//
//*********************************************************************************************//
//*********************************************************************************************//


bool FindSetsAndLoops(IsletGraph* I, string setType, string fprefix,string sthresh){
  edge_t e;
  std::pair<out_edge_iter, out_edge_iter> oep;
  std::pair<vertex_iter,vertex_iter>vp;
  bool verbose=0;
  int vectorcount = 0, graphcount = 0;
  Graph * setG,* loopG;
  vector<vector<vertex_t> > * loopComp, *sets;
  vector<vector<vertex_t> >*setloops;
  vector<int>*loopcompnum;
  string loopType;

  double thresh;

  if(setType == "b"){
    setG = &(*I).bGraph;
    sets = &(*I).betasets;
    setloops = &(*I).betasetLoops;
    loopG = &(*I).adGraph;
    loopComp = &(*I).adComponents;
    loopType = "ad";
    loopcompnum = &(*I).betasetloopcomp;

    thresh = (*I).admthresh;
    
  }
  if(setType == "ad"){
    setG = &(*I).adGraph;
    sets = &(*I).adsets;
    setloops = &(*I).adsetLoops;
    loopG = &(*I).bGraph;
    loopComp = &(*I).bComponents;
    loopType = "b";
    loopcompnum = &(*I).adsetloopcomp;

    thresh = (*I).bbthresh;

  }


  /////////////////////////////////
  // Printing Islet info, Manu
  /////////////////////////////////
  
  //for(vp = vertices(*setG);vp.first!=vp.second;vp.first++){
  //    cout << "Islet cell is:" << (*setG)[*vp.first].IsletCellNum
  //                << " type:" << (*setG)[*vp.first].type
  //                << " at x:" << (*setG)[*vp.first].x
  //                << " at y:" << (*setG)[*vp.first].y
  //                << endl;
  //}
  //exit(1);
  
  /////////////////////////////////



  //Graph * allgraph;
  //allgraph = &(*I).bGraph;
  //for(vp = vertices(*allgraph);vp.first!=vp.second;vp.first++){
  //    cout << "Islet cell is:" << (*allgraph)[*vp.first].IsletCellNum
  //                << " type:" << (*allgraph)[*vp.first].type
  //                << " at x:" << (*allgraph)[*vp.first].x
  //                << " at y:" << (*allgraph)[*vp.first].y
  //                << endl;
  //}

  //allgraph = &(*I).adGraph;
  //for(vp = vertices(*allgraph);vp.first!=vp.second;vp.first++){
  //    cout << "Islet cell is:" << (*allgraph)[*vp.first].IsletCellNum
  //                << " type:" << (*allgraph)[*vp.first].type
  //                << " at x:" << (*allgraph)[*vp.first].x
  //                << " at y:" << (*allgraph)[*vp.first].y
  //                << endl;
  //}

  //exit(1);
  

  /////////////////////////////////
  
  //cout << "Press key to find " << setType << " sets... "  << endl;

  //if(verbose)
  //cout << "Islet has " << (*loopComp).size() << " components of type " << loopType << endl;

  for(int c = 0; c < (*loopComp).size();c++){

    if((*loopComp)[c].size()>2){

        Graph Tree,loopCompG;graphcount +=2;

        if(verbose)cout << "****Islet = " << (*I).IsletNum  << " and C =  " << c << "************" << endl;

        //Add vertices to Trees and GpadComp
        CopyComponentVerticesFromGraph(I,c,&Tree,loopType);
        CopyComponentVerticesFromGraph(I,c,&loopCompG,loopType);
        //Printout vertices in Tree and loopcomp (Gpadcomp)

        //// Print ad component vertices
	      //for(vp = vertices(loopCompG);vp.first!=vp.second;vp.first++)
	      //  cout << "loop component:" << c<<  " Vertex:" << (loopCompG)[*vp.first].IsletCellNum << endl;

        if(verbose){
	          cout << "For component " << c << ":";
	          for(int i =0; i < num_vertices(Tree);i++)cout << Tree[i].IsletCellNum << " ";
        }
        //Add edges to GpadComp
        CopyComponentEdgesFromGraph(I,c,&loopCompG,loopType);

        if(verbose){
	          cout << "loopCompGraph contains: " << endl;
	          for(vp = vertices(loopCompG);vp.first!=vp.second;vp.first++)
	            cout << "Vertex " << (loopCompG)[*vp.first].IsletCellNum << " with degree " << out_degree(*vp.first,loopCompG) << endl;
            PrintoutGnuplotIslet(I,"tempIslet","./","labels");
        }
        

        //Find a spanning tree of Gp
        SpanningTree(&loopCompG,&Tree);

        if(verbose)cout << "Created spanning tree" << endl;
        if(c==1)verbose=0;
        //if(verbose){
	          //PrintoutGnuplot(&Tree,"tempTree", "./", "abd");
	          //PrintoutGnuplot(&loopCompG,"tempComp","./","abd");
        //}
        verbose=0;
        //Find loops around beta cells of original spanning tree
        // int origsetsize = (*sets).size();
        
        FindCellSets(&Tree,setG,&loopCompG,sets,loopcompnum, thresh);


        //int newsetsize = (*sets).size();
        //if(newsetsize > origsetsize){
        //for(int i = origsetsize; i < newsetsize;i++)
        //  (*loopcompnum)[i] = c;
        //}

        //Printout Prelim beta set
        if((*sets).size()>0 and verbose){
	          cout << "The betasets are " << endl;
	          for(int r = 0; r < (*sets).size();r++){
	            for(int rr = 0; rr < (*sets)[r].size();rr++){
	              cout << (*setG)[(*sets)[r][rr]].IsletCellNum << " ";
	            }
	            cout << endl;
	          }
	          cout << endl;
        }
        RemoveGraph(&Tree);graphcount--;
        RemoveGraph(&loopCompG);graphcount--;
    }
  }

  //if((*sets).size()>0 ){
	//    cout << "The betasets are:" << endl;
	//    for(int r = 0; r < (*sets).size();r++){
	//      for(int rr = 0; rr < (*sets)[r].size();rr++){
	//        cout << (*setG)[(*sets)[r][rr]].IsletCellNum << ",";
	//      }
	//      cout << endl;
	//    }
	//    cout << endl;
  //}

  //exit(1);

  verbose=0;
  if(verbose){
    cout << "Finished Finding sets!" << endl;
    cout.flush();
    cout << "The sets are " << endl;
    for(int i = 0; i < (*sets).size();i++){
      cout << "Comp " << (*loopcompnum)[i] << ":";
      for(int j = 0; j < (*sets)[i].size();j++){
	        cout << (*setG)[(*sets)[i][j]].IsletCellNum << " ";
      }
    }
  }

  verbose=0;
  //cout << "Finding " << setType << " loops..." << endl;
  //cout << "Press key to FindSetLoops" << endl;
  //getchar();
  if((*sets).size() >0 ){
    //Find betaset loops
    //    vector<vertex_t> LeaveIPs;vectorcount++;
    bool FBSL_Failed=FindSetLoops(I,setG,loopG,sets,setloops,loopcompnum,fprefix,sthresh);
    if(FBSL_Failed)return 1;
  }


  ofstream gnupb;
  string gnupbfile = "../Loops/" + fprefix + "." + loopType + "mantles";
  gnupb.open(gnupbfile.c_str());
  
  //if(verbose){
    //cout << "Finished Finding sets and setloops!" << endl;
    //cout.flush();
    //cout << "The sets are " << endl;
    for(int i = 0; i < (*sets).size();i++){
      for(int j = 0; j < (*sets)[i].size();j++){
	        //cout << (*setG)[(*sets)[i][j]].IsletCellNum << " ";
          gnupb << (*setG)[(*sets)[i][j]].IsletCellNum << ",";
      }
      //cout << endl << " with loop " ;
      gnupb << ":";
      for(int j = 0; j < (*setloops)[i].size();j++){
          // MANU : Do not know how to get rid of IP cells easily.
	        if((*setloops)[i][j] < 100000){
              //cout << (*loopG)[(*setloops)[i][j]].IsletCellNum << ",";
              gnupb << (*loopG)[(*setloops)[i][j]].x << "," << (*loopG)[(*setloops)[i][j]].y << ",";
          }
	        else{
              //cout << (*loopG)[(*setloops)[i][j]].x << "," << (*loopG)[(*setloops)[i][j]].y << ",";
              gnupb << (*loopG)[(*setloops)[i][j]].x << "," << (*loopG)[(*setloops)[i][j]].y << ",";
          }
      }
      //cout << endl << endl;
      gnupb << endl;
    }
  //}

    gnupb.close();

  if(graphcount>0 || vectorcount < 0)cout << "Problem with vectors and/or graphs! vectorcount = " << vectorcount << " and graphcount = " << graphcount << endl;
  return 0;
}

//*********************************************************************************************//
//*********************************************************************************************//
//*********************************************************************************************//
 
double isLeft( double P0x,double P0y, double P1x,double P1y, double P2x, double P2y )
 {
   //isLeft(): tests if a point is Left|On|Right of an infinite line.
   //    Input:  three points P0, P1, and P2
   //    Return: >0 for P2 left of the line through P0 and P1
   //            =0 for P2  on the line
   //            <0 for P2  right of the line
   //    See: Algorithm 1 "Area of Triangles and Polygons"
   //algorithm from http://geomalgorithms.com/a03-_inclusion.html
   return ( (P1x - P0x) * (P2y - P0y)
            - (P2x -  P0x) * (P1y - P0y) );
 }

//*********************************************************************************************//
//*********************************************************************************************//
//*********************************************************************************************//

void SpanningTree(Graph * Gp,  Graph * Tree){ //Creates a Spanning Tree of a 1-component graph
  std::pair<vertex_iter, vertex_iter> vp;
  std::pair<out_edge_iter, out_edge_iter> oep;
  edge_t e;
  std::pair<edge_iter, edge_iter> ep;
  vertex_t u,v;
  int vectorcount = 0;
  bool verbose=0,verbosemem=0;
  
  vector<vertex_t> queue;vectorcount++;if(verbosemem)cout << "Adding vector queue: vectorcount =" << vectorcount << endl;cout.flush();
  //Add 0 to queue
  vp = vertices((*Gp));
  queue.push_back(*vp.first);
  while(queue.size() >0){
    u = queue[0];
    if(verbose)cout << "For cell " << (*Gp)[u].IsletCellNum << endl;
    for(oep = out_edges(u,(*Gp));oep.first!=oep.second;oep.first++){
      v = target(*oep.first,(*Gp));
      if(verbose)cout << " with edge to cell " << (*Gp)[v].IsletCellNum << " with degree " << out_degree(v,(*Gp))<<  endl;
      //if degree of out edge target on tree = 0, add edge between target and source
      //bool inTree = 1;
      //      if((*Gp)[v].CellNum > 99999){
       
      //if(out_degree((*Gp)[v].CellNum,(*Tree)) ==0){
      // EDITED BY MANU
      if(out_degree(v,(*Tree)) ==0){

	        //inTree = 0;
 	        if(verbose)cout << " Adding edge " << endl;cout.flush();

 	        //update tree
 	        bool b;
 	        //tie(e,b) = add_edge((*Gp)[u].CellNum,(*Gp)[v].CellNum,(*Tree));
          
          // EDITED BY MANU
 	        tie(e,b) = add_edge(u,v,(*Tree));

	        double dist = Distancexy((*Gp)[u].x,(*Gp)[u].y,(*Gp)[v].x,(*Gp)[v].y);
	        InitializeEdge(Tree,e,dist);
          	//add to queue
 	        queue.push_back(v);
      }
    }
    queue.erase(queue.begin());
    if(verbose){
      cout << "queue contains" << endl;
      for(int i = 0; i < queue.size();i++) cout << (*Gp)[queue[i]].IsletCellNum << " ";
      cout << endl;
    }
  }

  if(verbose)cout << "Finished creating spanning tree" << endl;cout.flush();
  
  //Clear vectors
  RemoveVector(&queue);vectorcount--;if(verbosemem)cout << "Removing vector queue: vectorcount =" << vectorcount << endl;cout.flush();
  
  if(vectorcount > 0) cout << "Problem with vectors! vectorcount = " << vectorcount << " in Spanning Tree " << endl;
  return;
}

///*********************************************************************************************//
//*********************************************************************************************//
//*********************************************************************************************//

int WindingNumber(Graph* Tree, vector<vertex_t>* Loop,double Px,double Py ){
  //Gp = Tree;
  // wn_PnPoly(): winding number test for a point in a polygon
  //      Input:   P = a point,
  //               V[] = vertex points of a polygon V[n+1] with V[n]=V[0]
  //      Return:  wn = the winding number (=0 only when P is outside)
  //algorithm from http://geomalgorithms.com/a03-_inclusion.html
  int    wn = 0;    // the  winding number counter
  
    // loop through all edges of the polygon
    //cout << "Inside winding number " << endl;
    for (int i=0; i<(*Loop).size()-1; i++) {   // edge from V[i] to  V[i+1]
      if ((*Tree)[(*Loop)[i]].y <= Py) {          // start y <= P.y
	if ((*Tree)[(*Loop)[i+1]].y  > Py)      // an upward crossing
	  if (isLeft((*Tree)[(*Loop)[i]].x ,(*Tree)[(*Loop)[i]].y,(*Tree)[(*Loop)[i+1]].x ,(*Tree)[(*Loop)[i+1]].y,Px,Py) > 0)  // P left of  edge
	    ++wn;            // have  a valid up intersect
      }
      else {                        // start y > P.y (no test needed)
	if ((*Tree)[(*Loop)[i+1]].y  <= Py)     // a downward crossing
	  if (isLeft((*Tree)[(*Loop)[i]].x ,(*Tree)[(*Loop)[i]].y,(*Tree)[(*Loop)[i+1]].x ,(*Tree)[(*Loop)[i+1]].y,Px,Py ) < 0)  // P right of  edge
	    --wn;            // have  a valid down intersect
      }
    }
    //cout << "Finished winding number " << endl;
  return wn;
}

///*********************************************************************************************//
//*********************************************************************************************//
//*********************************************************************************************//








