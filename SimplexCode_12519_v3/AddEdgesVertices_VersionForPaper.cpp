  #include "./GraphNN_VersionForPaper.h"
 #define PI 3.14159

using namespace boost;

bool AddEdges(Graph* Gp,  string celltype, float thresh, double shadowthresh){
  ////This program first adds all neighboring cells of given cell type that are within radius///
  ////Then checks each neighboring cell to see if there is a cell in between.  
  ////If so, it removes the newly added edge//

  std::pair<vertex_iter, vertex_iter>vp1,vp2;
  edge_t e;
  std::pair<out_edge_iter, out_edge_iter> oep;
  std::pair<edge_iter, edge_iter> ep;
  vertex_t u,v;
  
  int vectorcount = 0;
  bool verbose=0,verbosemem=0;//checked mem
  if(verbose ||verbosemem)cout << "*****Inside AddEdges: Celltype is " << celltype << endl;
  
  vector<vector<vertex_t> > EdgesErased;vectorcount++;if(verbosemem)cout << "Adding vector EdgesErased vectorcount = " << vectorcount << endl;cout.flush();
  
  for (vp1 = vertices(*Gp);vp1.first != vp1.second; vp1.first++){
    if((celltype=="abd" || celltype == "adb" )||
       (celltype =="aa" && (*Gp)[*vp1.first].type == "a") ||
       (celltype =="ab" && ((*Gp)[*vp1.first].type == "a" ||(*Gp)[*vp1.first].type == "b")) ||
       ((celltype =="ad" || celltype == "adm") && ((*Gp)[*vp1.first].type == "a" ||(*Gp)[*vp1.first].type == "d")) ||
       (celltype =="bb" && (*Gp)[*vp1.first].type == "b") ||
       (celltype =="bd" && ((*Gp)[*vp1.first].type == "b" ||(*Gp)[*vp1.first].type == "d"))||
       (celltype =="dd" && (*Gp)[*vp1.first].type == "d")){

      vp2 = vertices(*Gp);
      
      //Create list of neighbors of vp1 based on neighborhood threshold
      std::vector <VertexDistance> Neighbors;vectorcount++;if(verbosemem)cout << "Adding vector Neighbors vectorcount = " << vectorcount << endl;cout.flush();
      while(vp2.first != vp2.second){
	        if(*vp1.first!=*vp2.first){
	          double Celldist = Distancexy((*Gp)[*vp1.first].x,(*Gp)[*vp1.first].y,(*Gp)[*vp2.first].x,(*Gp)[*vp2.first].y);
	          	  
	          //Create Neighbors  
	          if(Celldist < thresh){ 
	            VertexDistance A;
	            A.vertex = *vp2.first;
	            A.distance = Celldist;
	            Neighbors.push_back(A);
	          }
	        }//done creating individual neighbor
	        ++vp2.first; 
      }//while vp2
      
      //Sort Neighbors
      std::sort(Neighbors.begin(), Neighbors.end(), compareByDistance);
      
      //Check if existing edges is greater than the number of neighbors found
      if(verbose && out_degree(*vp1.first,(*Gp)) > Neighbors.size()){ 
	        cout << "\n***Problem with edges and Neighbors***\nThe edges are" << endl;
	        for(oep = out_edges(*vp1.first,(*Gp));oep.first!=oep.second;oep.first++){
	          cout << (*Gp)[source(*oep.first,(*Gp))].IsletCellNum << "\t" << (*Gp)[target(*oep.first,(*Gp))].IsletCellNum << endl;
	        }
	        cout << "The neighbors are " << endl;
	        for(int j = 0; j < Neighbors.size();j++)cout <<(*Gp)[Neighbors[j].vertex].IsletCellNum << " " ;
	        cout << endl;
	        return 1;
      }
      
      //Printout Neighbors to check
      if(verbose && Neighbors.size() > 0){
	        cout << "For cell " << (*Gp)[*vp1.first].IsletCellNum << " (" << (*Gp)[*vp1.first].x << "," << (*Gp)[*vp1.first].y << ") there are " << Neighbors.size() << " sorted neighbors. They are " << endl;
	        for(int j = 0; j < Neighbors.size();j++)cout <<(*Gp)[Neighbors[j].vertex].IsletCellNum << " (" <<(*Gp)[Neighbors[j].vertex].x << "," << (*Gp)[Neighbors[j].vertex].y << ") with distance " << Neighbors[j].distance << endl;
	//PrintoutGnuplot(Gp,"GraphVerts");
      }
      
      //Add a temporary edge between vp1 and all neighbors 
      for(std::vector<VertexDistance>::iterator it = Neighbors.begin(); it != Neighbors.end(); it++){
	        bool b;
	        boost::tie(e,b) = edge(*vp1.first,(*it).vertex,(*Gp));
	        if(!b){
	          if(verbose)cout << "Celltype = " << celltype << " and it.vertex is of type " << (*Gp)[(*it).vertex].type << endl;

	          if(celltype=="abd" ||
	             (celltype =="bb" && (*Gp)[*vp1.first].type == "b" && (*Gp)[(*it).vertex].type == "b")  ||
 	             (celltype == "adm"&& ((*Gp)[*vp1.first].type == "a" ||(*Gp)[*vp1.first].type == "d") &&
 	              ((*Gp)[(*it).vertex].type == "a" ||(*Gp)[(*it).vertex].type == "d"))  || 
	             (celltype == "adb"&& (((*Gp)[*vp1.first].type == "a" &&(*Gp)[(*it).vertex].type == "b") ||
	        			   ((*Gp)[*vp1.first].type == "d" &&(*Gp)[(*it).vertex].type == "b") ||
	        			   ((*Gp)[*vp1.first].type == "b" &&(*Gp)[(*it).vertex].type == "d") ||
	        			   ((*Gp)[*vp1.first].type == "b" &&(*Gp)[(*it).vertex].type == "a")) )  ||
	             (celltype =="aa" && (*Gp)[*vp1.first].type == "a" && (*Gp)[(*it).vertex].type == "a") ||
	             (celltype =="ab" && ((*Gp)[*vp1.first].type == "a" ||(*Gp)[*vp1.first].type == "b") &&
	              ((*Gp)[(*it).vertex].type == "a" ||(*Gp)[(*it).vertex].type == "b") && 
	              ((*Gp)[*vp1.first].type !=(*Gp)[(*it).vertex].type))||
	             (celltype =="ad" && ((*Gp)[*vp1.first].type == "a" ||(*Gp)[*vp1.first].type == "d") &&
	              ((*Gp)[(*it).vertex].type == "a" ||(*Gp)[(*it).vertex].type == "d") && 
	              ((*Gp)[*vp1.first].type !=(*Gp)[(*it).vertex].type))||
	             (celltype =="bd" && ((*Gp)[*vp1.first].type == "b" ||(*Gp)[*vp1.first].type == "d") &&
	              ((*Gp)[(*it).vertex].type == "b" ||(*Gp)[(*it).vertex].type == "d") && 
	              ((*Gp)[*vp1.first].type !=(*Gp)[(*it).vertex].type))||
	             (celltype == "dd" && (*Gp)[*vp1.first].type == "d" && (*Gp)[(*it).vertex].type == "d") 
	             ){
	            boost::tie(e,b) = add_edge(*vp1.first,(*it).vertex,(*Gp));
	            if(verbosemem)cout << "Adding temporary edge between " << (*Gp)[*vp1.first].IsletCellNum << " (type " << (*Gp)[*vp1.first].type << ") and " << (*Gp)[(*it).vertex].IsletCellNum << "(type " << (*Gp)[(*it).vertex].type << ")" << endl; 
	            InitializeEdge(Gp,e,(*it).distance);
	            
	          }
	        }
      }

      if(shadowthresh > 0){
	        //Edit list by angles
	        int ii =0;
	        bool end_loop = false;
	        int Num_erased = 0;
	        while (end_loop == false && Neighbors.size()>1){

	          std::vector<VertexDistance>::iterator it1 = Neighbors.begin();
	          it1+=ii;
	          	  
	          std::vector<VertexDistance>::iterator itemp = it1;
	          itemp++;
	          std::vector<VertexDistance>::iterator it2 =itemp; 

	          while(it2 !=Neighbors.end()){
	            float x1 = (*Gp)[*vp1.first].x, x2 = (*Gp)[(*it1).vertex].x, x3 = (*Gp)[(*it2).vertex].x;
	            float y1 = (*Gp)[*vp1.first].y, y2 = (*Gp)[(*it1).vertex].y, y3 = (*Gp)[(*it2).vertex].y;
	            if(verbose)cout<< "Checking " << (*Gp)[(*it1).vertex].IsletCellNum << " with " << (*Gp)[(*it2).vertex].IsletCellNum << "\n";
	            
	            //Calculate Angles
	            float alpha2 = atan(shadowthresh/(*it1).distance);
	            float alpha3 = atan(shadowthresh/(*it2).distance);
	            float beta2 = 0.0,beta3 = 0.0;
	            if(x2>x1)beta2 = atan((y2-y1)/(x2-x1));//Q1 and 4
	            else if(y2>y1)beta2 = atan((y2-y1)/(x2-x1))+ PI;//Q2
	            else beta2 = atan((y2-y1)/(x2-x1))- PI;//Q3
	            if(x3>x1) beta3 = atan((y3-y1)/(x3-x1));
	            else if(y3>y1)beta3 = atan((y3-y1)/(x3-x1))+PI; 
	            else beta3 = atan((y3-y1)/(x3-x1))-PI;
	            beta2 = mod(beta2,2*PI);
	            beta3 = mod(beta3,2*PI);
	            vector<double> beta2shad(2,0),beta3shad(3,0);vectorcount+=2;//if(verbosemem)cout << "Adding vectors beta2shad, beta3shad: vectorcount = " << vectorcount << endl;cout.flush();
	            beta2shad[0] = mod(beta2-alpha2,2*PI);beta2shad[1] = mod(beta2+alpha2,2*PI);
	            beta3shad[0] = mod(beta3-alpha3,2*PI);beta3shad[1] = mod(beta3+alpha3,2*PI);
	            
	            if(verbose)cout << "v2's shadow is (" << beta2shad[0] << "," <<  beta2shad[1]  << ") and v3's shadow is (" << beta3shad[0] << "," <<  beta3shad[1] << ")" << endl; 
	            bool checkerase = 0;
	            if(beta2shad[0]<beta2shad[1]){//segment doesn't wrap
	              if(beta3shad[0]<beta3shad[1]){//segment doesn't wrap
	        	        if(beta3shad[0]> beta2shad[0] && beta3shad[0] < beta2shad[1]){//shadows intersect 
	        	          checkerase = 1;
	        	          if(verbose)cout << "notchecked1 " << endl;
	        	        }
	        	        if(beta3shad[1]> beta2shad[0] && beta3shad[1] < beta2shad[1]){//shadows intersect 
	        	          checkerase = 1;
	        	          if(verbose)cout << "notchecked2 " << endl;
	        	        }
	              }
	              else{//beta3 wraps but beta2 doesn't
	        	        if(beta3shad[0] <beta2shad[1]){
	        	          checkerase = 1;
	        	          if(verbose)cout << "notchecked3 " << endl;
	        	        }
	        	        if(beta3shad[1] > beta2shad[0]){
	        	          checkerase = 1;
	        	          if(verbose)cout << "notchecked4 " << endl;
	        	        }
	              }
	            }
	            else{//beta2 wraps
	              if(beta3shad[0]<beta3shad[1]){//beta3 doesn't wrap
	        	        if(beta3shad[0]<beta2shad[1]){
	        	          checkerase = 1;
	        	          // cout << "notchecked5 " << endl;
	        	        } 
	        	        if(beta3shad[1]>beta2shad[0]){
	        	          checkerase = 1;
	        	          // cout << "notchecked6 " << endl;
	        	        }
	              }
	              else{//both beta2 and beta3 wrap
	        	        if(beta3shad[0]>beta2shad[0]){
	        	          checkerase = 1;
	        	          // cout << "notchecked7 " << endl;
	        	        } 
	        	        if(beta3shad[1]<beta2shad[1]){
	        	          checkerase = 1;
	        	          // cout << "notchecked8 " << endl;
	        	        } 
	              }
	            }
	            
	            //remove beta2 and beta3shad vectors
	            RemoveVector(&beta2shad);vectorcount--;//if(verbosemem)cout << "Removing vector beta2shad vectorcount = " << vectorcount << endl;cout.flush();
	            RemoveVector(&beta3shad);vectorcount--;//if(verbosemem)cout << "Removing vector beta3shad vectorcount = " << vectorcount << endl;cout.flush();
	            
	            if(checkerase ==1){
	              vector<vertex_t> temp(2);vectorcount++;if(verbosemem)cout << "Adding vector temp: vectorcount = " << vectorcount << endl;cout.flush();
	              temp[0] = *vp1.first;
	              temp[1] = (*it2).vertex;
	              if(verbose)cout<< "Adding edge between " << (*Gp)[temp[0]].IsletCellNum << " and " << (*Gp)[temp[1]].IsletCellNum << "to be erased \n";
	              EdgesErased.push_back(temp);
	              RemoveVector(&temp); vectorcount--;
                if(verbosemem)cout << "Removing vector temp: vectorcount = " << vectorcount << endl;cout.flush();
	              it2++;
	            }
	            else{
	              it2++;
	            }
	          }
	          ii++;
	          if((*it1).vertex ==(Neighbors.back()).vertex){
	            end_loop = true;
	          }
	          else if(Neighbors.size()<2 || (*(++it1)).vertex ==(Neighbors.back()).vertex){
	            end_loop=true;
	          }
	        }
      }

      RemoveVector(&Neighbors);vectorcount--;

      if(verbosemem)cout << "Removing vector Neighbors: vectorcount = " << vectorcount << endl;cout.flush();

      if(verbose){
	        cout << "There are " <<out_degree(*vp1.first,(*Gp)) << " edges remaining. They are:" << endl;
	        for(oep = out_edges(*vp1.first,(*Gp));oep.first != oep.second;oep.first++){
	          vertex_t u = source(*oep.first,(*Gp));
	          vertex_t v = target(*oep.first,(*Gp));
	          cout << (*Gp)[u].IsletCellNum << " " << (*Gp)[v].IsletCellNum << endl;
	        }
      }
    }
  }
  //Erase edges
  for(int i = 0; i < EdgesErased.size();i++){
    bool b;
    tie(e,b) = edge(EdgesErased[i][0],EdgesErased[i][1],(*Gp));
    if(b){
      if(verbose)cout<< "Erasing edge between " << (*Gp)[EdgesErased[i][0]].IsletCellNum << " and " << (*Gp)[EdgesErased[i][1]].IsletCellNum << "\n";
      ClearEdge(Gp,e);
      remove_edge(e,(*Gp));
    }
  }
  //Remove EdgesErased vector
  RemoveVector2d(&EdgesErased);vectorcount--;if(verbosemem)cout << "Removing vector EdgesErased: vectorcount = " << vectorcount << endl;cout.flush();
  
  //Calculate maximal number of edges based on threshold
  double cell_degrees = 2.0*tan((double)shadowthresh/(double)thresh);
  int npossedges = (int)(2.0*PI/cell_degrees) +1;
  for(vp1 = vertices(*Gp);vp1.first!=vp1.second;vp1.first++){
    if(out_degree(*vp1.first,(*Gp))> npossedges){
      cout <<  "Warning with edges" << endl;
      cout << " The max number of edges is " << npossedges << ": cell " <<(*Gp)[*vp1.first].IsletCellNum << "  has " << out_degree(*vp1.first,(*Gp)) << " edges \n" ;
      cout << " The cell degrees are " << cell_degrees << " with shadowthresh = " << shadowthresh << " and thresh = " << thresh << endl;
    }
  }
  if(verbose){
    //Printout edges
    cout << "Edges: "  << endl;
    for(ep = edges(*Gp); ep.first != ep.second; ep.first++)
      cout << (*Gp)[source(*ep.first,(*Gp))].IsletCellNum << "\t" <<(*Gp)[target(*ep.first,(*Gp))].IsletCellNum <<  endl;
  }
  if(verbosemem)cout << "*****Finished AddEdges" << endl;
  if(vectorcount!= 0)cout << "Problem with vectors in AddEdges! vectorcount = " << vectorcount << " in AddEdges()" << endl;
  return 0;
}

// //*********************************************************************************************//
// //*********************************************************************************************//
// //*********************************************************************************************//

void AddVertices(vector<IsletGraph>* Islets, string filename, string celltype, string AVopt){
  bool verbose=0,verbosemem=0;//checked
  vertex_t u;
  int graphcount = 0,isletcount = 0;
  
  if(verbosemem)cout << "*****Inside AddVertices" << endl;

  //*****************************************************************************//
  //                         Read in each line of data                           //
  //*****************************************************************************//
  //Add vertices to graph
  ifstream fin(filename.c_str());

  int num = 0, IsletNumber = 0,IsletCellNumber=0;
  while(!fin.eof()){
    int IsletNumberin = 0,CellNumberin = 0;
    double xin = 0,yin = 0,IsletAreain = 0;
    string ctypein;
    if(AVopt == "rgc"  || AVopt == "121113" ||AVopt == "abd")
      fin >>IsletNumberin >> CellNumberin >> xin >> yin >> ctypein;
    else
      fin >>IsletNumberin >> CellNumberin >> xin >> yin >> ctypein  >> IsletAreain;
    if(fin.eof())break;

    //cout << IsletNumberin << CellNumberin << xin << yin << ctypein  << IsletAreain << endl;

    //*****************************************************************************//
    //                        Add cells to new islets                              //
    //*****************************************************************************//
    //If new Islet
    if(IsletNumberin != IsletNumber){
      IsletCellNumber = 0;
      //beta cells
      if (((AVopt == "rgc" && ctypein =="gfp") || (AVopt == "121113"&& ctypein == "11") ||(AVopt == "abd"&& ctypein == "b") ||
	   (AVopt == "321Area"&& ctypein == "2") || (AVopt == "rgcArea" && ctypein == "gfp")|| (AVopt == "abdArea" && ctypein == "b") ) &&
	  (celltype == "ab"||celltype =="bb" || celltype == "bd" ||celltype == "abd" )){
	
	      num = 0;
	      Graph gnew;graphcount++;if(verbosemem)cout << "Adding graph gnew: graphcount = " << graphcount << endl; cout.flush();
	      //add vertex
	      u=add_vertex(gnew);
	      InitializeVertex(&gnew,u,xin,yin,num,"b",CellNumberin);
	      num++;

 	      IsletGraph G;isletcount++;if(verbosemem)cout << "Adding isletgraph G: isletcount = " << isletcount << endl; cout.flush();
	      InitializeIslet(&G,&gnew,IsletNumberin,IsletAreain);
	      RemoveGraph(&gnew);graphcount--;if(verbosemem)cout << "Removing graph gnew: graphcount = " << graphcount << endl; cout.flush();
	      G.nbeta = 1;
        	(*Islets).push_back(G);
	      RemoveIsletGraph(&G);isletcount--;if(verbosemem)cout << "Removing isletgraph G: isletcount = " << isletcount << endl; cout.flush();
 	      IsletNumber = IsletNumberin;
      }
      
      //alpha cells
      else if (((AVopt == "rgc" && ctypein =="rfp") || (AVopt == "121113"&& ctypein == "12") ||(AVopt == "abd"&& ctypein == "a") ||
		(AVopt == "321Area"&& ctypein == "3") || (AVopt == "rgcArea" && ctypein == "rfp")|| (AVopt == "abdArea" && ctypein == "a") ) &&
	       (celltype == "aa"||celltype =="ab" || celltype == "ad" ||celltype == "abd" || celltype == "adm")){ 
	      num = 0;
	      Graph gnew;graphcount++;if(verbosemem)cout << "Adding graph gnew: graphcount = " << graphcount << endl; cout.flush();
	      u=add_vertex(gnew);
	      InitializeVertex(&gnew,u,xin,yin,num,"a",CellNumberin);
	      num++;
	      
	      IsletGraph G;isletcount++;if(verbosemem)cout << "Adding isletgraph G: isletcount = " << isletcount << endl; cout.flush();
	      InitializeIslet(&G,&gnew,IsletNumberin,IsletAreain);
	      RemoveGraph(&gnew);graphcount--;if(verbosemem)cout << "Removing graph gnew: graphcount = " << graphcount << endl; cout.flush();
	      G.nalpha = 1;
	      (*Islets).push_back(G);
	      RemoveIsletGraph(&G);isletcount--;if(verbosemem)cout << "Removing isletgraph G: isletcount = " << isletcount << endl; cout.flush();
	      IsletNumber = IsletNumberin;
      }

      //delta cells
      else if (((AVopt == "rgc" && ctypein =="cy5") || (AVopt == "121113"&& ctypein == "13") ||(AVopt == "abd"&& ctypein == "d") ||
		(AVopt == "321Area"&& ctypein == "1") || (AVopt == "rgcArea" && ctypein == "cy5")|| (AVopt == "abdArea" && ctypein == "d")) &&
	       (celltype == "ad"||celltype =="bd" || celltype == "dd" ||celltype == "abd" || celltype == "adm" )){//delta cells
	num = 0;
	Graph gnew;graphcount++;if(verbosemem)cout << "Adding graph gnew: graphcount = " << graphcount << endl; cout.flush();
	u=add_vertex(gnew);
	InitializeVertex(&gnew,u,xin,yin,num,"d",CellNumberin);
	num++;
	
	IsletGraph G;isletcount++;if(verbosemem)cout << "Adding isletgraph G: isletcount = " << isletcount << endl; cout.flush();
	InitializeIslet(&G,&gnew,IsletNumberin,IsletAreain);
	RemoveGraph(&gnew);graphcount--;if(verbosemem)cout << "Removing graph gnew: graphcount = " << graphcount << endl; cout.flush();
	G.ndelta = 1;
	(*Islets).push_back(G);
	RemoveIsletGraph(&G);isletcount--;if(verbosemem)cout << "Removing isletgraph G: isletcount = " << isletcount << endl; cout.flush();
	IsletNumber = IsletNumberin;
      }
    }
    //*****************************************************************************//
    //                        Add cells to old islets                              //
    //*****************************************************************************//
    
    //If old islet
    else{
      //beta cells
      if (((AVopt == "rgc" && ctypein =="gfp") || (AVopt == "121113"&& ctypein == "11") ||(AVopt == "abd"&& ctypein == "b") ||
	   (AVopt == "321Area"&& ctypein == "2") || (AVopt == "rgcArea" && ctypein == "gfp") || (AVopt == "abdArea" && ctypein == "b")) &&
	  (celltype == "ab"||celltype =="bb" || celltype == "bd" ||celltype == "abd" )){
	u=add_vertex((*Islets).back().allGraph);
	InitializeVertex(&(*Islets).back().allGraph,u,xin,yin,num,"b",CellNumberin);
	num++;
	((*Islets).back()).nbeta++;
	IsletNumber = IsletNumberin;
      }
      //alpha cells
      else if (((AVopt == "rgc" && ctypein =="rfp") || (AVopt == "121113"&& ctypein == "12") ||(AVopt == "abd"&& ctypein == "a") ||
		(AVopt == "321Area"&& ctypein == "3") || (AVopt == "rgcArea" && ctypein == "rfp")|| (AVopt == "abdArea" && ctypein == "a") ) &&
	       (celltype == "aa"||celltype =="ab" || celltype == "ad" ||celltype == "abd" || celltype == "adm")){ 
      	u=add_vertex((*Islets).back().allGraph);
	InitializeVertex(&(*Islets).back().allGraph,u,xin,yin,num,"a",CellNumberin);
	num++;
	((*Islets).back()).nalpha++;
	IsletNumber = IsletNumberin;
      }
      //delta cells
      else if (((AVopt == "rgc" && ctypein =="cy5") || (AVopt == "121113"&& ctypein == "13") ||(AVopt == "abd"&& ctypein == "d") ||
		(AVopt == "321Area"&& ctypein == "1") || (AVopt == "rgcArea" && ctypein == "cy5") || (AVopt == "abdArea" && ctypein == "d")) &&
	       (celltype == "ad"||celltype =="bd" || celltype == "dd" ||celltype == "abd" || celltype == "adm" )){//delta cells
      	u=add_vertex((*Islets).back().allGraph);
	InitializeVertex(&(*Islets).back().allGraph,u,xin,yin,num,"d",CellNumberin);
	num++;
	((*Islets).back()).ndelta++;
	IsletNumber = IsletNumberin;
      }
    }
  }
  fin.close();

  if(verbose){
    cout << "Cells added are " << endl;
    for(int i = 0; i < num_vertices((*Islets).back().allGraph);i++)
      cout << (*Islets).back().allGraph[i].IsletCellNum << "\t" << (*Islets).back().allGraph[i].type << endl;
  }
  if(graphcount > 0 || isletcount > 0) cout << "Problem with islets and/or graphs! graphcount = " << graphcount << " and isletcount = " << isletcount << " in AddVertices " << endl; 
  if(verbosemem)cout << "*****Finished AddVertices" << endl;cout.flush();
  return;
}

//*********************************************************************************************//
//*********************************************************************************************//
//*********************************************************************************************//

void CopyComponentEdgesFromGraph(IsletGraph* IGpold, int c,Graph * Gpnew, string type){
  std::pair<vertex_iter, vertex_iter> vp1,vp2;
  edge_t eold,enew;
  Graph * Gpold;
  vector<vector<vertex_t> > *Cpold;
  bool verbose=0;

  //Set graph and component pointers
  if(type == "all"){
    Gpold = &(*IGpold).allGraph;
    Cpold = &(*IGpold).allComponents;
  }
  else if(type == "ad"){
    Gpold = &(*IGpold).adGraph;
    Cpold = &(*IGpold).adComponents;
  }
  else if(type == "b"){
    Gpold = &(*IGpold).bGraph;
    Cpold = &(*IGpold).bComponents;
  }
  
  //Add edges
  bool b = 0;
  for(vp1 = vertices(*Gpnew);vp1.first!=vp1.second;vp1.first++){//vertices of component
    if(verbose)cout << "Adding edges for vertex " << (*Gpnew)[*vp1.first].IsletCellNum << ":" << endl; cout.flush();
    vertex_t vp1old;
    for(int old1 = 0; old1 < (*Cpold)[c].size(); old1++){//finding vp1 in original graph component
      if((*Gpnew)[*vp1.first].IsletCellNum == (*Gpold)[(*Cpold)[c][old1]].IsletCellNum){
	vp1old = (*Cpold)[c][old1];
	old1 = (*Cpold)[c].size();
      }
    }
    if(verbose)cout << "  vp1old is " << (*Gpold)[vp1old].IsletCellNum << " with degree " << out_degree(vp1old, *Gpold) << endl;cout.flush();
    for(vp2 =vertices(*Gpnew);vp2.first!=vp2.second;vp2.first++){//vertices of component
      vertex_t vp2old;
      for(int old2 = 0; old2 < (*Cpold)[c].size() ;old2++){//vertices of original graph component
	if((*Gpnew)[*vp2.first].IsletCellNum == (*Gpold)[(*Cpold)[c][old2]].IsletCellNum){
	  vp2old = (*Cpold)[c][old2];
	  old2 = (*Cpold)[c].size();
	}
      }
      if(verbose) cout << "  vp2old is " << (*Gpold)[vp2old].IsletCellNum << " with degree " << out_degree(vp2old, *Gpold) << endl;cout.flush();
      if(vp1old!=vp2old){
	//Check if edge already exists in comp
	tie(enew,b)=edge(*vp1.first,*vp2.first,(*Gpnew));
	if(!b){
	  if(verbose)cout << "edge does not exist in new graph" << endl;cout.flush();
	  //Check if edge is in Gpold
	  tie(eold,b) = edge(vp1old,vp2old,(*Gpold));
	  if(b){
	    if(verbose)cout << " but edge exists in old graph" << endl; cout.flush();
	    tie(enew,b)=add_edge(*vp1.first,*vp2.first,(*Gpnew));
	    InitializeEdge(Gpnew,enew,(*Gpold)[eold].distance);
	    if(verbose)cout <<  "  added edge" << endl;cout.flush();
	  }
	}
      }
    }
  }

   return;
 }

//*********************************************************************************************//
//*********************************************************************************************//
//*********************************************************************************************//

void CopyComponentVerticesFromGraph(IsletGraph* IGpold, int c,Graph * Gpnew, string type, string votherType){
  vertex_t u;
  Graph * Gpold;
  vector<vector<vertex_t> > *Cpold;
  bool verbose=0;

  //Set graph and component pointers
  if(type == "all"){
    Gpold = &(*IGpold).allGraph;
    Cpold = &(*IGpold).allComponents;
  }
  else if(type == "ad"){
    Gpold = &(*IGpold).adGraph;
    Cpold = &(*IGpold).adComponents;
  }
  else if(type == "b"){
    Gpold = &(*IGpold).bGraph;
    Cpold = &(*IGpold).bComponents;
  }
  
  //Add vertices
  int num = 0;
  for(int j = 0; j <(*Cpold)[c].size();j++){
    u=add_vertex(*Gpnew);
    if(votherType == "")InitializeVertexFromVertex(Gpold, Gpnew, (*Cpold)[c][j] , u, num);
    else if(votherType == "all")InitializeVertexFromVertex(Gpold, Gpnew, (*Cpold)[c][j] , u, num, (*Gpold)[(*Cpold)[c][j]].vother);
    num++;
    (*Gpnew)[u].CompNum = c;
    if(verbose)cout << "Adding cell " << (*Gpold)[(*Cpold)[c][j]].IsletCellNum << " to component" << endl;
  }
  return;
}

//*********************************************************************************************//
//*********************************************************************************************//
//*********************************************************************************************//

void CopyGraph(Graph* Gpold,Graph* Gpnew){
  vertex_t u,v;
  std::pair<vertex_iter, vertex_iter> vp;
  std::pair<edge_iter, edge_iter> ep;
  edge_t enew;
  bool b = 0;
  int num = 0;
  bool verbose=0;

  //Add vertices
  for(vp = vertices(*Gpold);vp.first!=vp.second;vp.first++){
    u=add_vertex((*Gpnew));
    InitializeVertexFromVertex(Gpold, Gpnew, *vp.first, u, num);
    num++;
    if(verbose)cout << "Adding vertex " << (*Gpnew)[u].IsletCellNum << endl;
  }
  
  //Add edges
  for(ep = edges(*Gpold);ep.first!=ep.second;ep.first++){
    u = source(*ep.first,*Gpold);
    v = target(*ep.first,*Gpold);
    if(verbose)cout << "Checking edge between " << (*Gpold)[u].IsletCellNum << " and " << (*Gpold)[v].IsletCellNum << endl;cout.flush();
    tie(enew,b) = edge(u,v,(*Gpnew));
    if(!b){
      if(verbose)cout << "Adding edge between " << (*Gpnew)[u].IsletCellNum << " and " << (*Gpnew)[v].IsletCellNum;
      tie(enew,b)=add_edge(u,v,(*Gpnew));
      InitializeEdge(Gpnew,enew,(*Gpold)[*ep.first].distance);
      if(verbose)cout << " finished" << endl;cout.flush();
    }
  }
  
  return;
}

//*********************************************************************************************//
//*********************************************************************************************//
//*********************************************************************************************//


void CopyGraphEdges(Graph* Gpold, Graph * Gpnew, string type){
  edge_t eold,enew;
  bool b = 0;
  bool verbose=0;

  for(int i = 0; i < num_vertices(*Gpnew);i++){
    if(verbose)cout << "For i = " << i << " and cell number " << (*Gpnew)[i].IsletCellNum << " with type " << (*Gpnew)[i].type << endl;
    //if((*Gpnew)[i].CellNum < 99999){
      for(int j = i+1; j < num_vertices(*Gpnew);j++){
	if(verbose)cout << "  and j = " << j << " and cell number " << (*Gpnew)[j]. IsletCellNum << " with type " << (*Gpnew)[i].type <<  endl;
	//if((*Gpnew)[j].CellNum < 99999 &&
	if(
	   (type == "abd" ||
	    (type == "adb" && (((*Gpnew)[i].type =="a" && (*Gpnew)[j].type == "b") ||
			       ((*Gpnew)[i].type =="d" && (*Gpnew)[j].type == "b") ||
			     ((*Gpnew)[i].type =="b" && (*Gpnew)[j].type == "a") ||
			       ((*Gpnew)[i].type =="b" && (*Gpnew)[j].type == "d"))) ||
	    (type == "bb" && ((*Gpnew)[i].type =="b" && (*Gpnew)[j].type == "b")) ||
	    (type == "adm" &&(((*Gpnew)[i].type =="a" && (*Gpnew)[j].type == "a") ||
			      ((*Gpnew)[i].type =="d" && (*Gpnew)[j].type == "d") ||
			      ((*Gpnew)[i].type =="d" && (*Gpnew)[j].type == "a") ||
			      ((*Gpnew)[i].type =="a" && (*Gpnew)[j].type == "d"))))){
	  tie(eold,b)=edge((*Gpnew)[i].vother,(*Gpnew)[j].vother,*Gpold);
	  if(b){
	    if(verbose)cout << "    An edge exists between cells (from vother) " << (*Gpold)[(*Gpnew)[i].vother].IsletCellNum << " and " << (*Gpold)[(*Gpnew)[j].vother].IsletCellNum << endl;
	    tie(enew,b)=edge(i,j,(*Gpnew));
	    if(!b){
	      tie(enew,b)=add_edge(i,j,(*Gpnew));
	      InitializeEdge(Gpnew,enew,(*Gpold)[eold].distance);
	    }
	  }
	}
      }
  //}
    // else{//For IPs: add parent edges
    //   if(verbose)cout << "For cell " << (*Gpold)[i].IsletCellNum << " adding edges between " << endl;
    //   tie(enew,b)=edge((*Gpold)[i].IPparents[0],(*Gpold)[i].IPparents[1],(*Gpnew));
    //   if(!b){
    // 	if(verbose)cout << (*Gpold)[(*Gpold)[i].IPparents[0]].IsletCellNum << " and " << (*Gpold)[(*Gpold)[i].IPparents[1]].IsletCellNum << endl;
    // 	tie(enew,b)=add_edge((*Gpold)[i].IPparents[0],(*Gpold)[i].IPparents[1],(*Gpnew));
    // 	double dist = Distancexy((*Gpnew)[(*Gpold)[i].IPparents[0]].x,(*Gpnew)[(*Gpold)[i].IPparents[0]].y,(*Gpnew)[(*Gpold)[i].IPparents[1]].x,(*Gpnew)[(*Gpold)[i].IPparents[1]].y);
    // 	InitializeEdge(Gpnew,enew,dist,dist,0);
    //   }
    //   tie(enew,b)=edge((*Gpold)[i].IPparents[2],(*Gpold)[i].IPparents[3],(*Gpnew));
    //   if(!b){
    // 	if(verbose)cout << (*Gpold)[(*Gpold)[i].IPparents[0]].IsletCellNum << " and " << (*Gpold)[(*Gpold)[i].IPparents[1]].IsletCellNum << endl;
    // 	tie(enew,b)=add_edge((*Gpold)[i].IPparents[2],(*Gpold)[i].IPparents[3],(*Gpnew));
    // 	double dist = Distancexy((*Gpnew)[(*Gpold)[i].IPparents[2]].x,(*Gpnew)[(*Gpold)[i].IPparents[2]].y,(*Gpnew)[(*Gpold)[i].IPparents[3]].x,(*Gpnew)[(*Gpold)[i].IPparents[3]].y);
    // 	InitializeEdge(Gpnew,enew,dist,dist,0);
    //   }
    // }
  }
  
  return;
}

//*********************************************************************************************//
//*********************************************************************************************//
//*********************************************************************************************//

void CopyGraphVertices(Graph* Gpold, Graph * Gpnew, string type, bool setvother){
  vertex_t u;
  std::pair<vertex_iter, vertex_iter> vp;

  int num = 0;
  for(vp = vertices(*Gpold);vp.first!=vp.second;vp.first++){
    if((*Gpold)[*vp.first].CellNum < 99999 && 
       (type == "abd" ||
	(type == "b" && (*Gpold)[*vp.first].type == "b") ||
	(type == "adm" && ((*Gpold)[*vp.first].type == "a" ||
			   (*Gpold)[*vp.first].type == "d")))){
      u=add_vertex((*Gpnew));
      InitializeVertexFromVertex(Gpold, Gpnew, *vp.first, u, num);
      (*Gpnew)[u].vother = *vp.first;
      if(setvother)(*Gpold)[*vp.first].vother = num;
      num++;
    }
  }
  
  return;
}

//*********************************************************************************************//
//*********************************************************************************************//
//*********************************************************************************************//

void CopyIsletGraph(IsletGraph* Iold, IsletGraph* Inew){
  bool verbose=0;
  CopyGraph(&(*Iold).allGraph,&(*Inew).allGraph);
  if(verbose)cout << "There are " << num_vertices((*Inew).allGraph) << " vertices and " << num_edges((*Inew).allGraph) << " edges in the new graph" << endl;
  (*Inew).IsletNum = (*Iold).IsletNum;
  (*Inew).nalpha=(*Iold).nalpha;
  (*Inew).nbeta=(*Iold).nbeta;
  (*Inew).ndelta=(*Iold).ndelta;
  (*Inew).IsletArea=(*Iold).IsletArea; 
  (*Inew).admthresh=(*Iold).admthresh;
  (*Inew).bbthresh=(*Iold).bbthresh;
  (*Inew).adbthresh=(*Iold).adbthresh;
  for(int i = 0; i < (*Iold).bComponents.size();i++){
    vector<vertex_t> temp;
    for(int j = 0; j < (*Iold).bComponents[i].size();j++)
      temp.push_back((*Iold).bComponents[i][j]);
    (*Inew).bComponents.push_back(temp);
    RemoveVector(&temp);
  }
  for(int i = 0; i < (*Iold).bExteriorPaths.size();i++){
    vector<vertex_t> temp;
    for(int j = 0; j < (*Iold).bExteriorPaths[i].size();j++)
      temp.push_back((*Iold).bExteriorPaths[i][j]);
    (*Inew).bExteriorPaths.push_back(temp);
    RemoveVector(&temp);
  }
  
  for(int i = 0; i < (*Iold).adExteriorPaths.size();i++){
    vector<vertex_t> temp;
    for(int j = 0; j < (*Iold).adExteriorPaths[i].size();j++)
      temp.push_back((*Iold).adExteriorPaths[i][j]);
    (*Inew).adExteriorPaths.push_back(temp);
    RemoveVector(&temp);
  }
  for(int i = 0; i < (*Iold).adComponents.size();i++){
    vector<vertex_t> temp;
    for(int j = 0; j < (*Iold).adComponents[i].size();j++)
      temp.push_back((*Iold).adComponents[i][j]);
    (*Inew).adComponents.push_back(temp);
    RemoveVector(&temp);
  }
  return;
}

//*********************************************************************************************//
//*********************************************************************************************//
//*********************************************************************************************//


void Create_adAnd_bGraphs(IsletGraph * I){
  vertex_t u,v;
  edge_t e;
  std::pair<edge_iter, edge_iter> ep;
  std::pair<vertex_iter, vertex_iter> vp1;
  bool verbose=0;

  //Just in case: clear out old ad and b Graphs
  RemoveGraph(&(*I).bGraph);
  RemoveGraph(&(*I).adGraph);
  
  //add vertices and set vother
  bool setvother = 1;
  CopyGraphVertices(&(*I).allGraph,&(*I).bGraph,"b",setvother);
  CopyGraphVertices(&(*I).allGraph,&(*I).adGraph,"adm",setvother);
  if(verbose)cout << "Added " << num_vertices((*I).bGraph) << " betas and " << num_vertices((*I).adGraph) << " alphas and deltas" << endl;
  //add edges
  CopyGraphEdges(&(*I).allGraph,&(*I).bGraph, "bb");
  CopyGraphEdges(&(*I).allGraph,&(*I).adGraph, "adm");

  if(verbose)cout << "Finished add edges" << endl;

  return;
}

//*********************************************************************************************//
//*********************************************************************************************//
//*********************************************************************************************//

void InitializeEdge(Graph * G, edge_t e, double distance){
  (*G)[e].distance = distance;
  
  return;
}

// //*********************************************************************************************//
// //*********************************************************************************************//
// //*********************************************************************************************//

void InitializeIslet(IsletGraph * I, Graph* gnew, int IsletNumberin, double IsletAreain){
  
  (*I).allGraph = (*gnew);
  (*I).IsletNum = IsletNumberin;
  (*I).IsletArea = IsletAreain;
  (*I).nalpha = 0;
  (*I).nbeta = 0;
  (*I).ndelta = 0;
  (*I).admthresh = 0;
  (*I).adbthresh = 0;
  (*I).bbthresh = 0;
  return;

}

// //*********************************************************************************************//
// //*********************************************************************************************//
// //*********************************************************************************************//

void InitializeVertex(Graph * g, vertex_t u,double  xin,double yin, int num, string type, int CellNumberin, int vother){
  //default: vother = 0

  (*g)[u].x = xin;
  (*g)[u].y = yin;
  (*g)[u].CellNum = num;
  (*g)[u].type = type;
  (*g)[u].IsletCellNum = CellNumberin;
  // (*g)[u].ncellcomp = 0;
  (*g)[u].CompNum = 0;
  //(*g)[u].componentlabel = 0;
  (*g)[u].vother = vother;
 
  return;
}

//*****************************************************************************************
//*****************************************************************************************

void InitializeVertexFromVertex(Graph * gold, Graph* gnew, vertex_t uold , vertex_t unew, int num, int vother){
  bool verbose=0;

  (*gnew)[unew].x =(*gold)[uold].x ;
  (*gnew)[unew].y = (*gold)[uold].y ;
  (*gnew)[unew].CellNum =(*gold)[uold].CellNum ;
  (*gnew)[unew].IsletCellNum =(*gold)[uold].IsletCellNum ;
  //(*gnew)[unew].CompNum =(*gold)[uold].CompNum ;
  (*gnew)[unew].type =(*gold)[uold].type ;
  if((*gnew)[unew].IsletCellNum < 100000){
    (*gnew)[unew].CellNum = num;
    //(*gnew)[unew].ncellcomp = (*gold)[uold].ncellcomp;
    
    //(*gnew)[unew].componentlabel =(*gold)[uold].componentlabel ;
    
  }
  
  if(verbose){
    cout << "Gpnew[u] contains: " << endl;
    cout << "x:" << (*gnew)[unew].x << endl;
    cout << "y:" << (*gnew)[unew].y << endl;
    cout << "CellNum:" << (*gnew)[unew].CellNum << endl;
    cout << "IsletCellNum:" << (*gnew)[unew].IsletCellNum << endl;
    //cout << "CompNum:" << (*gnew)[unew].CompNum << endl;
    cout << "Type:" << (*gnew)[unew].type << endl;
  }
  if(verbose)cout << "Setting vother for cell " << (*gnew)[unew].IsletCellNum << " to " << (*gold)[(*gnew)[unew].vother].IsletCellNum << endl;
  return;
}


//*********************************************************************************************//
//*********************************************************************************************//
//*********************************************************************************************//

// bool AddEdgesSaturated(Graph* Gp, string celltype, double shadowthresh, int label){
//   //Calculated differently from AddEdgesThresh and AddEdges_gr
  
//   std::pair<vertex_iter, vertex_iter> vp1,vp2;//vp,vtemp2;
//   edge_t e;
  
//   bool verbose=0,verbosemem=0;
//   int vectorcount = 0;
//   if(verbose || verbosemem)cout <<  "*****Inside AddEdgesSaturated:CellType is " << celltype << endl;
//   for (vp1 = vertices(*Gp);vp1.first != vp1.second; vp1.first++){
//     if(celltype=="abd" ||
//        (celltype == "aa" && (*Gp)[*vp1.first].type == "a") ||
//        (celltype == "bb" && (*Gp)[*vp1.first].type == "b") ||
//        (celltype == "dd" && (*Gp)[*vp1.first].type == "d") || 
//        (celltype == "ab" && ((*Gp)[*vp1.first].type == "a" ||(*Gp)[*vp1.first].type == "b")) ||
//        (celltype == "ad" && ((*Gp)[*vp1.first].type == "a" ||(*Gp)[*vp1.first].type == "d")) ||
//        (celltype == "bd" && ((*Gp)[*vp1.first].type == "b" ||(*Gp)[*vp1.first].type == "d")) ||
//        (celltype == "adm"&& ((*Gp)[*vp1.first].type == "a" ||(*Gp)[*vp1.first].type == "d")) || 
//        (celltype == "adb" )){
//       //find distance to all vertices and sort
//       vector<VertexDistance> Dists;vectorcount++;if(verbosemem)cout << "Adding vector Dists: vectorcount = " << vectorcount << endl;cout.flush();
//       for(vp2 = vertices(*Gp);vp2.first!=vp2.second;vp2.first++){ //Include all vertices in list of shadows
// 	if(*vp1.first!=*vp2.first){
// 	  double dist = Distancexy((*Gp)[*vp1.first].x,(*Gp)[*vp1.first].y,(*Gp)[*vp2.first].x,(*Gp)[*vp2.first].y);
// 	  VertexDistance temp;
// 	  temp.vertex=*vp2.first;
// 	  temp.distance=dist;
// 	  Dists.push_back(temp);
// 	}
//       }
//       sort(Dists.begin(),Dists.end(),compareByDistance);
      
//       //Remove intervals and add edges to graph
//       bool intervalfull = 0;
//       vector<vector<double> > interval; vectorcount++;if(verbosemem)cout << "Adding vector interval: vectorcount = " << vectorcount << endl;cout.flush();
//       vector<double> temp(2,0);vectorcount++;if(verbosemem)cout << "Adding vector temp: vectorcount = " << vectorcount << endl;cout.flush();
//       temp[0]=0;temp[1]=2*PI;
//       interval.push_back(temp);
//       int dit = 0;
//       while(!intervalfull && dit < Dists.size()){
// 	//Determine if dit is able to be added
// 	//Calculate Angles
// 	float alpha2 = atan(shadowthresh/Dists[dit].distance);
// 	float beta2 = 0.0;
// 	double x1 = (*Gp)[*vp1.first].x,y1 = (*Gp)[*vp1.first].y;
// 	double x2 = (*Gp)[Dists[dit].vertex].x,y2 = (*Gp)[Dists[dit].vertex].y;
// 	if(x2>=x1)beta2 = atan((y2-y1)/(x2-x1));//Q1 and 4
// 	else if(y2>y1)beta2 = atan((y2-y1)/(x2-x1))+ PI;//Q2
// 	else beta2 = atan((y2-y1)/(x2-x1))- PI;//Q3
// 	beta2 = mod(beta2,2*PI);
// 	vector<double> beta2shad(2,0);vectorcount++;if(verbosemem)cout << "Adding vector beta2shad: vectorcount = " << vectorcount << endl;cout.flush();
// 	beta2shad[0] = mod(beta2-alpha2,2*PI);beta2shad[1] = mod(beta2+alpha2,2*PI);
// 	if(verbose)cout<< "Checking " << (*Gp)[*vp1.first].IsletCellNum << " with " << (*Gp)[Dists[dit].vertex].IsletCellNum << "\n";
// 	if(verbose)cout << "v2's shadow is (" << beta2shad[0] << "," <<  beta2shad[1]  << ")" << endl;
	
// 	bool betagood = 0;
// 	if(beta2shad[0]<beta2shad[1]){//new segment doesn't wrap
// 	  //loop thru segments to see where beta2shad[0] falls
// 	  for(int sn = 0;sn< interval.size(); sn++){
// 	    if(beta2shad[0]> interval[sn][0] && beta2shad[0] < interval[sn][1]){// beta2shad0 is located in this segment
// 	      if(beta2shad[1] < interval[sn][1]){// beta2shad1 is located in this segment
// 		betagood=1;
// 		temp[0] = beta2shad[1];
// 		temp[1] = interval[sn][1];
// 		interval[sn][1] = beta2shad[0];
// 		interval.insert(interval.begin()+sn+1,temp);
// 		sn = interval.size();
// 		if(verbose)cout << "Case 1: not checked" << endl;
// 	      }
// 	      else{//beta2shad1 is located in another segment
// 		interval[sn][1] = beta2shad[0];
// 		while(sn < interval.size()-1 && beta2shad[1]>interval[sn+1][1]){
// 		  interval.erase(interval.begin()+sn+1);
// 		}
// 		if(sn < interval.size()-1 && beta2shad[1]>interval[sn+1][0]){
// 		  interval[sn+1][0] = beta2shad[1];
// 		}
// 		if(verbose)cout << "Case 2: not checked" << endl;
// 	      }
// 	      sn = interval.size();
// 	    }
	    
// 	    else if((interval.size() > 1 && sn < interval.size()-1)&& (beta2shad[0]> interval[sn][1] && beta2shad[0] < interval[sn+1][0])){//beta2shad0 is located between two segments
// 	      while(sn < interval.size()-1 && beta2shad[1] > interval[sn+1][1]){//remove segments that are in shadow
// 		interval.erase(interval.begin()+sn+1);
// 	      }
// 	      if(sn < interval.size()-1 && beta2shad[1] > interval[sn+1][0]){//update segment that contains beta2shad1
// 		// cout << "check here" << endl;
// 		interval[sn+1][0] = beta2shad[1];
// 	      }
// 	      if(verbose)cout << "Case 3: not checked" << endl;
// 	    }
	    
// 	    else if(beta2shad[0] < interval[sn][0] && sn == 0){//beta2shad0 is located before first segment
// 	      while(beta2shad[1] > interval[sn][1] && sn < interval.size()){//remove segments that are in shadow
// 		interval.erase(interval.begin()+sn);
// 	      }
// 	      if(beta2shad[1] > interval[sn][0] && beta2shad[1] < interval[sn][1]){
// 		interval[sn][0] = beta2shad[1];
// 	      }
// 	      if(verbose)cout << "Case 4: not checked" << endl;
// 	    }
// 	  }
// 	}
// 	else{//new segment wraps
// 	  if(beta2shad[0] < interval[interval.size()-1][1] && beta2shad[0] > interval[interval.size()-1][0] && interval[interval.size()-1][1] == 2*PI && interval[0][0]==0 && beta2shad[1] < interval[0][1]){//good wrap
// 	    betagood=1;
// 	    interval[interval.size()-1][1]=beta2shad[0];
// 	    interval[0][0]=beta2shad[1];
// 	    if(verbose)cout << "Case 5: not checked" << endl;
// 	  }
// 	  else{
// 	    for(int sn = 0;sn < interval.size(); sn++){
// 	      while(beta2shad[0]< interval[sn][0] && sn < interval.size()){//remove segments that are in shadow
// 		interval.erase(interval.begin()+sn);
// 		if(verbose)cout << "Case 6: not checked" << endl;
// 	      }
// 	      while(beta2shad[1] > interval[sn][1] && sn < interval.size()){
// 		interval.erase(interval.begin()+sn);
// 		if(verbose)cout << "Case 7: not checked" << endl;
// 	      }
// 	      if( sn < interval.size() && beta2shad[1]> interval[sn][0] && beta2shad[1] < interval[sn][1]){// beta2shad1 is located in this segment
// 		interval[sn][0] = beta2shad[1];
// 		if(verbose)cout << "Case 8: not checked" << endl;
// 	      }
// 	      if( sn < interval.size() && beta2shad[0]>interval[sn][0] && beta2shad[0] < interval[sn][1]){//beta2shad0 is located in this segment
// 		interval[sn][1] = beta2shad[0];
// 		if(verbose)cout << "Case 9: not checked" << endl;
// 	      }
// 	    }
// 	  }
// 	}
// 	//check intervals 
// 	for(int sn = 0; sn < interval.size();sn++){
// 	  if(interval[sn][0]>interval[sn][1])cout << "Problem with intervals"  << endl;
// 	}
// 	if(betagood){
// 	  if(verbose)cout << "Beta is good!" << endl;
// 	  bool b;
// 	  boost::tie(e,b) = edge(*vp1.first,Dists[dit].vertex,(*Gp));
// 	  if(!b){ 
// 	    if(celltype=="abd"  ||
//  	       (celltype =="aa" && (*Gp)[*vp1.first].type == "a" && (*Gp)[Dists[dit].vertex].type == "a") ||
//  	       (celltype =="ab" && ((*Gp)[*vp1.first].type == "a" ||(*Gp)[*vp1.first].type == "b") &&
//  		((*Gp)[Dists[dit].vertex].type == "a" ||(*Gp)[Dists[dit].vertex].type == "b") && 
//  		((*Gp)[*vp1.first].type !=(*Gp)[Dists[dit].vertex].type))||
//  	       (celltype =="ad" && ((*Gp)[*vp1.first].type == "a" ||(*Gp)[*vp1.first].type == "d") &&
//  		((*Gp)[Dists[dit].vertex].type == "a" ||(*Gp)[Dists[dit].vertex].type == "d") && 
//  		((*Gp)[*vp1.first].type !=(*Gp)[Dists[dit].vertex].type))||
//  	       (celltype =="bb" && (*Gp)[*vp1.first].type == "b" && (*Gp)[Dists[dit].vertex].type == "b") ||
//  	       (celltype =="bd" && ((*Gp)[*vp1.first].type == "b" ||(*Gp)[*vp1.first].type == "d") &&
//  		((*Gp)[Dists[dit].vertex].type == "b" ||(*Gp)[Dists[dit].vertex].type == "d") && 
//  		((*Gp)[*vp1.first].type !=(*Gp)[Dists[dit].vertex].type))||
//  	       (celltype == "dd" && (*Gp)[*vp1.first].type == "d" && (*Gp)[Dists[dit].vertex].type == "d") || 
//  	       (celltype == "adm"&& ((*Gp)[*vp1.first].type == "a" ||(*Gp)[*vp1.first].type == "d") &&
//  		((*Gp)[Dists[dit].vertex].type == "a" ||(*Gp)[Dists[dit].vertex].type == "d"))|| 
//  	       (celltype == "adb"&& (((*Gp)[*vp1.first].type == "a" &&(*Gp)[Dists[dit].vertex].type == "d") ||
// 				     ((*Gp)[*vp1.first].type == "d" &&(*Gp)[Dists[dit].vertex].type == "a") ||
// 				     ((*Gp)[*vp1.first].type == "a" &&(*Gp)[Dists[dit].vertex].type == "b") ||
// 				     ((*Gp)[*vp1.first].type == "b" &&(*Gp)[Dists[dit].vertex].type == "a")) )){
// 	      if(verbosemem)cout << "Adding edge between " << (*Gp)[*vp1.first].IsletCellNum << " of type " <<  (*Gp)[*vp1.first].type << " with cell " << (*Gp)[Dists[dit].vertex].IsletCellNum << " of type " << (*Gp)[Dists[dit].vertex].type << endl;
// 	      boost::tie(e,b) = add_edge(*vp1.first,Dists[dit].vertex,(*Gp));
// 	      InitializeEdge(Gp,e,Dists[dit].distance,Dists[dit].distance,label);
// 	    }
// 	  }
	  
// 	}
// 	dit++;
// 	if(verbose){
// 	  cout << "The new intervals are: " ;
// 	  for(int i = 0; i < interval.size();i++){
// 	    cout << "[" << interval[i][0] << "," << interval[i][1] << "] ";
// 	  }
// 	  cout << endl;
// 	}
// 	if(interval.size()==0)
// 	  intervalfull=1;
// 	RemoveVector(&beta2shad);vectorcount--;if(verbosemem)cout << "Removing vector beta2shad: vectorcount = " << vectorcount << endl;cout.flush();

//       }
//       RemoveVector(&temp);vectorcount--;if(verbosemem)cout << "Removing vector temp: vectorcount = " << vectorcount << endl;cout.flush();
//       RemoveVector(&Dists);vectorcount--;if(verbosemem)cout << "Removing vector Dists: vectorcount = " << vectorcount << endl;cout.flush();
//       RemoveVector2d(&interval);vectorcount--;if(verbosemem)cout << "Removing vector interval: vectorcount = " << vectorcount << endl;cout.flush();
//     }
//   }
//   if(verbosemem)cout << "*****Finished AddEdgesSaturated"<<endl;
//   if(vectorcount > 0 ){
//     cout << "Problem with memory in AddEdgesSaturated !!! Too many vectors!!!" << endl;
//     cout << "vector count = " << vectorcount << endl;
//     return 1;
//   }
//   return 0;
// }

// //*********************************************************************************************//
// //*********************************************************************************************//
// //*********************************************************************************************//

// void AddIPEdges(Graph* G, string filename){
//   bool verbose=0;
//   std::pair<vertex_iter, vertex_iter> vp;
//   if(verbose)cout << "Reading in " << filename << endl;
//   ifstream fin(filename.c_str());
  
 
//   while(!fin.eof()){
//     int n1 = 0,n2 = 0,in1=0,in2=0; //n1 and 2 are the cellnums; in1 and 2 are the isletcellnumbers
//     double dist = 0;
    
//     fin >>n1 >> n2 >> in1 >> in2 >> dist;
//     if(fin.eof())break;
    
//     //Find v1 and v2 in G
//     int numInG = 0;
//     vertex_t v1,v2;
//     for(vp=vertices(*G);vp.first!=vp.second;vp.first++){
//       if((*G)[*vp.first].CellNum == n1){
// 	v1 = *vp.first;
// 	numInG++;
//       }
//       if((*G)[*vp.first].CellNum == n2){
// 	v2 = *vp.first;
// 	numInG++;
//       }
//       if(numInG == 2)break;
//     }
    
//     edge_t e;bool b;
//     tie(e,b) = edge(v1,v2, *G);
//     if(!b)tie(e,b) = add_edge(v1,v2,*G);
//     InitializeEdge(G,e,dist,dist,0);
//   }
  
//   return;
// }   


// //*********************************************************************************************//
// //*********************************************************************************************//
// //*********************************************************************************************//

// void AddIPVertices(Graph* G, string filename){
//   bool verbose=0,verbosemem = 0;
//   int vectorcount = 0;
//   vertex_t u;
//   ifstream fin(filename.c_str());
//   if(verbosemem)cout << "***** Inside AddIPVertices" << endl;
//   if(verbose)cout << "Reading in " << filename << endl;
//   while(!fin.eof()){
//     int CellNumberin = 0,IsletCellNumberin = 0,CompNumberin = 0;
//     double xin = 0,yin = 0;
//     vector<int> IPparents(4,0);vectorcount++;if(verbosemem)cout << "Adding vector IPparents:vectorcount = " << vectorcount << endl;cout.flush();
    
//     fin  >> CompNumberin >> CellNumberin  >> xin >> yin >> IPparents[0] >> IPparents[1]>> IPparents[2]>> IPparents[3];
//     if(fin.eof()){
//       RemoveVector(&IPparents);vectorcount--;if(verbosemem)cout << "Removing vector IPparents:vectorcount = " << vectorcount << endl;cout.flush();
//       break;
//     }
//     u=add_vertex(*G);
//     (*G)[u].x = xin;
//     (*G)[u].y = yin;
//     (*G)[u].CellNum = CellNumberin;
//     (*G)[u].IsletCellNum = CellNumberin;
//     (*G)[u].CompNum = CompNumberin;
//     for(int i = 0; i < 4; i++){
//       (*G)[u].IPparents.push_back(IPparents[i]);
//     }
//     RemoveVector(&IPparents);vectorcount--;if(verbosemem)cout << "Removing vector IPparents:vectorcount = " << vectorcount << endl;cout.flush();
//   }
//   if(vectorcount!=0)cout << "Problem with AddIPVertices! vectorcount = " << vectorcount << endl; cout.flush();
//   if(verbosemem)cout << "*****Finished AddIPVertices" << endl;
//   return;
// }   

// //*********************************************************************************************//
// //*********************************************************************************************//
// //*********************************************************************************************//



// //*********************************************************************************************//
// //*********************************************************************************************//
// //*********************************************************************************************//



// void CopyComponentEdgesFromGraphWithIPs(Graph* Gpold, int c,Graph * Gpnew){
//   bool verbose=0;
//   std::pair<vertex_iter, vertex_iter> vpn,vpo;
//   std::pair<out_edge_iter, out_edge_iter> oep;
//   edge_t enew;
//   bool b = 0;

//   //PrintoutGnuplotWithIPs(Gpold, "temp.oldcopyCompEdges");
//   if(verbose)cout << "looking for component " << c << endl;cout.flush();
//   for(vpo = vertices(*Gpold);vpo.first!=vpo.second;vpo.first++){
    
//     if((*Gpold)[*vpo.first].CompNum == c && out_degree(*vpo.first, *Gpold) > 0){
//       if(verbose)cout << "Checking cell " << (*Gpold)[*vpo.first].IsletCellNum << " which has component " << (*Gpold)[*vpo.first].CompNum << " and degree " << out_degree(*vpo.first, *Gpold) << endl;cout.flush();
//       vertex_t v1n;
//       //Find source in Gpnew
//       for(vpn =vertices(*Gpnew);vpn.first!=vpn.second;vpn.first++)
// 	if((*Gpold)[*vpo.first].IsletCellNum == (*Gpnew)[*vpn.first].IsletCellNum)
// 	  v1n = *vpn.first;
//       for(oep = out_edges(*vpo.first,*Gpold);oep.first!=oep.second;oep.first++){
// 	vertex_t vtarget = target(*oep.first,*Gpold);
// 	if(verbose)cout << "Found edge between " << (*Gpnew)[v1n].IsletCellNum  << " and " << (*Gpold)[vtarget].IsletCellNum << endl;
// 	vertex_t v2n;
// 	//Find target in Gpnew
// 	for(vpn =vertices(*Gpnew);vpn.first!=vpn.second;vpn.first++)
// 	  if((*Gpold)[vtarget].IsletCellNum == (*Gpnew)[*vpn.first].IsletCellNum)
// 	    v2n = *vpn.first;
// 	//Check if edge exists in Gpnew
// 	tie(enew,b) = edge(v1n,v2n,*Gpnew);
// 	if(!b){
// 	  if(verbose)cout << "...Adding edge between " << (*Gpnew)[v1n].IsletCellNum << " and " << (*Gpnew)[v2n].IsletCellNum << endl;
// 	  tie(enew,b) = add_edge(v1n,v2n,*Gpnew);
// 	  InitializeEdge(Gpnew,enew,(*Gpold)[*oep.first].distance,(*Gpold)[*oep.first].edge_weight,(*Gpold)[*oep.first].label);
// 	}
//       }
//     }
//   }
//   //PrintoutGnuplot(Gpnew, "temp.copyCompEdges");
//    return;
//  }

// //*********************************************************************************************//
// //*********************************************************************************************//
// //*********************************************************************************************//



// void CopyComponentVerticesFromGraphWithIPs(Graph* Gpold, int c,Graph * Gpnew){
//   vertex_t u;
//   bool verbose=0;
  
//   // if(c==3)PrintoutGnuplotWithIPs(Gpold, "temp.comp3");
//   //if(c==4)PrintoutGnuplotWithIPs(Gpold, "temp.comp4");
//   //Add vertices
//   std::pair<vertex_iter, vertex_iter> vp,vp2;
//   int num = 0;
  
//   for(vp = vertices(*Gpold);vp.first!=vp.second;vp.first++){
//     if(verbose)cout << "Checking cell " << (*Gpold)[*vp.first].IsletCellNum << " with component " << (*Gpold)[*vp.first].CompNum << endl;cout.flush();
//    if((*Gpold)[*vp.first].CompNum == c){
//      u=add_vertex(*Gpnew);
//        InitializeVertexFromVertex(Gpold, Gpnew, *vp.first , u,num );
//        (*Gpnew)[u].vother = (*Gpold)[*vp.first].CellNum;
//        for(int i = 0; i <(*Gpold)[*vp.first].IPparents.size() ;i++)
//     	(*Gpnew)[u].IPparents.push_back((*Gpold)[*vp.first].IPparents[i]);
//       num++;
//       if(verbose)
//     	cout << "Adding cell " << (*Gpold)[*vp.first].IsletCellNum << " to component " << c << " with old IP parents " << endl ;
//       if(verbose && (*Gpold)[*vp.first].CellNum > 99999){
//     	for(int i = 0; i < 4; i++)cout << (*Gpold)[(*Gpold)[*vp.first].IPparents[i]].IsletCellNum << " ";
//     	cout << endl;cout.flush();
//       }
//      }
//   }
//   //fix IPparents
//   for(vp = vertices(*Gpnew);vp.first!=vp.second;vp.first++){
//     if((*Gpnew)[*vp.first].CellNum > 99999){
//       for(int i = 0; i < 4; i++){
//   	vertex_t oldIP = (*Gpnew)[*vp.first].IPparents[i];
//   	for(vp2 = vertices(*Gpnew);vp2.first!=vp2.second;vp2.first++){
//   	  if(verbose)cout << "Checking " <<(*Gpold)[oldIP].IsletCellNum << " with " <<  (*Gpnew)[*vp2.first].IsletCellNum << endl;
//   	  if((*Gpold)[oldIP].IsletCellNum == (*Gpnew)[*vp2.first].IsletCellNum){
	    
//   	    (*Gpnew)[*vp.first].IPparents[i] = (*Gpnew)[*vp2.first].CellNum;
//   	    break;
//   	  }
//   	}
//       }
//       if(verbose){
//   	cout << "updated parents to " ;
//   	for(int i = 0; i < 4; i++)cout << (*Gpnew)[(*Gpnew)[*vp.first].IPparents[i]].IsletCellNum << " ";
//   	cout << endl;cout.flush();
//       }
//     }
//   }
	  
//   return;
// }

// //*********************************************************************************************//
// //*********************************************************************************************//
// //*********************************************************************************************//

// void CopyGraphMinusIP(Graph* Gpold,Graph* Gpnew){
//   vertex_t u;
//   std::pair<vertex_iter, vertex_iter> vp;
//   edge_t enew;
//   std::pair<out_edge_iter, out_edge_iter> oep;
//   bool b = 0;
//   int num = 0;
//   bool verbose=0;

//   //Add vertices
//   for(vp = vertices(*Gpold);vp.first!=vp.second;vp.first++){
//     if((*Gpold)[*vp.first].CellNum < 99999){
//       u=add_vertex((*Gpnew));
//       InitializeVertexFromVertex(Gpold, Gpnew, *vp.first, u, num);
//       num++;
//     }
//   }
//   //Add edges
//   for(vp = vertices(*Gpold);vp.first!=vp.second;vp.first++){
//     if((*Gpold)[*vp.first].CellNum < 99999){
//       for(oep = out_edges(*vp.first,(*Gpold));oep.first!=oep.second;oep.first++){
// 	u = target(*oep.first,(*Gpold));
// 	if((*Gpold)[u].CellNum > 99999){
// 	  if(verbose)cout << "For cell " << (*Gpold)[u].IsletCellNum << " with IPparents " << (*Gpold)[(*Gpold)[u].IPparents[0]].IsletCellNum << " " << (*Gpold)[(*Gpold)[u].IPparents[1]].IsletCellNum << " " << (*Gpold)[(*Gpold)[u].IPparents[2]].IsletCellNum << " " << (*Gpold)[(*Gpold)[u].IPparents[3]].IsletCellNum << " " << "adding edges between ";
// 	  vertex_t otherend;
// 	  if(*vp.first == (*Gpold)[u].IPparents[0]) otherend =(*Gpold)[u].IPparents[1]; 
// 	  else if(*vp.first == (*Gpold)[u].IPparents[1]) otherend =(*Gpold)[u].IPparents[0]; 
// 	  else if(*vp.first == (*Gpold)[u].IPparents[2]) otherend =(*Gpold)[u].IPparents[3]; 
// 	  else if(*vp.first == (*Gpold)[u].IPparents[3]) otherend =(*Gpold)[u].IPparents[2];
// 	  if(verbose)cout <<  (*Gpold)[*vp.first].IsletCellNum;cout.flush();
// 	  if(verbose)cout << " and " << (*Gpold)[otherend].IsletCellNum << endl;cout.flush();
// 	  tie(enew,b) = edge(*vp.first,otherend,(*Gpnew));
// 	  if(!b){
// 	    if(verbose)cout << "Adding edge between " << (*Gpold)[*vp.first].IsletCellNum << " and " << (*Gpold)[otherend].IsletCellNum << endl;cout.flush();
// 	    tie(enew,b)=add_edge(*vp.first,otherend,(*Gpnew));
// 	    double dist = Distancexy((*Gpold)[*vp.first].x,(*Gpold)[*vp.first].y,(*Gpold)[otherend].x,(*Gpold)[otherend].y);
// 	    InitializeEdge(Gpnew,enew,dist,dist,0);
// 	  }
// 	}
// 	else{
// 	  tie(enew,b) = edge(*vp.first,u,(*Gpnew));
// 	  if(!b){
// 	    tie(enew,b)=add_edge(*vp.first,u,(*Gpnew));
// 	    InitializeEdge(Gpnew,enew,(*Gpold)[*oep.first].distance,(*Gpold)[*oep.first].edge_weight,(*Gpold)[*oep.first].label);
// 	  }
// 	}
//       }
//     }
//   }
//   return;
// }

// //*********************************************************************************************//
// //*********************************************************************************************//
// //*********************************************************************************************//






// void CopyOppositeComponents(IsletGraph* IGpold, Graph * Gpcopy,string type,int c1, int c2, Graph *Gnew ){
//   vertex_t u;
//   Graph * Gp1, *Gp2;
//   vector<vertex_t> *Cp1, *Cp2;
//   vector<int>  *Ep1, *Ep2;
//   bool verbose=0;
//   std::pair<vertex_iter, vertex_iter> vp1,vp2,vp;
//   std::pair<out_edge_iter, out_edge_iter> oep; 

//   if(type == "ad"){
//     Gp1 = &(*IGpold).adGraph;
//     Gp2 = &(*IGpold).bGraph;
//     Cp1 = &((*IGpold).adComponents[c1]);
//     Cp2 = &((*IGpold).bComponents[c2]);
//     Ep1 = &((*IGpold).adExteriorPaths[c1]);
//     Ep2 = &((*IGpold).bExteriorPaths[c2]);
//   }
//   else if(type == "b"){
//     Gp1 = &(*IGpold).bGraph;
//     Gp2 = &(*IGpold).adGraph;
//     Cp1 = &((*IGpold).bComponents[c1]);
//     Cp2 = &((*IGpold).adComponents[c2]);
//     Ep1 = &((*IGpold).bExteriorPaths[c1]);
//     Ep2 = &((*IGpold).adExteriorPaths[c2]);
//   }
  
//   //Add vertices
//   int num = 0;
//   for(int j = 0; j <(*Cp1).size();j++){
//     //check if j is in Exterior path
//     bool inpath = 0; 
//     for(int ep = 0; ep < (*Ep1).size();ep++){
//       if((*Gp1)[(*Cp1)[j]].CellNum == (*Ep1)[ep]){
// 	inpath = 1;
// 	break;
//       }
//     }
//     if(inpath){
//       u=add_vertex(*Gnew);
//       InitializeVertexFromVertex(Gp1, Gnew, (*Cp1)[j] , u, num);
//       if(verbose)cout << "Adding cell " << (*Gp1)[(*Cp1)[j]].IsletCellNum << " to component" << endl;
//     }    
//   }
//   for(int j = 0; j <(*Cp2).size();j++){
//     //check if j is in Exterior path
//     bool inpath = 0; 
//     for(int ep = 0; ep < (*Ep2).size();ep++){
//       if((*Gp2)[(*Cp2)[j]].CellNum == (*Ep2)[ep]){
// 	inpath = 1;
// 	break;
//       }
//     }
//     if(inpath){
//       u=add_vertex(*Gnew);
//       InitializeVertexFromVertex(Gp2, Gnew, (*Cp2)[j] , u, num);
//       if(verbose)cout << "Adding cell " << (*Gp2)[(*Cp2)[j]].IsletCellNum << " to component" << endl;
//     }    
//   }

//   //get edges from gpcopy
//   for(vp1 = vertices(*Gnew);vp1.first!=vp1.second;vp1.first++){
//     int gpcopy1=0;
//     for(int i = 0; i < num_vertices(*Gpcopy);i++){
//       if((*Gnew)[*vp1.first].IsletCellNum == (*Gpcopy)[i].IsletCellNum){
// 	gpcopy1 = i;
// 	break;
//       }
//     }
//     for(vp2 = vertices(*Gnew);vp2.first!=vp2.second;vp2.first++){
//       int gpcopy2=0;
//       for(int i = 0; i < num_vertices(*Gpcopy);i++){
// 	if((*Gnew)[*vp2.first].IsletCellNum == (*Gpcopy)[i].IsletCellNum){
// 	  gpcopy2 = i;
// 	  break;
// 	}
//       }
//       edge_t eold, enew; bool b;
//       tie(eold,b) = edge(gpcopy1,gpcopy2, (*Gpcopy));
//       cout << "Checking edge between " << (*Gpcopy)[gpcopy1].IsletCellNum << " and " << (*Gpcopy)[gpcopy2].IsletCellNum << endl;
//       if(b){
// 	cout << " an edge exists: adding edge between  " <<(*Gnew)[*vp1.first].IsletCellNum << " and " << (*Gnew)[*vp2.first].IsletCellNum <<  endl; 
// 	tie(enew,b) = edge(*vp1.first,*vp2.first,(*Gnew));
// 	if(!b){
// 	  cout << ": adding edge between  " <<(*Gnew)[*vp1.first].IsletCellNum << " and " << (*Gnew)[*vp2.first].IsletCellNum <<  endl; 
// 	  tie(enew,b) = add_edge(*vp1.first,*vp2.first,(*Gnew));
// 	  InitializeEdge(Gnew,enew, (*Gpcopy)[eold].distance,(*Gpcopy)[eold].edge_weight, (*Gpcopy)[eold].label);
// 	}
//       }
//     }
//   }

//   cout << "Outer Edges for gpcopy are" <<endl;
//   for(vp = vertices(*Gnew);vp.first!=vp.second;vp.first++){
//     cout << (*Gnew)[*vp.first].IsletCellNum << ": ";
//     for(oep = out_edges(*vp.first,(*Gnew));oep.first!=oep.second;oep.first++)
//       cout<< (*Gnew)[target(*oep.first,(*Gnew))].IsletCellNum << " ";
//     cout << endl;
//   }
//   return;
// }

// //*********************************************************************************************//
// //*********************************************************************************************//
// //*********************************************************************************************//



// void InitializeIP(Graph * gold, vertex_t uold,IsletGraph *I){
 
//   int Cellcount = 0;
//   bool verbosemem=0;
//   if(verbosemem)cout << "****Inside Initialize: Initializing cell " << (*gold)[uold].IsletCellNum ;
//   Cell IPnew;Cellcount++;if(verbosemem)cout << "Adding Cell IPnew: Cellcount = " << Cellcount << endl; cout.flush();
  
//   IPnew.x = (*gold)[uold].x;
//   IPnew.y =(*gold)[uold].y ;
//   IPnew.CellNum = (*gold)[uold].CellNum ;
//   IPnew.type =(*gold)[uold].type;
//   IPnew.IsletCellNum =(*gold)[uold].IsletCellNum ;
//   IPnew.ncellcomp =(*gold)[uold].ncellcomp ;
//   IPnew.CompNum =(*gold)[uold].CompNum ;
//   //IPnew.CompCellNum =(*gold)[uold].CompCellNum ;
//   IPnew.componentlabel =(*gold)[uold].componentlabel ;
//   IPnew.vother =(*gold)[uold].vother ;
//   for(int i = 0; i < (*gold)[uold].IPparents.size();i++)
//     IPnew.IPparents.push_back((*gold)[uold].IPparents[i]);
//   (*I).IPs.push_back(IPnew);
//   RemoveCell(&IPnew);Cellcount--;if(verbosemem)cout << "Removing Cell IPnew: Cellcount = " << Cellcount << endl; cout.flush();
//   cout << " finished with type "<<  (*I).IPs[(*I).IPs.size()-1].type << " with parents " << (*gold)[(*I).IPs[(*I).IPs.size()-1].IPparents[0]].IsletCellNum << " " << (*gold)[(*I).IPs[(*I).IPs.size()-1].IPparents[1]].IsletCellNum << " " << (*gold)[(*I).IPs[(*I).IPs.size()-1].IPparents[2]].IsletCellNum << " " << (*gold)[(*I).IPs[(*I).IPs.size()-1].IPparents[3]].IsletCellNum << " " <<   endl;
  
//   if(Cellcount!=0)cout << "Problem with InitializeIP! Cellcount = " << Cellcount << endl; cout.flush(); 
//   if(verbosemem)cout << "*****Finished InitializeIP"<< endl;cout.flush();
//   return;
// }

// //*********************************************************************************************//
// //*********************************************************************************************//
// //*********************************************************************************************//





// void RemoveInteriorVertices( Graph * G){
//   std::pair<vertex_iter, vertex_iter> vp;
//   bool verbose=0;
  
//   if(verbose){
//     cout << "Before removing cells" << endl;
//     for(int i = 0; i < num_vertices(*G); i++)
//       cout << (*G)[i].IsletCellNum << " ";
//     cout << endl;
//   }
//   for(int v = 0; v < num_vertices(*G);v++){
//     if((*G)[v].componentlabel == 1){
//       if(verbose)cout << "Cell " << (*G)[v].IsletCellNum << " is in the interior" <<  endl; 
//       vp = vertices(*G);
//       for(int i = 0; i < v; i++)vp.first++;
//       if(verbose)cout << "Removing cell " << (*G)[*vp.first].IsletCellNum << endl;
//       ClearVertex(G,*vp.first);
//       remove_vertex(*vp.first,*G);
//     }
//   }
//   if(verbose){
//     cout << "After removing cells" << endl;
//     for(int i = 0; i < num_vertices(*G); i++)
//       cout << (*G)[i].IsletCellNum << " ";
//     cout << endl;
//   }
//   return;
// }


// //*********************************************************************************************//
// //*********************************************************************************************//
// //*********************************************************************************************//


