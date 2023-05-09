#include "./GraphNN_VersionForPaper.h"
#define PI 3.14159

using namespace boost;

//*********************************************************************************************//
//*********************************************************************************************//
//*********************************************************************************************//

void PrintoutEdges(Graph *Gp, string fprefix,string directory ,bool append){
  std::pair<edge_iter, edge_iter> ep;
  bool verbose=0;

  //Printout edges
  ofstream gnupl;
  if(append)gnupl.open((directory+"/figures/" + fprefix+".edges").c_str(),ios::app);
  else gnupl.open((directory+"/figures/" + fprefix+".edges").c_str());
  
  if(verbose)cout << "There are " << num_edges((*Gp)) << " edges in Gp" << endl;
  //Print out lines on Graphs 
  
  for(ep = edges((*Gp));ep.first!=ep.second;ep.first++){
    if((*Gp)[source(*ep.first,(*Gp))].CellNum <(*Gp)[target(*ep.first,(*Gp))].CellNum){
	      vertex_t v1 = source(*ep.first,(*Gp));
	      vertex_t v2 = target(*ep.first,(*Gp)); 
	      if(verbose)cout << "printing out edge between " << (*Gp)[v1].IsletCellNum << " and " << (*Gp)[v2].IsletCellNum << endl;
	      gnupl << (*Gp)[v1].x << "\t" << (*Gp)[v1].y <<endl;
	      gnupl << (*Gp)[v2].x << "\t" << (*Gp)[v2].y <<endl;
	      gnupl << endl;
    }
  }

  gnupl.close();

  // ADDED BY MANU
  if(append)gnupl.open((directory+"/figures/" + fprefix+".edges.format2").c_str(),ios::app);
  else gnupl.open((directory+"/figures/" + fprefix+".edges.format2").c_str());
  
  //Print out lines on Graphs 
  
  for(ep = edges((*Gp));ep.first!=ep.second;ep.first++){
    if((*Gp)[source(*ep.first,(*Gp))].CellNum <(*Gp)[target(*ep.first,(*Gp))].CellNum){
	      vertex_t v1 = source(*ep.first,(*Gp));
	      vertex_t v2 = target(*ep.first,(*Gp)); 
	      if(verbose)cout << "printing out edge between " << (*Gp)[v1].IsletCellNum << " and " << (*Gp)[v2].IsletCellNum << endl;
	      gnupl << (*Gp)[v1].IsletCellNum << "," << (*Gp)[v2].IsletCellNum <<endl;
    }
  }

  gnupl.close();

  return;
}

//*********************************************************************************************//
//*********************************************************************************************//
//*********************************************************************************************//

void PrintoutExteriorPaths(IsletGraph *I, string fprefix,string directory, string type ){
  bool verbose=0;
  Graph *Gp;
  vector<vector<vertex_t> > * Pp;
  if(verbose)cout << "Inside PrintoutExteriorPaths with type " << type << endl;
  if(type == "b"){
    Gp = &(*I).bGraph;
    Pp = &(*I).bExteriorPaths;
  }
  else if(type == "ad"){
    Gp = &(*I).adGraph;
    Pp = &(*I).adExteriorPaths;
  }
  else if(type == "all"){
    Gp = &(*I).allGraph;
    Pp = &(*I).allExteriorPaths;
  }
  
  ofstream gnupl((directory+"figures/" + fprefix+".exterioredges").c_str());
  if(type == "b" || type == "ad" || type == "all"){
    for(int i = 0; i < (*Pp).size();i++){
      for(int j = 0; j < (*Pp)[i].size()-1;j++){
	if(verbose){ 
	  cout << "Comp Num " << i << " Pos " << j  << " is Cell " << (*Gp)[(*Pp)[i][j]].CellNum << endl;
	  cout << (*Gp)[(*Pp)[i][j]].x << "\t" << (*Gp)[(*Pp)[i][j]].y << endl;
	  cout  << (*Gp)[(*Pp)[i][j+1]].x << "\t" << (*Gp)[(*Pp)[i][j+1]].y << endl;
	}
	cout.flush();
	gnupl << (*Gp)[(*Pp)[i][j]].x << "\t" << (*Gp)[(*Pp)[i][j]].y << endl;
	gnupl << (*Gp)[(*Pp)[i][j+1]].x << "\t" << (*Gp)[(*Pp)[i][j+1]].y << endl;
	gnupl << endl;
      }
    }
  }
  else{
    //printout b exterior paths
    if(verbose)cout << "printing b exterior paths...";
    Gp = &(*I).bGraph;
    Pp = &(*I).bExteriorPaths;
    for(int i = 0; i < (*Pp).size();i++){
      for(int j = 0; j < (*Pp)[i].size()-1;j++){
  	//if((*Pp)[i][j] <100000)
	gnupl << (*Gp)[(*Pp)[i][j]].x << "\t" << (*Gp)[(*Pp)[i][j]].y << endl;
  	// else{
  	//   int ippos = 0;
  	//   for(int ip = 0; ip < (*I).IPs.size();ip++){
  	//     if((*Pp)[i][j] == (*I).IPs[ip].CellNum){
  	//       ippos = ip;
  	//       break;
  	//     }
  	//   }
  	//   gnupl << (*I).IPs[ippos].x << "\t" << (*I).IPs[ippos].y << endl;
  	// }
  	// if((*Pp)[i][j+1] <100000)
	gnupl << (*Gp)[(*Pp)[i][j+1]].x << "\t" << (*Gp)[(*Pp)[i][j+1]].y << endl;
  	// else{
  	//   int ippos = 0;
  	//   for(int ip = 0; ip < (*I).IPs.size();ip++){
  	//     if((*Pp)[i][j+1] == (*I).IPs[ip].CellNum){
  	//       ippos = ip;
  	//       break;
  	//     }
  	//   }
  	//   gnupl << (*I).IPs[ippos].x << "\t" << (*I).IPs[ippos].y << endl;
  	// }
  	gnupl << endl;
      }
    }
  //   cout << "Finished\n" ;
     //printout ad exterior paths
     if(verbose)cout << "printing ad exterior paths...";
    Gp = &(*I).adGraph;
    Pp = &(*I).adExteriorPaths;
    for(int i = 0; i < (*Pp).size();i++){
      for(int j = 0; j < (*Pp)[i].size()-1;j++){
  // 	if((*Pp)[i][j] <100000)
	gnupl << (*Gp)[(*Pp)[i][j]].x << "\t" << (*Gp)[(*Pp)[i][j]].y << endl;
  // 	else{
  // 	  int ippos = 0;
  // 	  for(int ip = 0; ip < (*I).IPs.size();ip++){
  // 	    if((*Pp)[i][j] == (*I).IPs[ip].CellNum){
  // 	      ippos = ip;
  // 	      break;
  // 	    }
  // 	  }
  // 	  gnupl << (*I).IPs[ippos].x << "\t" << (*I).IPs[ippos].y << endl;
  // 	}
  // 	if((*Pp)[i][j+1] <100000)
	gnupl << (*Gp)[(*Pp)[i][j+1]].x << "\t" << (*Gp)[(*Pp)[i][j+1]].y << endl;
  // 	else{
  // 	  int ippos = 0;
  // 	  for(int ip = 0; ip < (*I).IPs.size();ip++){
  // 	    if((*Pp)[i][j+1] == (*I).IPs[ip].CellNum){
  // 	      ippos = ip;
  // 	      break;
  // 	    }
  // 	  }
  // 	  gnupl << (*I).IPs[ippos].x << "\t" << (*I).IPs[ippos].y << endl;
  // 	}
   	gnupl << endl;

      }
    }
  //   cout << "Finished\n";
  }
  gnupl.close();
  
  return;
}

//*********************************************************************************************//
//*********************************************************************************************//
//*********************************************************************************************//


void PrintoutGnuplot(Graph * Gp, string fprefix, string directory,string celltype){
  cout << "Check " << directory << fprefix << ".pl for a figure created." << endl;  
  
  PrintoutEdges(Gp,fprefix,directory);
  PrintoutVertices(Gp,fprefix,directory,celltype);  
  
  //create pl file
  ofstream pl((directory+"figures/"+fprefix+".pl").c_str());

  pl << "set term post eps enhanced color \"Arial\" 26 " << endl;
  pl << "set output '" <<  fprefix << ".eps'" << endl;
  
  pl << "set size 3,3" << endl;
  pl << "set multiplot" << endl;
  
  pl << "set style line 1 lt 1 pt 7 lc rgb \"blue\" lw 4 " << endl;
  pl << "set style line 2  lt 1 lc rgb \"red\" pt 5 " << endl;
  pl << "set style line 3  lt 1 lc rgb \"green\" pt 5 " << endl;
  pl << "set style line 4  lt 1 lc rgb \"brown\" pt 5 " << endl;
  pl << "set style line 5 lt 1 pt 7 lc rgb \"black\" lw 4 " << endl; 
  
  pl << "set origin 0,0" << endl;
  pl << "set size 3,3 " << endl;
  
  pl << "unset border" << endl;
  
  pl << "set title \"" <<  fprefix << "\"" << endl;
  pl << "plot '" << fprefix << ".edges' w l ls 5 notitle ,\\" << endl;
  pl << "    \"<awk '{if($4 == \\\"a\\\")print}' " <<  fprefix << ".abd.vertices\" w p ls 2 notitle,\\" << endl;
  pl << "    \"<awk '{if($4 == \\\"b\\\")print}' " <<  fprefix << ".abd.vertices\" w p ls 3 notitle,\\" << endl;
  pl << "    \"<awk '{if($4 == \\\"d\\\")print}' " <<  fprefix << ".abd.vertices\" w p ls 4 notitle,\\" << endl;
  pl << "    \"" <<  fprefix << ".abd.vertices\" u 1:2:($3) w labels font \"Arial,25\" offset 2  notitle" << endl;
  pl.close();
  return; 
}

//*********************************************************************************************//
//*********************************************************************************************//
//*********************************************************************************************//

void PrintoutGnuplotExteriorPaths(IsletGraph * I, string fprefix, string directory, string type){
  bool verbose=0;
  cout << "Check " << directory << fprefix << ".pl for a figure created." << endl; 

  //Printout edges and vertices
  if(verbose)cout << "Printing out edges... ";
  PrintoutEdges(&(*I).adGraph,fprefix,directory);
  bool append = 1;
  PrintoutEdges(&(*I).bGraph,fprefix,directory,append);
  if(verbose)cout << "Finished" << endl;
  if(verbose)cout << "Printing out vertices... ";
  PrintoutVertices(&(*I).allGraph,fprefix,directory,type);  
  if(verbose)cout << "Finished" << endl;
  if(verbose)cout << "Printing out exterior edges... ";
  PrintoutExteriorPaths(I, fprefix,directory);
  if(verbose)cout << "Finished" << endl;

  //create pl file
  cout << "Printing out pl file... ";
  ofstream pl((directory+"figures/"+fprefix+".pl").c_str());

  pl << "set term post eps enhanced color \"Arial\" 26 " << endl;
  pl << "set output '" <<  fprefix << ".eps'" << endl;
  
  pl << "set size 3,3" << endl;
  pl << "set multiplot" << endl;
  
  pl << "set style line 1 lt 1 pt 7 lc rgb \"blue\" lw 4 " << endl;
  pl << "set style line 2  lt 1 lc rgb \"red\" pt 5 " << endl;
  pl << "set style line 3  lt 1 lc rgb \"green\" pt 5 " << endl;
  pl << "set style line 4  lt 1 lc rgb \"brown\" pt 5 " << endl;
  pl << "set style line 5 lt 1 pt 7 lc rgb \"black\" lw 4 " << endl; 
  pl << "set style line 6 lt 1 pt 7 lc rgb \"black\" lw 8 " << endl;
  
  pl << "set origin 0,0" << endl;
  pl << "set size 3,3 " << endl;
  
  pl << "unset border" << endl;
  
  pl << "set title \"" <<  fprefix << "\"" << endl;
  pl << "plot '" << fprefix << ".edges' w l ls 5 notitle ,\\" << endl;
  pl << "    '" << fprefix << ".exterioredges' w l ls 6 notitle,\\" << endl;
  pl << "    \"<awk '{if($4 == \\\"a\\\")print}' " <<  fprefix << ".abd.vertices\" w p ls 2 notitle,\\" << endl;
  pl << "    \"<awk '{if($4 == \\\"b\\\")print}' " <<  fprefix << ".abd.vertices\" w p ls 3 notitle,\\" << endl;
  pl << "    \"<awk '{if($4 == \\\"d\\\")print}' " <<  fprefix << ".abd.vertices\" w p ls 4 notitle";//,\\" << endl;
  // pl << "    \"" <<  fprefix << ".abd.vertices\" u 1:2:($3) w labels font \"Arial,25\" offset 2  notitle" << endl;
  pl.close();
  if(verbose)cout << "Finished ";
  return; 
}

//*********************************************************************************************//
//*********************************************************************************************//
//*********************************************************************************************//

void PrintoutGnuplotIslet(IsletGraph* I, string fprefix, string directory,string option){
  cout << "Check " << directory << fprefix << ".pl for a figure created." << endl;  

  PrintoutVertices(&((*I).allGraph),fprefix,directory);
  PrintoutEdges(&((*I).adGraph),fprefix,directory);
  PrintoutEdges(&((*I).bGraph),fprefix,directory,1);
  
  //create pl file
  ofstream pl((directory+"figures/"+fprefix+".pl").c_str());

  pl << "set term post eps enhanced color \"Arial\" 26 " << endl;
  pl << "set output '" <<  fprefix << ".eps'" << endl;
  
  pl << "set size 3,3" << endl;
  pl << "set multiplot" << endl;
  
  pl << "set style line 1 lt 1 pt 7 lc rgb \"blue\" lw 4 " << endl;
  pl << "set style line 2  lt 1 lc rgb \"red\" pt 5 " << endl;
  pl << "set style line 3  lt 1 lc rgb \"green\" pt 5 " << endl;
  pl << "set style line 4  lt 1 lc rgb \"brown\" pt 5 " << endl;
  pl << "set style line 5 lt 1 pt 7 lc rgb \"black\" lw 4 " << endl; 
  
  pl << "set origin 0,0" << endl;
  pl << "set size 3,3 " << endl;
  
  pl << "unset border" << endl;
  
  pl << "set title \"" <<  fprefix << "\"" << endl;
  pl << "plot '" << fprefix << ".edges' w l ls 5 notitle ,\\" << endl;
  pl << "    \"<awk '{if($4 == \\\"a\\\")print}' " <<  fprefix << ".abd.vertices\" w p ls 2 notitle,\\" << endl;
  pl << "    \"<awk '{if($4 == \\\"b\\\")print}' " <<  fprefix << ".abd.vertices\" w p ls 3 notitle,\\" << endl;
  pl << "    \"<awk '{if($4 == \\\"d\\\")print}' " <<  fprefix << ".abd.vertices\" w p ls 4 notitle";
  if(option=="labels"){
    pl << ",\\" << endl;
    pl << "    \"" <<  fprefix << ".abd.vertices\" u 1:2:($3) w labels font \"Arial,25\" offset 2  notitle" << endl;
  }
  pl.close();
  return; 
}
//*********************************************************************************************//
//*********************************************************************************************//
//*********************************************************************************************//

void PrintoutGnuplotSetLoops(IsletGraph* I, int Gno, string fprefix, string directory, string type){
  bool verbose=0;
  cout << "Check " << directory << fprefix << ".pl for a figure created." << endl;  

  ////Printout edges and vertices
  //string sIsletNum = static_cast<ostringstream*>( &(ostringstream() << (*I).IsletNum) )->str();
  //if(verbose)cout << "Printing out edges... ";
  //PrintoutEdges(&(*I).adGraph,fprefix+".Islet"+sIsletNum,directory);
  //bool append = 1;
  //PrintoutEdges(&(*I).bGraph,fprefix+".Islet"+sIsletNum,directory,append);
  //if(verbose)cout << "Finished" << endl;
  //if(verbose)cout << "Printing out vertices... ";
  //PrintoutVertices(&(*I).allGraph,fprefix+".Islet"+sIsletNum,directory);  
  //if(verbose)cout << "Finished" << endl;
  //if(verbose)cout << "Printing out exterior edges... ";
  //PrintoutSetLoops(I, fprefix+".Islet"+sIsletNum,directory, type);
  //if(verbose)cout << "Finished" << endl;
  
  
  ////create pl file
  //int fignum = Gno/12, pos = Gno%12;
  //string sfignum = static_cast<ostringstream*>( &(ostringstream() << fignum) )->str();
  //
  //cout << "Printing out pl file... ";
  //ofstream pl((directory+"figures/"+fprefix+".fig" + sfignum +".pl").c_str(),ios::app);
  //
  //if(pos==0){
  //  //Printout top matter
  //  pl << "set term post eps enhanced color \"Arial\" 26" << endl;
  //  pl << "set output '" << fprefix << ".fig" << sfignum << ".eps'" << endl;
  //  pl << "" << endl;
  //  pl << "set size 3,4" << endl;
  //  pl << "set multiplot" << endl;
  //  pl << "" << endl;
  //  pl << "set style line 1 lt 1 pt 7 lc rgb \"blue\" lw 4 " << endl;
  //  pl << "set style line 2  lt 1 lc rgb \"red\" pt 7 ps 0.7" << endl;
  //  pl << "set style line 3  lt 1 lc rgb \"green\" pt 7 ps 0.7" << endl;
  //  //pl << "set style line 3  lt 1 lc rgb \"#9932CC\" pt 7 ps 0.7" << endl;//indigo
  //  pl << "set style line 4  lt 1 lc rgb \"brown\" pt 7 ps 0.7" << endl;
  //  pl << "set style line 5 lt 1 pt 7 lc rgb \"orange\" lw 4 " << endl;
  //  pl << "set style line 6 lt 1 pt 7 lc rgb \"black\" lw 2" << endl;
  //  pl << endl << endl;
  //}
  //
  //pl << "set origin " << pos%3 << "," << 3-pos/3 << endl;
  //pl << "set size 1,1 " << endl;
  //pl << "unset border" << endl;
  //pl << "set title \"Islet " << (*I).IsletNum <<": " << (*I).nalpha+(*I).nbeta+ (*I).ndelta << " cells\"" << endl;
  //
  //pl << "plot '" << fprefix <<".Islet" << sIsletNum <<  ".edges' w l ls 6 notitle ,\\" << endl;
  //pl << "    '" << fprefix <<".Islet" << sIsletNum <<  ".adsetLoops' w l ls 1 notitle,\\" << endl;
  //pl << "    '" << fprefix <<".Islet" << sIsletNum <<  ".betasetLoops' w l ls 5 notitle,\\" << endl;
  //pl << "    \"<awk '{if($4 == \\\"a\\\")print}' " <<  fprefix <<".Islet" << sIsletNum <<  ".abd.vertices\" w p ls 2 notitle,\\" << endl;
  //pl << "    \"<awk '{if($4 == \\\"b\\\")print}' " <<  fprefix <<".Islet" << sIsletNum <<  ".abd.vertices\" w p ls 3 notitle,\\" << endl;
  //pl << "    \"<awk '{if($4 == \\\"d\\\")print}' " <<  fprefix <<".Islet" << sIsletNum <<  ".abd.vertices\" w p ls 4 notitle,\\" << endl;
  //
  //pl << endl << endl;
  //pl.close();
  
  return;
}

//*********************************************************************************************//
//*********************************************************************************************//
//*********************************************************************************************//

void PrintoutVertices(Graph* Gp, string fprefix, string directory,string celltype,bool append){
   std::pair<vertex_iter, vertex_iter> vp;
   bool verbose = 0;
   ofstream gnupb;
   string gnupbfile = directory+"/figures/" + fprefix+"." + celltype+".vertices";
   if (verbose)cout << "Creating file " << gnupbfile << " to printout vertices" << endl;
   if(append)gnupb.open(gnupbfile.c_str(),ios::app);
   else gnupb.open(gnupbfile.c_str());

   //Printout points
   for(vp=vertices(*Gp); vp.first!=vp.second; vp.first++){
     if((*Gp)[*vp.first].CellNum < 99999){
       if((celltype == "abd") ||
	  (((*Gp)[*vp.first].type == "a") && (celltype == "aa" || celltype == "ab"|| celltype == "ad" || celltype == "adm" || celltype == "adb")) ||
	  (((*Gp)[*vp.first].type == "b") && (celltype == "ab" || celltype == "bb"|| celltype == "bd" || celltype == "adb")) ||
	  (((*Gp)[*vp.first].type == "d") && (celltype == "ad" || celltype == "bd"|| celltype == "dd" || celltype == "adm" || celltype == "adb"))){
	 gnupb << (*Gp)[*vp.first].x << "\t" << (*Gp)[*vp.first].y  << "\t"<< (*Gp)[*vp.first].IsletCellNum << "\t" << (*Gp)[*vp.first].type << "\n";
       }
     }
   }
   gnupb.close();
   return; 
}
//*********************************************************************************************//
//*********************************************************************************************//
//*********************************************************************************************//

void PrintoutGnuplotRandomlyChooseIsletsBySizeandFraction(vector<IsletGraph> *I,string fileprefix){
  int tmseed = 0;
  //tmseed = (int)(time0+simno);
  bool verbose = 0;
  srand(tmseed);
  
  //create pl file
  ofstream pl (("figures/"+fileprefix+".RandomlyChosenIslets.pl").c_str()); 
  if(verbose)cout << "Printing out pl file... ";cout.flush();
  
  //Printout top matter
  pl << "set term post eps enhanced color \"Arial\" 26" << endl;
  pl << "set output '" << fileprefix << ".RandomlyChosenIslets.eps'" << endl;
  pl << "" << endl;
  pl << "set size 4,5" << endl;
  pl << "set multiplot" << endl;
  pl << "" << endl;
  pl << "set style line 1 lt 1 pt 7 lc rgb \"blue\" lw 4 " << endl;
  pl << "set style line 2  lt 1 lc rgb \"red\" pt 7 ps 0.7" << endl;
  pl << "set style line 3  lt 1 lc \"green\" pt 7 ps 0.7" << endl;
  //pl << "set style line 3  lt 1 lc \"#9932CC\" pt 7 ps 0.5" << endl; //indigo
  pl << "set style line 4  lt 1 lc rgb \"brown\" pt 7 ps 0.7" << endl;
  pl << "set style line 5 lt 1 pt 7 lc rgb \"orange\" lw 4 " << endl;
  pl << "set style line 6 lt 1 pt 7 lc rgb \"black\" lw 1" << endl;
  pl << "set style line 7 lt 1 pt 7 lc rgb \"black\" lw 4" << endl;
  pl << endl << endl;
  
  vector<int> Islets;
  //Loop through Islet cell count
  for(int isletsize = 5;isletsize < 10;isletsize++){
    //Loop through Islet cell fraction
    for(double isletfrac = 0.2;isletfrac < 1.0;isletfrac+=0.2){
      for(int i = 0; i < (*I).size();i++){
	int num_verts = (*I)[i].nalpha+(*I)[i].nbeta+(*I)[i].ndelta;
	double cell_frac = (double)(*I)[i].nbeta/(double)num_verts;
	if((log2(num_verts) > isletsize && log2(num_verts) < isletsize+1) ||
	   (isletsize == 9 && log2(num_verts) > 9)||
	   (isletsize == 5 && log2(num_verts) >4 && log2(num_verts) < 6)){
	  if(cell_frac > isletfrac && cell_frac < isletfrac+0.2){
	    if(verbose)cout << " inside 2nd if" << endl; cout.flush();
	    Islets.push_back(i);
	  }
	}
      }
      if(verbose)cout << "finished loop" << endl; cout.flush();
      //Randomly choose an islet
      if(Islets.size() == 0)
	cout << "isletsize " << isletsize << " and isletfrac " << isletfrac << " has no entries" << endl; 
      else{
	int irand =Islets[rand()%Islets.size()];
	if(verbose)cout << "irand is " << irand << endl;cout.flush();
	int irandcount = (*I)[irand].nalpha+(*I)[irand].nbeta+ (*I)[irand].ndelta;
	if(verbose)cout << "irandcount is " << irandcount << endl;cout.flush();
	double irandfrac = (double)(*I)[irand].nbeta/(double)irandcount;
	if(verbose)cout << "irandfrac is " << irandfrac << endl;cout.flush();

	//string sIsletNum = static_cast<ostringstream*>( &(ostringstream() << (*I)[irand].IsletNum) )->str();
	//string directory = ".";
	//PrintoutVertices(&((*I)[irand].allGraph),fileprefix+".Islet"+sIsletNum,directory);
	//PrintoutEdges(&((*I)[irand].adGraph),fileprefix+".Islet"+sIsletNum,directory);
	//PrintoutEdges(&((*I)[irand].bGraph),fileprefix+".Islet"+sIsletNum,directory,1);
	//
	//pl << "set origin " << (int)(isletfrac/0.2)-1 << "," << isletsize-5 << endl;
	//pl << "set size 1,1 " << endl;
	//pl << "unset border" << endl;
	//pl << "unset xtics" << endl;
	//pl << "unset ytics" << endl;
	//pl << setprecision(3) << "set title \"Count = " <<irandcount  << ", Fraction = " << irandfrac << "\"" << endl;
	//
	//pl << "plot '" << fileprefix <<".Islet" << sIsletNum <<  ".edges' w l ls 6 notitle ,\\" << endl;
	//// pl << "    '" << fileprefix  << ".Islet" << sIsletNum<< ".exterioredges' w l ls 7 notitle,\\" << endl;
	//pl << "    \"<awk '{if($4 == \\\"a\\\")print}' " <<  fileprefix <<".Islet" << sIsletNum <<  ".abd.vertices\" w p ls 2 notitle,\\" << endl;
	//pl << "    \"<awk '{if($4 == \\\"b\\\")print}' " <<  fileprefix <<".Islet" << sIsletNum <<  ".abd.vertices\" w p ls 3 notitle,\\" << endl;
	//pl << "    \"<awk '{if($4 == \\\"d\\\")print}' " <<  fileprefix <<".Islet" << sIsletNum <<  ".abd.vertices\" w p ls 4 notitle" << endl;
	//pl << endl << endl;
      }
  //    //Clear out Islets vector
  //    int Isletvecsize = Islets.size();
  //    for(int j = 0;j < Isletvecsize;j++)Islets.pop_back();
    }
  }
  //
  //pl.close();
  Islets.clear(); vector<int> ().swap(Islets);
  
  return;
}

//*********************************************************************************************//
//*********************************************************************************************//
//*********************************************************************************************//

void PrintoutPartialLoopResults(IsletGraph * I, string fprefix, string directory){
  string resfile = directory + "figures/"+fprefix + ".PartialLoopsResults.txt";
  ofstream res;
  bool verbose=1;
  if(file_exists(resfile))res.open((resfile).c_str(),ios::app);
  else{
    res.open((resfile).c_str());
    res << "IsletNum\tIsletCellCount\tIsletbetaCellFraction\tComponentType\tComponentCellCount\tExteriorPathCount\tOrigSet?\tDoubleSet?\tMantlePercentage\tNumCellsinLoop\tLoop" << endl;
  }
  
  for(int c = 0; c < (*I).bComponents.size();c++){
    res <<  (*I).IsletNum << "\t" << (*I).nalpha+(*I).nbeta+(*I).ndelta << "\t" << (double)((*I).nbeta)/(double)((*I).nalpha+(*I).nbeta+(*I).ndelta) << "\tb\t" ;
    //add up non-IPs in component;
    int nComp = 0; 
    for(int cc = 0; cc < (*I).bComponents[c].size(); cc++)
      if((*I).bGraph[(*I).bComponents[c][cc]].IsletCellNum< 100000)nComp++;
    res << nComp << "\t";
    //add up non-IPs in exterior path
    int nExtPath = 0;
    for(int cc = 0; cc < (*I).bExteriorPaths[c].size();cc++)
      if((*I).bGraph[(*I).bExteriorPaths[c][cc]].IsletCellNum< 100000)nExtPath++;
    if((*I).bGraph[(*I).bExteriorPaths[c][0]].IsletCellNum < 100000)nExtPath--;
    res << nExtPath<< "\t" ;
    //if mantle is from an original set
    if((*I).bMantle_type[c]==1){
      res << "1\t0\t1\t" <<  "\t" << (*I).bMantle[c].size() << "\t" ;
      for(int bl = 0; bl < (*I).bMantle[c].size(); bl++)res <<(*I).allGraph[(*I).bMantle[c][bl]].IsletCellNum << " ";
      res << endl;
    }
    //if mantle is from a double set
    else if((*I).bMantle_type[c]==2){
      res << "0\t1\t1\t" << "\t" << (*I).bMantle[c].size() << "\t" ;
      for(int bl = 0; bl < (*I).bMantle[c].size(); bl++)res <<(*I).allGraph[(*I).bMantle[c][bl]].IsletCellNum << " ";
      res<< endl;
    }
    //if mantle is a Partial loop
    else if((*I).bMantle_type[c]==3){
      //add up percentages and average
      double PercSum = 0;
      if(nExtPath>0){
	for(int cc = 0; cc < (*I).bMantlePercent[c].size();cc++)
	  if((*I).bGraph[(*I).bExteriorPaths[c][cc]].IsletCellNum < 100000)
	    PercSum += (*I).bMantlePercent[c][cc];
	PercSum /= (double)nExtPath;
      }
      res << "0\t0\t" << PercSum <<"\t0" <<  endl;
    }
    if(verbose){
      //Printout results to screen
      cout << "BComponent " << c << " consists of " << endl;
      for(int i = 0; i < (*I).bComponents[c].size();i++)
	cout << (*I).bGraph[(*I).bComponents[c][i]].IsletCellNum << " ";
      cout << endl;
      cout << "BExteriorPath " << c << " consists of " << endl;
      for(int i = 0; i < (*I).bExteriorPaths[c].size();i++){
	if((*I).bExteriorPaths[c][i]<100000)cout << (*I).bGraph[(*I).bExteriorPaths[c][i]].IsletCellNum << " ";
	else cout << (*I).bExteriorPaths[c][i] << " ";
      }
      cout << endl;
      cout << "BMantle " << c << " consists of " << endl;
      for(int i = 0; i < (*I).bMantle[c].size();i++)
	cout << (*I).adGraph[(*I).bMantle[c][i]].IsletCellNum << " ";
      cout << endl;
      cout << "The bMantle Percentage is "  << endl;
      for(int i = 0; i < (*I).bMantlePercent[c].size();i++)
	cout << (*I).bMantlePercent[c][i] << " ";
      cout << endl;
      cout << "The beta mantle type is " << (*I).bMantle_type[c] << endl;
      
     
    }
  }
  
  for(int c = 0; c < (*I).adComponents.size();c++){
    res <<  (*I).IsletNum << "\t" << (*I).nalpha+(*I).nbeta+(*I).ndelta << "\t" << (double)((*I).nbeta)/(double)((*I).nalpha+(*I).nbeta+(*I).ndelta) << "\tad\t" ;
    //add up non-IPs in component;
    int nComp = 0; 
    for(int cc = 0; cc < (*I).adComponents[c].size(); cc++)
      if((*I).adGraph[(*I).adComponents[c][cc]].IsletCellNum< 100000)nComp++;
    res << nComp << "\t";
    //add up non-IPs in exterior path
    int nExtPath = 0;
    for(int cc = 0; cc < (*I).adExteriorPaths[c].size();cc++)
      if((*I).adGraph[(*I).adExteriorPaths[c][cc]].IsletCellNum< 100000)nExtPath++;
    if((*I).adGraph[(*I).adExteriorPaths[c][0]].IsletCellNum < 100000)nExtPath--;
    res << nExtPath<< "\t" ;
    //if mantle is from an original set
    if((*I).adMantle_type[c]==1){
      res << "1\t0\t1\t" <<  "\t" << (*I).adMantle[c].size() << "\t" ;
      for(int adl = 0; adl < (*I).adMantle[c].size(); adl++)res <<(*I).allGraph[(*I).adMantle[c][adl]].IsletCellNum << " ";
      res << endl;
    }
    //if mantle is from a double set
    else if((*I).adMantle_type[c]==2){
      res << "0\t1\t1\t" << "\t" << (*I).adMantle[c].size() << "\t" ;
      for(int adl = 0; adl < (*I).adMantle[c].size(); adl++)res <<(*I).allGraph[(*I).adMantle[c][adl]].IsletCellNum << " ";
      res<< endl;
    }
    //if mantle is a Partial loop
    else if((*I).adMantle_type[c]==3){
      //add up percentages and average
      double PercSum = 0;
      if(nExtPath>0){
	for(int cc = 0; cc < (*I).adMantlePercent[c].size();cc++)
	  if((*I).adGraph[(*I).adExteriorPaths[c][cc]].IsletCellNum < 100000)
	    PercSum += (*I).adMantlePercent[c][cc];
	PercSum /= (double)nExtPath;
      }
      res << "0\t0\t" << PercSum <<"\t0" <<  endl;
    }
    if(verbose){
      //Printout results to screen
      cout << "ADComponent " << c << " consists of " << endl;
      for(int i = 0; i < (*I).adComponents[c].size();i++)
	cout << (*I).adGraph[(*I).adComponents[c][i]].IsletCellNum << " ";
      cout << endl;
      cout << "ADExteriorPath " << c << " consists of " << endl;
      for(int i = 0; i < (*I).adExteriorPaths[c].size();i++){
	if((*I).adExteriorPaths[c][i]<100000)cout << (*I).adGraph[(*I).adExteriorPaths[c][i]].IsletCellNum << " ";
	else cout << (*I).adExteriorPaths[c][i] << " ";
      }
      cout << endl;
      cout << "ADMantle " << c << " consists of " << endl;
      for(int i = 0; i < (*I).adMantle[c].size();i++)
	cout << (*I).adGraph[(*I).adMantle[c][i]].IsletCellNum << " ";
      cout << endl;
      cout << "The adMantle Percentage is "  << endl;
      for(int i = 0; i < (*I).adMantlePercent[c].size();i++)
	cout << (*I).adMantlePercent[c][i] << " ";
      cout << endl;
      cout << "The adeta mantle type is " << (*I).adMantle_type[c] << endl;
      
      
    }
  }
  res.close();
  return;
}

//*********************************************************************************************//
//*********************************************************************************************//
//*********************************************************************************************//

void PrintoutSetLoops(IsletGraph *I, string fprefix,string directory, string type ){
 bool verbose=0;
  Graph *Gp;
  vector<vector<vertex_t> > * Pp;
  vector<vertex_t> *Hp;
  if(verbose)cout << "Inside PrintoutExteriorPaths with type " << type << endl;
  if(type == "beta"){
    Gp = &(*I).adGraph;
    Pp = &(*I).betasetLoops;
  }
  else if(type == "ad"){
    Gp = &(*I).bGraph;
    Pp = &(*I).adsetLoops;
  }
   
  if(type == "beta" || type == "ad" ){
    ofstream gnupl((directory+"figures/" + fprefix+"." + type+"setLoops").c_str());
    for(int i = 0; i < (*Pp).size();i++){
      for(int j = 0; j < (*Pp)[i].size()-1;j++){
	      gnupl << (*Gp)[(*Pp)[i][j]].x << "\t" << (*Gp)[(*Pp)[i][j]].y << endl;
	      gnupl << (*Gp)[(*Pp)[i][j+1]].x << "\t" << (*Gp)[(*Pp)[i][j+1]].y << endl;
	      gnupl << endl;
      }
    }
    gnupl.close();
  }
  
  else{
    //printout betaset loops
    ofstream gnupl((directory+"figures/" + fprefix+".betasetLoops").c_str());
    if(verbose)cout << "printing beta set loops...";
    Gp = &(*I).adGraph;
    Pp = &(*I).betasetLoops;
    for(int i = 0; i < (*Pp).size();i++){
      for(int j = 0; j < (*Pp)[i].size()-1;j++){
	      gnupl << (*Gp)[(*Pp)[i][j]].x << "\t" << (*Gp)[(*Pp)[i][j]].y << endl;
	      gnupl << (*Gp)[(*Pp)[i][j+1]].x << "\t" << (*Gp)[(*Pp)[i][j+1]].y << endl;
	      gnupl << endl;
      }
    }
    gnupl.close();
    cout << "Finished\n" ;
    //printout ad set loops
    if(verbose)cout << "printing ad set loops...";
    gnupl.open((directory+"figures/" + fprefix+".adsetLoops").c_str());
    Gp = &(*I).bGraph;
    Pp = &(*I).adsetLoops;
    for(int i = 0; i < (*Pp).size();i++){
      for(int j = 0; j < (*Pp)[i].size()-1;j++){
	      gnupl << (*Gp)[(*Pp)[i][j]].x << "\t" << (*Gp)[(*Pp)[i][j]].y << endl;
	      gnupl << (*Gp)[(*Pp)[i][j+1]].x << "\t" << (*Gp)[(*Pp)[i][j+1]].y << endl;
	      gnupl << endl;
      }
    }
    cout << "Finished\n";
    gnupl.close();
  }
  return;
}




//*********************************************************************************************//
//*********************************************************************************************//
//*********************************************************************************************//



 void PrintoutGnuplotGofr4x5(vector<vector<double> > *Gofr, IsletGraph* I, string fig_directory, string fileprefix, string type, double min){

   bool verbose=0;

   //string ssubjno = static_cast<ostringstream*>( &(ostringstream() << subjno) )->str();

   string this_file = fileprefix+"."+type+".gr";

   string tfile = fig_directory+"/"+this_file + ".pl"; 

   cout << "Check " << tfile << " for a figure created." << endl;  

   string tfile_gr = fig_directory+"/"+this_file; 

   //string group = getGroup(fprefix);

   //Printout Gr file
   //string subjprefix = fprefix+".subj"+ssubjno;
   
   PrintoutSubjectGofr(Gofr, tfile_gr);

   ////Set xend to be the last value above one
   //double xendf =(*Gofr).size()-1;
   //while ((*Gofr)[xendf][1] < 1) xendf--;
   //double xend = (*Gofr)[xendf][0]+4;


   //create pl file
   //int fignum = inum/20, pos = inum%20;
   
   //string sfignum = static_cast<ostringstream*>( &(ostringstream() << fignum) )->str();

   cout << "Printing out pl file " << tfile << "... ";

   ofstream pl(tfile.c_str(), ios::app);

   //if(pos==0){
   
   //Printout top matter
   pl << "set term post eps enhanced color \"Arial\" 26" << endl;
   pl << "set output '" << this_file << ".eps'" << endl;
   //pl << "" << endl;
   //pl << "set size 4,5" << endl;
   //pl << "set multiplot" << endl;
   pl << "" << endl;
   pl << "set style line 1 lt 1 pt 7 lc rgb \"black\" lw 8" << endl;
   pl << "set style line 2 lt 1 pt 7 lc rgb \"red\" lw 8" << endl;
   pl << endl << endl;

   //}

   //pl << "set origin " << pos%4 << "," << 4-pos/4 << endl;
   pl << "set origin " << 0 << "," << 0 << endl;

   pl << "set size 1,1 " << endl;
   pl << "unset arrow" << endl;
   //pl << "set arrow from 0,1 to " << xend << ",1 nohead lt -1 lw 2 #creates line at y=1" << endl;
   pl << "set arrow from 0,1 to 100,1 nohead lt -1 lw 2 #creates line at y=1" << endl;
   pl << "set arrow from " << min << ",0.5 to " << min << ", 1.5 nohead ls 2  " << endl;
   pl << "set xtics font \"Arial-Bold,30\"" << endl;
   pl << "set ytics font \"Arial-Bold,30\"" << endl;
   pl << "unset border" << endl;
   //pl << "set xrange [0:" << xend <<"]" << endl;
   pl << "set xrange [0:100]" << endl;
   pl << "set title \"Islet # " << (*I).IsletNum<< "(" ;
   if(type == "bb") pl << (*I).nbeta << " b cells)\"" <<  endl;
   if(type == "adb") pl << (*I).nbeta+(*I).nalpha+(*I).ndelta << " total cells)\"" <<  endl;
   if(type == "adm") pl << (*I).nalpha+(*I).ndelta << " ad cells)\"" <<  endl;

   //pl << "plot '" << subjprefix <<".subj" << (*I).IsletNum << "." << type << ".gr' u 1:2  w l ls 1 notitle" << endl; 
   pl << "plot '" << this_file << "' u 1:2  w l ls 1 notitle" << endl; 
   pl << endl << endl;
   pl.close();

   return;

 }

//*********************************************************************************************//
//*********************************************************************************************//
//*********************************************************************************************//





// // void PrintoutTsvfile(vector<IsletGraph>* I, string fileprefix){
// //   ofstream fout((fileprefix+".tsv").c_str());
// //   std::pair<vertex_iter, vertex_iter> vp;
// //   for(int i = 0; i < (*I).size();i++){
// //     for(vp = vertices((*I)[i].allGraph); vp.first!=vp.second;vp.first++)
// //       fout << (*I)[i].IsletNum << "\t" << (*I)[i].allGraph[*vp.first].IsletCellNum << "\t" << (*I)[i].allGraph[*vp.first].x << "\t" << (*I)[i].allGraph[*vp.first].y << "\t" << (*I)[i].allGraph[*vp.first].type<< "\t" << (*I)[i].IsletArea << endl;
// //   }
// //   fout.close();
// //   return;
// // }

// // void PrintoutSubjectSpecificInfo(vector<IsletGraph>* I,string fileprefix,string subjprefix){
// //   int vectorcount = 0;
// //   //calculate beta cell-fraction and size of each islet
// //   vector<double> bfrac((*I).size(),0.0);vectorcount++;
// //   vector<double> log2size((*I).size(),0.0);vectorcount++;
// //   int nbetaIslets = 0, nadmIslets = 0,nadbIslets=0;
// //   for(int i = 0; i < (*I).size();i++){
// //     double numcells = (double)(num_vertices((*I)[i].allGraph));
// //     bfrac[i] = (double)(*I)[i].nbeta/numcells;
// //     log2size[i] = log2(numcells);
// //     if((*I)[i].nbeta > 5)nbetaIslets++;
// //     if((*I)[i].nalpha + (*I)[i].ndelta>5)nadmIslets++;
// //     if((*I)[i].nbeta > 5 && (*I)[i].nalpha + (*I)[i].ndelta>5)nadbIslets++;
// //   }
// //   double bfrac_mean = 0, bfrac_stdev = 0, log2size_mean = 0, log2size_stdev = 0;
// //   CalculateMeanAndStdev(&bfrac,&bfrac_mean,&bfrac_stdev);
// //   CalculateMeanAndStdev(&log2size,&log2size_mean,&log2size_stdev);
// //   RemoveVector(&bfrac);vectorcount--;
// //   RemoveVector(&log2size);vectorcount--;

// //   //Printout results
// //   string foutname = "Dev.Subject_bcellfraction_log2size";
// //   ofstream fout;
// //   if(file_exists(foutname))
// //     fout.open(foutname.c_str(),ios::app);
// //   else{
// //     fout.open(foutname.c_str());
// //     fout << "Age\tSubject\t#Islets\t#Islets(b>5)\t#Islets(a+d>5)\t#Islets(both)\tbcell fraction mean\tbcell fractions stdev\tlog2size mean\tlog2size stdev" << endl;
// //   }
// //   fout << fileprefix << "\t" << subjprefix << "\t" << (*I).size() << "\t" <<nbetaIslets << "\t" << nadmIslets << "\t" << nadbIslets << "\t" <<  bfrac_mean << "\t" << bfrac_stdev << "\t" << log2size_mean << "\t" << log2size_stdev << endl;  
// //   fout.close();
  
// //   return;
// // }


// // void PrintoutIPVertices(Graph* g, string filename){
// //   ofstream fout (filename.c_str());
// //   //cout << "Printing out " << filename << endl;
// //   for(int i = 0; i < num_vertices(*g);i++){
// //     if((*g)[i].CellNum > 99999){
// //       fout << (*g)[i].CompNum << "\t" <<(*g)[i].CellNum <<  "\t" << (*g)[i].x << "\t" << (*g)[i].y <<"\t" << (*g)[i].IPparents[0] << "\t" << (*g)[i].IPparents[1]<< "\t" << (*g)[i].IPparents[2]<< "\t" << (*g)[i].IPparents[3]<< endl;
// //     }
// //   }
// //   fout.close();
// //   return;
// // }

// // void PrintoutIPEdges(Graph* g, string filename){
// //   ofstream fout (filename.c_str());
// //   std::pair<edge_iter, edge_iter> ep; 
// //   //cout << "Printing out " << filename << endl;
// //   for(ep = edges(*g);ep.first!=ep.second;ep.first++)
// //     fout << (*g)[source(*ep.first,*g)].CellNum << "\t" << (*g)[target(*ep.first,*g)].CellNum << "\t" <<(*g)[source(*ep.first,*g)].IsletCellNum << "\t" << (*g)[target(*ep.first,*g)].IsletCellNum <<"\t" << (*g)[*ep.first].distance << endl;
// //   fout.close();
// //   return;

// // }





// // void PrintoutLoopResults(vector<IsletGraph> * Islets,string fprefix, string directory ){
// //   int vectorcount = 0;
// //   bool verbose=0;
// //   if(verbose)cout << "inside Printout Label results..." << endl;
// //   ofstream ressize((directory + "figures/"+fprefix + ".LogLoopResults.txt").c_str());
// //   ofstream resfrac((directory + "figures/"+fprefix + ".CellFractionLoopResults.txt").c_str());
  

// //   ressize << "numcells(2^x)\tnumIslets\t#ad cells in sets\t#ad Loops (b cells) \t#bb cells in sets \t#b Loops (ad cells)" << endl;
// //   resfrac << "Cell Fraction(b/total)\tnumIslets\t#ad cells in sets\t#ad Loops (b cells) \t#bb cells in sets \t#b Loops (ad cells)" << endl;
// //   vector<int> nIsletssize(11,0),nIsletsfrac(11,0);vectorcount+=2;
// //   vector<vector<int> > nadsetloopssize, nbsetloopssize,ncellsinadsetssize,ncellsinbsetssize;vectorcount+=4;
// //   vector<vector<int> > nadsetloopsfrac, nbsetloopsfrac,ncellsinadsetsfrac,ncellsinbsetsfrac;vectorcount+=4;
// //   vector<vector<int> > nadsetloopssizefrac(10,vector<int>(11,0)),nbsetloopssizefrac(10,vector<int>(11,0)),nIsletssizefrac(10,vector<int>(11,0)),ncellsinadsetssizefrac(10,vector<int>(11,0)),ncellsinbsetssizefrac(10,vector<int>(11,0));vectorcount+=5;

// //   //Setup vectors
// //   for(int i = 0; i < 11;i++){
// //     vector<int>temp;vectorcount++;
// //     nadsetloopssize.push_back(temp); nbsetloopssize.push_back(temp);ncellsinadsetssize.push_back(temp);ncellsinbsetssize.push_back(temp);
// //     nadsetloopsfrac.push_back(temp); nbsetloopsfrac.push_back(temp);ncellsinadsetsfrac.push_back(temp);ncellsinbsetsfrac.push_back(temp);
// //     RemoveVector(&temp);vectorcount--;
// //   }
// //   for(int g = 0; g < (*Islets).size();g++){
// //     if((*Islets)[g].nbeta > 5  && (*Islets)[g].nalpha + (*Islets)[g].ndelta > 5){
// //       int lognumcells = (int)log2(num_vertices((*Islets)[g].allGraph));
// //       int CellFraction = (int)((double)(*Islets)[g].nbeta/(double)num_vertices((*Islets)[g].allGraph)*10.0);
// //       if(lognumcells > 9) lognumcells = 9;
// //       nIsletssize[lognumcells]++;
// //       nIsletsfrac[CellFraction]++;
// //       nIsletssizefrac[lognumcells][CellFraction]++;
      
// //       nadsetloopssizefrac[lognumcells][CellFraction]+=(*Islets)[g].adsetLoops.size();
// //       nadsetloopssize[lognumcells].push_back((*Islets)[g].adsetLoops.size());
// //       nadsetloopsfrac[CellFraction].push_back((*Islets)[g].adsetLoops.size());
// //       nbsetloopssizefrac[lognumcells][CellFraction]+=(*Islets)[g].betasetLoops.size();
// //       nbsetloopssize[lognumcells].push_back((*Islets)[g].betasetLoops.size());
// //       nbsetloopsfrac[CellFraction].push_back((*Islets)[g].betasetLoops.size());
// //       for(int i = 0; i < (*Islets)[g].adsets.size();i++){
// // 	ncellsinadsetssize[lognumcells].push_back((*Islets)[g].adsets[i].size());
// // 	ncellsinadsetsfrac[CellFraction].push_back((*Islets)[g].adsets[i].size());
// // 	ncellsinadsetssizefrac[lognumcells][CellFraction]+=(*Islets)[g].adsets[i].size();
// //       }
// //       for(int i = 0; i < (*Islets)[g].betasets.size();i++){
// // 	ncellsinbsetssize[lognumcells].push_back((*Islets)[g].betasets[i].size());
// // 	ncellsinbsetsfrac[CellFraction].push_back((*Islets)[g].betasets[i].size());
// // 	ncellsinbsetssizefrac[lognumcells][CellFraction]+=(*Islets)[g].betasets[i].size();
// //       }
// //     }
// //   }
  
// //   for(int i = 3; i < 10; i++){
// //     //Printout sizeresults
// //     if(i == 9) ressize << "9+\t" << nIsletssize[i];
// //     else ressize << i << "\t" << nIsletssize[i];
// //     int sum = 0;
// //     for(int j = 0; j < ncellsinadsetssize[i].size();j++)sum+=ncellsinadsetssize[i][j];
// //     ressize  << "\t" << sum;
// //     sum = 0;
// //     for(int j = 0; j < nadsetloopssize[i].size();j++)sum+=nadsetloopssize[i][j];
// //     ressize  << "\t" << sum;
// //     sum = 0;
// //     for(int j = 0; j < ncellsinbsetssize[i].size();j++)sum+=ncellsinbsetssize[i][j];
// //     ressize  << "\t" << sum;
// //     sum = 0;
// //     for(int j = 0; j < nbsetloopssize[i].size();j++)sum+=nbsetloopssize[i][j];
// //     ressize  << "\t" << sum << endl;
    
// //     //Printout files
// //     string si = static_cast<ostringstream*>( &(ostringstream() << i) )->str();
// //     string file = directory + "figures/"+fprefix + ".adsets.size"+si+".LogLoopResults.txt";
// //     ofstream fout(file.c_str());
// //     for(int j = 0; j < ncellsinadsetssize[i].size(); j++)fout << ncellsinadsetssize[i][j] << endl;
// //     fout.close();
// //     file = directory + "figures/"+fprefix + ".adsetloops.size"+si+".LogLoopResults.txt";
// //     fout.open(file.c_str());
// //     for(int j = 0; j < nadsetloopssize[i].size(); j++)fout << nadsetloopssize[i][j] << endl;
// //     fout.close();
// //     file = directory + "figures/"+fprefix + ".bsets.size"+si+".LogLoopResults.txt";
// //     fout.open(file.c_str());
// //     for(int j = 0; j < ncellsinbsetssize[i].size(); j++)fout << ncellsinbsetssize[i][j] << endl;
// //     fout.close();
// //     file = directory + "figures/"+fprefix + ".bsetloops.size"+si+".LogLoopResults.txt";
// //     fout.open(file.c_str());
// //     for(int j = 0; j < nbsetloopssize[i].size(); j++)fout << nbsetloopssize[i][j] << endl;
// //     fout.close();
// //   }
// //   ressize.close();
// //   for(int i = 0; i < 11; i++){
// //     //Printout frac results
// //     resfrac << i << "\t" << nIsletsfrac[i];
// //     int sum = 0;
// //     for(int j = 0; j < ncellsinadsetsfrac[i].size();j++)sum+=ncellsinadsetsfrac[i][j];
// //     resfrac  << "\t" << sum;
// //     sum = 0;
// //     for(int j = 0; j < nadsetloopsfrac[i].size();j++)sum+=nadsetloopsfrac[i][j];
// //     resfrac  << "\t" << sum;
// //     sum = 0;
// //     for(int j = 0; j < ncellsinbsetsfrac[i].size();j++)sum+=ncellsinbsetsfrac[i][j];
// //     resfrac  << "\t" << sum;
// //     sum = 0;
// //     for(int j = 0; j < nbsetloopsfrac[i].size();j++)sum+=nbsetloopsfrac[i][j];
// //     resfrac  << "\t" << sum << endl;
    
// //     //Printout files
// //     string si = static_cast<ostringstream*>( &(ostringstream() << i) )->str();
// //     string file = directory + "figures/"+fprefix + ".adsets.frac"+si+".LogLoopResults.txt";
// //     ofstream fout(file.c_str());
// //     for(int j = 0; j < ncellsinadsetsfrac[i].size(); j++)fout << ncellsinadsetsfrac[i][j] << endl;
// //     fout.close();
// //     file = directory + "figures/"+fprefix + ".adsetloops.frac"+si+".LogLoopResults.txt";
// //     fout.open(file.c_str());
// //     for(int j = 0; j < nadsetloopsfrac[i].size(); j++)fout << nadsetloopsfrac[i][j] << endl;
// //     fout.close();
// //     file = directory + "figures/"+fprefix + ".bsets.frac"+si+".LogLoopResults.txt";
// //     fout.open(file.c_str());
// //     for(int j = 0; j < ncellsinbsetsfrac[i].size(); j++)fout << ncellsinbsetsfrac[i][j] << endl;
// //     fout.close();
// //     file = directory + "figures/"+fprefix + ".bsetloops.frac"+si+".LogLoopResults.txt";
// //     fout.open(file.c_str());
// //     for(int j = 0; j < nbsetloopsfrac[i].size(); j++)fout << nbsetloopsfrac[i][j] << endl;
// //     fout.close();
// //   }
// //   resfrac.close();
 

// //   ofstream ressizefracnislets((directory + "figures/"+fprefix + ".IsletCount.CellFractionvsSize.txt").c_str());
// //    for(int i = 0; i < nIsletssizefrac.size();i++){
// //      for(int j = 0; j <nIsletssizefrac[i].size();j++){
       
// // 	ressizefracnislets  << (double)nIsletssizefrac[i][j] << "\t";
       
// //      }
// //     ressizefracnislets  << endl;
// //    }
// //    ressizefracnislets.close();

// //    ofstream ressizefracadsets((directory + "figures/"+fprefix + ".adsets.CellFractionvsSize.txt").c_str());
// //    for(int i = 0; i < ncellsinadsetssizefrac.size();i++){
// //      for(int j = 0; j < ncellsinadsetssizefrac[i].size();j++){
// //        if(nIsletssizefrac[i][j] > 0)
// // 	 ressizefracadsets << (double)ncellsinadsetssizefrac[i][j]/(double)nIsletssizefrac[i][j] << "\t";
// //        else ressizefracadsets << "0\t";
// //      }
// //      ressizefracadsets << endl;
// //    }
// //    ressizefracadsets.close();
// //    ofstream ressizefracadloops((directory + "figures/"+fprefix + ".adloops.CellFractionvsSize.txt").c_str());
// //    for(int i = 0; i < nadsetloopssizefrac.size();i++){
// //      for(int j = 0; j < nadsetloopssizefrac[i].size();j++){
// //        if(nIsletssizefrac[i][j] > 0)
// // 	 ressizefracadloops << (double)nadsetloopssizefrac[i][j]/(double)nIsletssizefrac[i][j] << "\t";
// //        else ressizefracadloops << "0\t";
// //      }
// //      ressizefracadloops << endl;
// //    }
// //    ressizefracadloops.close();
// //    ofstream ressizefracbsets((directory + "figures/"+fprefix + ".bsets.CellFractionvsSize.txt").c_str());
// //    for(int i = 0; i < ncellsinbsetssizefrac.size();i++){
// //      for(int j = 0; j < ncellsinbsetssizefrac[i].size();j++){
// //        if(nIsletssizefrac[i][j] > 0)
// // 	 ressizefracbsets << (double)ncellsinbsetssizefrac[i][j]/(double)nIsletssizefrac[i][j] << "\t";
// //        else ressizefracbsets << "0\t";
// //      }
// //      ressizefracbsets << endl;
// //    }
// //    ressizefracbsets.close();
// //    ofstream ressizefracbloops((directory + "figures/"+fprefix + ".bloops.CellFractionvsSize.txt").c_str());
// //    for(int i = 0; i < nbsetloopssizefrac.size();i++){
// //      for(int j = 0; j < nbsetloopssizefrac[i].size();j++){
// //        if(nIsletssizefrac[i][j] > 0)
// // 	 ressizefracbloops << (double)nbsetloopssizefrac[i][j]/(double)nIsletssizefrac[i][j] << "\t";
// //        else ressizefracbloops << "0\t";
// //      }
// //      ressizefracbloops << endl;
// //    }
// //    ressizefracbloops.close();
   
  

 
// //    RemoveVector2d(&nadsetloopssize); RemoveVector2d(&nbsetloopssize); RemoveVector(&nIsletssize); RemoveVector2d(&ncellsinadsetssize); RemoveVector2d(&ncellsinbsetssize);vectorcount-=5;
// //    RemoveVector2d(&nadsetloopsfrac); RemoveVector2d(&nbsetloopsfrac); RemoveVector(&nIsletsfrac); RemoveVector2d(&ncellsinadsetsfrac); RemoveVector2d(&ncellsinbsetsfrac);vectorcount-=5;
// //    RemoveVector2d(&nadsetloopssizefrac);RemoveVector2d(&nbsetloopssizefrac),RemoveVector2d(&nIsletssizefrac),RemoveVector2d(&ncellsinadsetssizefrac),RemoveVector2d(&ncellsinbsetssizefrac);vectorcount-=5;
// //    if(vectorcount!=0)cout << "Problem with PrintoutLoopResults: vectorcount = " << vectorcount << endl;cout.flush();
// //    return;
// //  }

// //  //*********************************************************************************************//
// // //*********************************************************************************************//
// // //*********************************************************************************************//

// // void PrintoutGnuplotIslet(vector<IsletGraph> * I, int ino, string fprefix, string directory,string option){
// //   cout << "Check " << directory << fprefix << ".pl for a figure created." << endl;  

// //   PrintoutVertices(&((*I)[ino].allGraph),fprefix,directory);
// //   PrintoutEdges(&((*I)[ino].adGraph),fprefix,directory);
// //   PrintoutEdges(&((*I)[ino].bGraph),fprefix,directory,1);
  
// //   //create pl file
// //   ofstream pl((directory+"figures/"+fprefix+".pl").c_str());

// //   pl << "set term post eps enhanced color \"Arial\" 26 " << endl;
// //   pl << "set output '" <<  fprefix << ".eps'" << endl;
  
// //   pl << "set size 3,3" << endl;
// //   pl << "set multiplot" << endl;
  
// //   pl << "set style line 1 lt 1 pt 7 lc rgb \"blue\" lw 4 " << endl;
// //   pl << "set style line 2  lt 1 lc rgb \"red\" pt 5 " << endl;
// //   pl << "set style line 3  lt 1 lc rgb \"green\" pt 5 " << endl;
// //   pl << "set style line 4  lt 1 lc rgb \"brown\" pt 5 " << endl;
// //   pl << "set style line 5 lt 1 pt 7 lc rgb \"black\" lw 4 " << endl; 
  
// //   pl << "set origin 0,0" << endl;
// //   pl << "set size 3,3 " << endl;
  
// //   pl << "unset border" << endl;
  
// //   pl << "set title \"" <<  fprefix << "\"" << endl;
// //   pl << "plot '" << fprefix << ".edges' w l ls 5 notitle ,\\" << endl;
// //   pl << "    \"<awk '{if($4 == \\\"a\\\")print}' " <<  fprefix << ".abd.vertices\" w p ls 2 notitle,\\" << endl;
// //   pl << "    \"<awk '{if($4 == \\\"b\\\")print}' " <<  fprefix << ".abd.vertices\" w p ls 3 notitle,\\" << endl;
// //   pl << "    \"<awk '{if($4 == \\\"d\\\")print}' " <<  fprefix << ".abd.vertices\" w p ls 4 notitle";
// //   if(option=="labels"){
// //     pl << ",\\" << endl;
// //     pl << "    \"" <<  fprefix << ".abd.vertices\" u 1:2:($3) w labels font \"Arial,25\" offset 2  notitle" << endl;
// //   }
// //   pl.close();
// //   return; 
// // }
// //  //*********************************************************************************************//
// // //*********************************************************************************************//
// // //*********************************************************************************************//



// // void PrintoutGnuplot(Graph * Gp, string fprefix, string directory,string celltype){
// //   cout << "Check " << directory << fprefix << ".pl for a figure created." << endl;  

// //   PrintoutEdges(Gp,fprefix,directory);
// //   PrintoutVertices(Gp,fprefix,directory,celltype);  
  
// //   //create pl file
// //   ofstream pl((directory+"figures/"+fprefix+".pl").c_str());

// //   pl << "set term post eps enhanced color \"Arial\" 26 " << endl;
// //   pl << "set output '" <<  fprefix << ".eps'" << endl;
  
// //   pl << "set size 3,3" << endl;
// //   pl << "set multiplot" << endl;
  
// //   pl << "set style line 1 lt 1 pt 7 lc rgb \"blue\" lw 4 " << endl;
// //   pl << "set style line 2  lt 1 lc rgb \"red\" pt 5 " << endl;
// //   pl << "set style line 3  lt 1 lc rgb \"green\" pt 5 " << endl;
// //   pl << "set style line 4  lt 1 lc rgb \"brown\" pt 5 " << endl;
// //   pl << "set style line 5 lt 1 pt 7 lc rgb \"black\" lw 4 " << endl; 
  
// //   pl << "set origin 0,0" << endl;
// //   pl << "set size 3,3 " << endl;
  
// //   pl << "unset border" << endl;
  
// //   pl << "set title \"" <<  fprefix << "\"" << endl;
// //   pl << "plot '" << fprefix << ".edges' w l ls 5 notitle ,\\" << endl;
// //   pl << "    \"<awk '{if($4 == \\\"a\\\")print}' " <<  fprefix << ".abd.vertices\" w p ls 2 notitle,\\" << endl;
// //   pl << "    \"<awk '{if($4 == \\\"b\\\")print}' " <<  fprefix << ".abd.vertices\" w p ls 3 notitle,\\" << endl;
// //   pl << "    \"<awk '{if($4 == \\\"d\\\")print}' " <<  fprefix << ".abd.vertices\" w p ls 4 notitle,\\" << endl;
// //   pl << "    \"" <<  fprefix << ".abd.vertices\" u 1:2:($3) w labels font \"Arial,25\" offset 2  notitle" << endl;
// //   pl.close();
// //   return; 
// // }

// // //*********************************************************************************************//
// // //*********************************************************************************************//
// // //*********************************************************************************************//

// // void PrintoutGnuplotWithIPs(Graph * Gp, string fprefix, string directory,string celltype){
// //   cout << "Check " << directory << fprefix << ".pl for a figure created." << endl;  

// //   PrintoutEdges(Gp,fprefix,directory);
// //   PrintoutVerticesWithIPs(Gp,fprefix,directory,celltype);  
  
// //   //create pl file
// //   ofstream pl((directory+"figures/"+fprefix+".pl").c_str());

// //   pl << "set term post eps enhanced color \"Arial\" 26 " << endl;
// //   pl << "set output '" <<  fprefix << ".eps'" << endl;
  
// //   pl << "set size 3,3" << endl;
// //   pl << "set multiplot" << endl;
  
// //   pl << "set style line 1 lt 1 pt 7 lc rgb \"blue\" lw 4 " << endl;
// //   pl << "set style line 2  lt 1 lc rgb \"red\" pt 5 " << endl;
// //   pl << "set style line 3  lt 1 lc rgb \"green\" pt 5 " << endl;
// //   pl << "set style line 4  lt 1 lc rgb \"brown\" pt 5 " << endl;
// //   pl << "set style line 5 lt 1 pt 7 lc rgb \"black\" lw 4 " << endl; 
  
// //   pl << "set origin 0,0" << endl;
// //   pl << "set size 3,3 " << endl;
  
// //   pl << "unset border" << endl;
  
// //   pl << "set title \"" <<  fprefix << "\"" << endl;
// //   pl << "plot '" << fprefix << ".edges' w l ls 5 notitle ,\\" << endl;
// //   pl << "    \"<awk '{if($4 == \\\"a\\\")print}' " <<  fprefix << ".abd.vertices\" w p ls 2 notitle,\\" << endl;
// //   pl << "    \"<awk '{if($4 == \\\"b\\\")print}' " <<  fprefix << ".abd.vertices\" w p ls 3 notitle,\\" << endl;
// //   pl << "    \"<awk '{if($4 == \\\"d\\\")print}' " <<  fprefix << ".abd.vertices\" w p ls 4 notitle,\\" << endl;
// //   pl << "    \"<awk '{if($4 != \\\"d\\ && $4 != \\\"a\\ && $4 != \\\"b\\\")print}' " <<  fprefix << ".abd.vertices\" w p ls 5 notitle,\\" << endl;
// //   pl << "    \"" <<  fprefix << ".abd.vertices\" u 1:2:($3) w labels font \"Arial,25\" offset 2  notitle" << endl;
// //   pl.close();
// //   return; 
// // }
// // //*********************************************************************************************//
// // //*********************************************************************************************//
// // //*********************************************************************************************//



// // void PrintoutGnuplotExteriorPaths4x5(IsletGraph* I, int Gno, string fprefix, string directory, string type){
// //   bool verbose=0;
// //   cout << "Check " << directory << fprefix << ".pl for a figure created." << endl;  

// //   //Printout edges and vertices
// //   string sIsletNum = static_cast<ostringstream*>( &(ostringstream() << (*I).IsletNum) )->str();
// //   if(verbose)cout << "Printing out edges... ";
// //   PrintoutEdges(&(*I).adGraph,fprefix+".Islet"+sIsletNum,directory);
// //   bool append = 1;
// //   PrintoutEdges(&(*I).bGraph,fprefix+".Islet"+sIsletNum,directory,append);
// //   if(verbose)cout << "Finished" << endl;
// //   if(verbose)cout << "Printing out vertices... ";
// //   PrintoutExteriorPaths(I, fprefix+".Islet"+sIsletNum,directory);
// //   PrintoutVertices(&(*I).allGraph,fprefix+".Islet"+sIsletNum,directory);  
// //   if(verbose)cout << "Finished" << endl;
// //   if(verbose)cout << "Printing out exterior edges... ";
 
  
  
// //   //create pl file
// //   int fignum = Gno/20, pos = Gno%20;
// //   string sfignum = static_cast<ostringstream*>( &(ostringstream() << fignum) )->str();
  
// //   cout << "Printing out pl file... ";
// //   ofstream pl((directory+"figures/"+fprefix+".fig" + sfignum +".pl").c_str(),ios::app);
  
// //   if(pos==0){
// //     //Printout top matter
// //     pl << "set term post eps enhanced color \"Arial\" 26" << endl;
// //     pl << "set output '" << fprefix << ".fig" << sfignum << ".eps'" << endl;
// //     pl << "" << endl;
// //     pl << "set size 4,5" << endl;
// //     pl << "set multiplot" << endl;
// //     pl << "" << endl;
// //     pl << "set style line 1 lt 1 pt 7 lc rgb \"blue\" lw 4 " << endl;
// //     pl << "set style line 2  lt 1 lc rgb \"red\" pt 7 ps 0.7" << endl;
// //     pl << "set style line 3  lt 1 lc \"green\" pt 7 ps 0.7" << endl;
// //     //pl << "set style line 3  lt 1 lc \"#9932CC\" pt 7 ps 0.5" << endl; //indigo
// //     pl << "set style line 4  lt 1 lc rgb \"brown\" pt 7 ps 0.7" << endl;
// //     pl << "set style line 5 lt 1 pt 7 lc rgb \"orange\" lw 4 " << endl;
// //     pl << "set style line 6 lt 1 pt 7 lc rgb \"black\" lw 1" << endl;
// //     pl << "set style line 7 lt 1 pt 7 lc rgb \"black\" lw 4" << endl;
// //     pl << endl << endl;
// //   }
  
// //   pl << "set origin " << pos%4 << "," << 4-pos/4 << endl;
// //   pl << "set size 1,1 " << endl;
// //   pl << "unset border" << endl;
// //   pl << "set title \"Subject " << (*I).IsletNum/10000 -1 << ",Islet " << (*I).IsletNum %10000 <<": " << (*I).nalpha+(*I).nbeta+ (*I).ndelta << " cells\"" << endl;
  
// //   pl << "plot '" << fprefix <<".Islet" << sIsletNum <<  ".edges' w l ls 6 notitle ,\\" << endl;
// //   pl << "    '" << fprefix  << ".Islet" << sIsletNum<< ".exterioredges' w l ls 7 notitle,\\" << endl;
// //   pl << "    \"<awk '{if($4 == \\\"a\\\")print}' " <<  fprefix <<".Islet" << sIsletNum <<  ".abd.vertices\" w p ls 2 notitle,\\" << endl;
// //   pl << "    \"<awk '{if($4 == \\\"b\\\")print}' " <<  fprefix <<".Islet" << sIsletNum <<  ".abd.vertices\" w p ls 3 notitle,\\" << endl;
// //   pl << "    \"<awk '{if($4 == \\\"d\\\")print}' " <<  fprefix <<".Islet" << sIsletNum <<  ".abd.vertices\" w p ls 4 notitle,\\" << endl;
  
// //   pl << endl << endl;
// //   pl.close();
  
// //   return;
// // }
// // //*********************************************************************************************//
// // //*********************************************************************************************//
// // //*********************************************************************************************//



// // void PrintoutGnuplotIslets4x5(IsletGraph* I, int Gno, string fprefix, string directory, string type){
// //   bool verbose=0;
// //   cout << "Check " << directory << fprefix << ".pl for a figure created." << endl;  

// //   //Printout edges and vertices
// //   string sIsletNum = static_cast<ostringstream*>( &(ostringstream() << (*I).IsletNum) )->str();
// //   if(verbose)cout << "Printing out edges... ";
// //   PrintoutEdges(&(*I).adGraph,fprefix+".Islet"+sIsletNum,directory);
// //   bool append = 1;
// //   PrintoutEdges(&(*I).bGraph,fprefix+".Islet"+sIsletNum,directory,append);
// //   if(verbose)cout << "Finished" << endl;
// //   if(verbose)cout << "Printing out vertices... ";
// //   PrintoutVertices(&(*I).allGraph,fprefix+".Islet"+sIsletNum,directory);  
// //   if(verbose)cout << "Finished" << endl;
// //   if(verbose)cout << "Printing out exterior edges... ";
 
  
  
// //   //create pl file
// //   int fignum = Gno/20, pos = Gno%20;
// //   string sfignum = static_cast<ostringstream*>( &(ostringstream() << fignum) )->str();
  
// //   cout << "Printing out pl file... ";
// //   ofstream pl((directory+"figures/"+fprefix+".fig" + sfignum +".pl").c_str(),ios::app);
  
// //   if(pos==0){
// //     //Printout top matter
// //     pl << "set term post eps enhanced color \"Arial\" 26" << endl;
// //     pl << "set output '" << fprefix << ".fig" << sfignum << ".eps'" << endl;
// //     pl << "" << endl;
// //     pl << "set size 4,5" << endl;
// //     pl << "set multiplot" << endl;
// //     pl << "" << endl;
// //     pl << "set style line 1 lt 1 pt 7 lc rgb \"blue\" lw 4 " << endl;
// //     pl << "set style line 2  lt 1 lc rgb \"red\" pt 7 ps 0.7" << endl;
// //     pl << "set style line 3  lt 1 lc \"green\" pt 7 ps 0.7" << endl;
// //     //pl << "set style line 3  lt 1 lc \"#9932CC\" pt 7 ps 0.5" << endl; //indigo
// //     pl << "set style line 4  lt 1 lc rgb \"brown\" pt 7 ps 0.7" << endl;
// //     pl << "set style line 5 lt 1 pt 7 lc rgb \"orange\" lw 4 " << endl;
// //     pl << "set style line 6 lt 1 pt 7 lc rgb \"black\" lw 2" << endl;
// //     pl << endl << endl;
// //   }
  
// //   pl << "set origin " << pos%4 << "," << 4-pos/4 << endl;
// //   pl << "set size 1,1 " << endl;
// //   pl << "unset border" << endl;
// //   pl << "set title \"Subject " << (*I).IsletNum/10000 -1 << ",Islet " << (*I).IsletNum %10000 <<": " << (*I).nalpha+(*I).nbeta+ (*I).ndelta << " cells\"" << endl;
  
// //   pl << "plot '" << fprefix <<".Islet" << sIsletNum <<  ".edges' w l ls 6 notitle ,\\" << endl;
 
// //   pl << "    \"<awk '{if($4 == \\\"a\\\")print}' " <<  fprefix <<".Islet" << sIsletNum <<  ".abd.vertices\" w p ls 2 notitle,\\" << endl;
// //   pl << "    \"<awk '{if($4 == \\\"b\\\")print}' " <<  fprefix <<".Islet" << sIsletNum <<  ".abd.vertices\" w p ls 3 notitle,\\" << endl;
// //   pl << "    \"<awk '{if($4 == \\\"d\\\")print}' " <<  fprefix <<".Islet" << sIsletNum <<  ".abd.vertices\" w p ls 4 notitle,\\" << endl;
  
// //   pl << endl << endl;
// //   pl.close();
  
// //   return;
// // }



// // void PrintoutGnuplotHull3x4(IsletGraph* I, int Gno, string fprefix, string directory){
// //   bool verbose=0;
// //   cout << "Check " << directory << fprefix << ".pl for a figure created." << endl;  
// //   string type = "hull";
// //   //Printout edges and vertices
// //   string sIsletNum = static_cast<ostringstream*>( &(ostringstream() << (*I).IsletNum) )->str();
// //   string sGno = static_cast<ostringstream*>( &(ostringstream() << Gno) )->str();
// //   if(verbose)cout << "Printing out edges... ";
// //   //PrintoutEdges(&(*I).allGraph,fprefix+".Islet"+sIsletNum,directory);
// //   PrintoutVertices(&(*I).allGraph,fprefix+".Islet"+sIsletNum,directory);  
// //   if(verbose)cout << "Finished" << endl;cout.flush();
// //   if(verbose)cout << "Printing out hull edges... ";
// //   PrintoutSetLoops(I, fprefix+".Islet"+sIsletNum+"."+sGno,directory, type);
// //   if(verbose)cout << "Finished" << endl;cout.flush();
  
  
// //   //create pl file
// //   int fignum = Gno/12, pos = Gno%12;
// //   string sfignum = static_cast<ostringstream*>( &(ostringstream() << fignum) )->str();
  
// //   cout << "Printing out pl file... ";cout.flush();
// //   ofstream pl((directory+"figures/"+fprefix+".hull.fig" + sfignum +".pl").c_str(),ios::app);
  
// //   if(pos==0){
// //     //Printout top matter
// //     pl << "set term post eps enhanced color \"Arial\" 26" << endl;
// //     pl << "set output '" << fprefix << ".hull.fig" << sfignum << ".eps'" << endl;
// //     pl << "" << endl;
// //     pl << "set size 3,4" << endl;
// //     pl << "set multiplot" << endl;
// //     pl << "" << endl;
// //     pl << "set style line 1 lt 1 pt 7 lc rgb \"blue\" lw 4 " << endl;
// //     pl << "set style line 2  lt 1 lc rgb \"red\" pt 7 ps 0.7" << endl;
// //     pl << "set style line 3  lt 1 lc rgb \"green\" pt 7 ps 0.7" << endl;
// //     //pl << "set style line 3  lt 1 lc rgb \"#9932CC\" pt 7 ps 0.7" << endl;//indigo
// //     pl << "set style line 4  lt 1 lc rgb \"brown\" pt 7 ps 0.7" << endl;
// //     pl << "set style line 5 lt 1 pt 7 lc rgb \"orange\" lw 4 " << endl;
// //     pl << "set style line 6 lt 1 pt 7 lc rgb \"black\" lw 2" << endl;
// //     pl << endl << endl;
// //   }
  
// //   pl << "set origin " << pos%3 << "," << 3-pos/3 << endl;
// //   pl << "set size 1,1 " << endl;
// //   pl << "unset border" << endl;
// //   pl << "set title \"Islet " << (*I).IsletNum <<": t = " << Gno << "\"" << endl;
  
// //   // pl << "plot '" << fprefix <<".Islet" << sIsletNum <<  ".edges' w l ls 6 notitle ,\\" << endl;
// //   pl << "plot '" << fprefix <<".Islet" << sIsletNum <<  "." << sGno << ".hullsetLoops' w l ls 1 notitle,\\" << endl;
// //   pl << "    \"<awk '{if($4 == \\\"a\\\")print}' " <<  fprefix <<".Islet" << sIsletNum <<  ".abd.vertices\" w p ls 2 notitle,\\" << endl;
// //   pl << "    \"<awk '{if($4 == \\\"b\\\")print}' " <<  fprefix <<".Islet" << sIsletNum <<  ".abd.vertices\" w p ls 3 notitle,\\" << endl;
// //   pl << "    \"<awk '{if($4 == \\\"d\\\")print}' " <<  fprefix <<".Islet" << sIsletNum <<  ".abd.vertices\" w p ls 4 notitle,\\" << endl;
  
// //   pl << endl << endl;
// //   pl.close();
  
// //   return;
// // }

// // //*********************************************************************************************//
// // //*********************************************************************************************//
// // //*********************************************************************************************//

// // void PrintoutGnuplotPartialLoops(IsletGraph* I, string fprefix, string directory, string type){
// //   bool verbose=0;

// //   vector<vector<int> > * Ep,* Mp;
// //   Graph * Gp;
 
// //   if(type == "b"){
// //     Ep = &((*I).bExteriorPaths);
// //     Gp = &((*I).bGraph);
// //     Mp = &((*I).bMantle);
// //   }
// //   if(type == "ad"){
// //     Ep = &((*I).adExteriorPaths);
// //     Gp = &((*I).adGraph);
// //     Mp = &((*I).bMantle);
// //   }

// //   cout << "Check " << directory <<"figures/" << fprefix << ".pl for a figure created." << endl;  
// //   //create pl file
// //   ofstream pl((directory+"figures/"+fprefix+".pl").c_str());
// //   string sIsletNum = static_cast<ostringstream*>( &(ostringstream() << (*I).IsletNum) )->str();
// //   //Printout vertices
// //   if(verbose)cout << "Printing out vertices... ";
// //   PrintoutVertices(&(*I).allGraph,fprefix+".Islet"+sIsletNum,directory);
// //   if(verbose)cout << "Finished" << endl;
// // if(verbose)cout << "Printing out edges... ";
// //   PrintoutEdges(&(*I).adGraph,fprefix+".Islet"+sIsletNum,directory);
// //   bool append = 1;
// //   PrintoutEdges(&(*I).bGraph,fprefix+".Islet"+sIsletNum,directory,append);
// //   if(verbose)cout << "Finished" << endl;
// //   if(verbose)cout << "Printing out exterior edges... ";
// //   cout.flush();
// //   PrintoutExteriorPaths(I, fprefix+".Islet"+sIsletNum,directory);
// //   if(verbose)cout << "Finished" << endl;
// //   if(verbose)cout << "Printing out mantle... ";
// //   PrintoutMantlePaths(I,fprefix+".Islet"+sIsletNum,directory,type);
// //   if(verbose)cout << "Finished" << endl;
// //   if(verbose)cout << "Printing out sets... ";
// //   PrintoutSetLoops(I, fprefix+".Islet"+sIsletNum,directory, type); 
// //   if(verbose)cout << "Finished" << endl;

  
// //   //Printout top matter
// //   pl << "set term post eps enhanced color \"Arial\" 26" << endl;
// //   pl << "set output '" << fprefix << ".eps'" << endl;
// //   pl << "" << endl;
// //   pl << "set size 3,3" << endl;
// //   pl << "set multiplot" << endl;
// //   pl << "" << endl;
// //   pl << "set style line 1  lt 1 lc rgb \"red\" pt 7 ps 0.7" << endl; //alpha vertices
// //   pl << "set style line 2  lt 1 lc rgb \"green\" pt 7 ps 0.7" << endl; //beta vertices
// //   pl << "set style line 3  lt 1 lc rgb \"brown\" pt 7 ps 0.7" << endl; //delta vertices
// //   pl << "set style line 4 lt 1 pt 7 lc rgb \"black\" lw 2" << endl; //edges
// //   pl << "set style line 5 lt 1 pt 7 lc rgb \"black\" lw 3" << endl; //exterior edges
// //   pl << "set style line 6 lt 1 pt 7 lc rgb \"blue\" lw 5  ps 1 " << endl; //ad set loops
// //   pl << "set style line 7 lt 1 pt 7 lc rgb \"orange\" lw 5  ps 1 " << endl; //beta set loops
// //   pl << "set style line 8 lt 6 pt 7 lc rgb \"red\" lw 4  ps 1.5 " << endl; //ad set partial loops
// //   pl << "set style line 9 lt 9 pt 7 lc rgb \"green\" lw 4  ps 1.5 " << endl; //beta set partial loops
  
  
// //   pl << endl << endl;
  
// //   pl << "set origin 0,0" << endl;
// //   pl << "set size 3,3 " << endl;
// //   pl << "set style fill solid 1.0 noborder " << endl;
// //   // pl << "set xrange[" << xmin-10 << ":" << xmax+10 << "]" << endl;
// //   //     pl << "set yrange[" << ymin-10 << ":" << ymax+10 << "]" << endl;
// //   pl << "unset border" << endl;
// //   pl << "set title \"Islet " << (*I).IsletNum <<": " << (*I).nalpha+(*I).nbeta+ (*I).ndelta << " cells\"" << endl;
    
// //   pl << "plot '" << fprefix <<".Islet" << sIsletNum <<  ".edges' w l ls 4 notitle ,\\" << endl;
// //   pl << "    '"  << fprefix <<  ".exterioredges' w l ls 5 notitle,\\" << endl;
// //   pl << "    \"<awk '{if($4 == \\\"a\\\")print}' " <<  fprefix <<".Islet" << sIsletNum <<  ".abd.vertices\" w p ls 1 notitle,\\" << endl;
// //   pl << "    \"<awk '{if($4 == \\\"b\\\")print}' " <<  fprefix <<".Islet" << sIsletNum <<  ".abd.vertices\" w p ls 2 notitle,\\" << endl;
// //   pl << "    \"<awk '{if($4 == \\\"d\\\")print}' " <<  fprefix <<".Islet" << sIsletNum <<  ".abd.vertices\" w p ls 3 notitle,\\" << endl;
// //   pl << "    \"" <<  fprefix <<".Islet" << sIsletNum << ".abd.vertices\" u 1:2:($3) w labels font \"Arial,25\" offset 2  notitle,\\" << endl;
// //   pl << "    '" << fprefix <<".Islet" << sIsletNum <<  ".admantleedges' w l ls 9 notitle,\\" << endl;
// //   pl << "    '" << fprefix <<".Islet" << sIsletNum <<  ".bmantleedges' w l ls 8 notitle,\\" << endl;
// //   pl << "    '" << fprefix <<".Islet" << sIsletNum <<  ".adsetLoops' w l ls 6 notitle,\\" << endl;
// //   pl << "    '" << fprefix <<".Islet" << sIsletNum <<  ".betasetLoops' w l ls 7 notitle,\\" << endl;
// //   pl << "    '" << fprefix <<".Islet" << sIsletNum <<  ".bmantlevertices' w p ls 8 notitle,\\" << endl;
// //   pl << "    '" << fprefix <<".Islet" << sIsletNum <<  ".admantlevertices' w p ls 9 notitle" << endl;
  
// //   pl << endl << endl;
// //   pl.close();
  
// //   return;
// // }

// // //*********************************************************************************************//
// // //*********************************************************************************************//
// // //*********************************************************************************************//
// // void PrintoutGnuplotPartialLoops3x4(IsletGraph* I, int Gno, string fprefix, string directory, string type){
// //   bool verbose=0;
 
// // //   vector<vector<int> > * Ep;
// // //   Graph * Gp;
// // //   vector<vector<vertex_t> > * Mp;

// //   // if(type == "b"){
// // //     Ep = &((*I).bExteriorPaths);
// // //     Gp = &((*I).bGraph);
// // //     Mp = &((*I).bMantle);
// // //   }
// // //   if(type == "ad"){
// // //     Ep = &((*I).adExteriorPaths);
// // //     Gp = &((*I).adGraph);
// // //     Mp = &((*I).bMantle);
// // //   }

// //   cout << "Check " << directory <<"figures/" << fprefix << ".pl for a figure created." << endl;  
  
// //   string sIsletNum = static_cast<ostringstream*>( &(ostringstream() << (*I).IsletNum) )->str();
// //   //Printout vertices
// //   if(verbose)cout << "Printing out vertices... ";cout.flush();
// //   PrintoutVertices(&(*I).allGraph,fprefix+".Islet"+sIsletNum,directory);
// //   if(verbose)cout << "Finished" << endl;cout.flush();
// //   if(verbose)cout << "Printing out edges... ";cout.flush();
// //   PrintoutEdges(&(*I).adGraph,fprefix+".Islet"+sIsletNum,directory);
// //   bool append = 1;
// //   PrintoutEdges(&(*I).bGraph,fprefix+".Islet"+sIsletNum,directory,append);
// //   if(verbose)cout << "Finished" << endl;cout.flush();
// //   if(verbose)cout << "Printing out exterior edges... ";cout.flush();
// //   cout.flush();
// //   PrintoutExteriorPaths(I, fprefix+".Islet"+sIsletNum,directory);
// //   if(verbose)cout << "Finished" << endl;cout.flush();
// //   if(verbose)cout << "Printing out mantle... ";cout.flush();
// //   PrintoutMantlePaths(I,fprefix+".Islet"+sIsletNum,directory,type);
// //   if(verbose)cout << "Finished" << endl;
// //   if(verbose)cout << "Printing out sets... ";cout.flush();
// //   PrintoutSetLoops(I, fprefix+".Islet"+sIsletNum,directory, type); 
// //   if(verbose)cout << "Finished" << endl;cout.flush();
  
// //   //create pl file
// //   int fignum = Gno/12, pos = Gno%12;
// //   string sfignum = static_cast<ostringstream*>( &(ostringstream() << fignum) )->str();
  
// //   cout << "Printing out pl file... ";
// //   ofstream pl((directory+"figures/"+fprefix+".fig" + sfignum +".pl").c_str(),ios::app);
  
// //   if(pos==0){
// //     //Printout top matter
// //     pl << "set term post eps enhanced color \"Arial\" 26" << endl;
// //     pl << "set output '" << fprefix << ".fig" << fignum << ".eps'" << endl;
// //     pl << "" << endl;
// //     pl << "set size 3,4" << endl;
// //     pl << "set multiplot" << endl;
// //     pl << "" << endl;
// //     pl << "set style line 1  lt 1 lc rgb \"red\" pt 7 ps 0.7" << endl; //alpha vertices
// //     pl << "set style line 2  lt 1 lc rgb \"green\" pt 7 ps 0.7" << endl; //beta vertices
// //     pl << "set style line 3  lt 1 lc rgb \"brown\" pt 7 ps 0.7" << endl; //delta vertices
// //     pl << "set style line 4 lt 1 pt 7 lc rgb \"black\" lw 2" << endl; //edges
// //     pl << "set style line 5 lt 1 pt 7 lc rgb \"black\" lw 3" << endl; //exterior edges
// //     pl << "set style line 6 lt 1 pt 7 lc rgb \"blue\" lw 5  ps 1 " << endl; //ad set loops
// //     pl << "set style line 7 lt 1 pt 7 lc rgb \"orange\" lw 5  ps 1 " << endl; //beta set loops
// //     pl << "set style line 8 lt 6 pt 7 lc rgb \"red\" lw 4  ps 1.5 " << endl; //ad set partial loops
// //     pl << "set style line 9 lt 9 pt 7 lc rgb \"green\" lw 4  ps 1.5 " << endl; //beta set partial loops
// //   }
  
// //   pl << endl << endl;
// //   pl << "set origin " << pos%3 << "," << 3-pos/3 << endl;
// //   pl << "set size 1,1 " << endl;
// //   pl << "unset border" << endl;
// //   pl << "set title \"Islet " << (*I).IsletNum <<": " << (*I).nalpha+(*I).nbeta+ (*I).ndelta << " cells\"" << endl;
  
// //   pl << "plot '" << fprefix <<".Islet" << sIsletNum <<  ".edges' w l ls 4 notitle ,\\" << endl;
// //   pl << "    '"  << fprefix <<".Islet" << sIsletNum <<  ".exterioredges' w l ls 5 notitle,\\" << endl;
// //   pl << "    \"<awk '{if($4 == \\\"a\\\")print}' " <<  fprefix <<".Islet" << sIsletNum <<  ".abd.vertices\" w p ls 1 notitle,\\" << endl;
// //   pl << "    \"<awk '{if($4 == \\\"b\\\")print}' " <<  fprefix <<".Islet" << sIsletNum <<  ".abd.vertices\" w p ls 2 notitle,\\" << endl;
// //   pl << "    \"<awk '{if($4 == \\\"d\\\")print}' " <<  fprefix <<".Islet" << sIsletNum <<  ".abd.vertices\" w p ls 3 notitle,\\" << endl;
// //   // pl << "    \"" <<  fprefix <<".Islet" << sIsletNum << ".abd.vertices\" u 1:2:($3) w labels font \"Arial,25\" offset 2  notitle,\\" << endl;
// //    pl << "    '" << fprefix <<".Islet" << sIsletNum <<  ".admantleedges' w l ls 9 notitle,\\" << endl;
// //   pl << "    '" << fprefix <<".Islet" << sIsletNum <<  ".bmantleedges' w l ls 8 notitle,\\" << endl;
// //   pl << "    '" << fprefix <<".Islet" << sIsletNum <<  ".adsetLoops' w l ls 6 notitle,\\" << endl;
// //   pl << "    '" << fprefix <<".Islet" << sIsletNum <<  ".betasetLoops' w l ls 7 notitle,\\" << endl;
// //   pl << "    '" << fprefix <<".Islet" << sIsletNum <<  ".bmantlevertices' w p ls 8 notitle,\\" << endl;
// //   pl << "    '" << fprefix <<".Islet" << sIsletNum <<  ".admantlevertices' w p ls 9 notitle" << endl;
  
// //   pl << endl << endl;
// //   pl.close();
  
// //   return;
// // }

// // //*********************************************************************************************//
// // //*********************************************************************************************//
// // //*********************************************************************************************//

// // void PrintoutGnuplotPartialLoopsShadowRange(IsletGraph* I, int c, int pathnum, vector<double> * shadowrange, string fprefix, string directory, string type,string color){
// //   bool verbose=0;

// //   vector<vector<int> > * Ep;
// //   Graph * Gp;
   
// //   if(type == "b"){
// //     Ep = &((*I).bExteriorPaths);
// //     Gp = &((*I).bGraph);
// //   }
// //   if(type == "ad"){
// //     Ep = &((*I).adExteriorPaths);
// //     Gp = &((*I).adGraph);
// //   }

// //   cout << "Check " << directory <<"figures/" << fprefix << ".pl for a figure created." << endl;  
// //   //create pl file
// //   ofstream pl;
// //   string sIsletNum = static_cast<ostringstream*>( &(ostringstream() << (*I).IsletNum) )->str();
// //   if(!file_exists(directory+"figures/"+fprefix+".LoopShadowRange.pl")){
// //     pl.open((directory+"figures/"+fprefix+".LoopShadowRange.pl").c_str());
// //     //Printout vertices
// //     if(verbose)cout << "Printing out vertices... ";
// //     PrintoutVertices(&(*I).allGraph,fprefix+".Islet"+sIsletNum,directory);
// //     PrintoutEdges(&(*I).adGraph,fprefix+".Islet"+sIsletNum,directory);
// //     bool append = 1;
// //     PrintoutEdges(&(*I).bGraph,fprefix+".Islet"+sIsletNum,directory,append);
// //     if(verbose)cout << "Finished" << endl;cout.flush();
// //     if(verbose)cout << "Printing out exterior edges... ";cout.flush();
// //     PrintoutExteriorPaths(I, fprefix,directory);
// //     if(verbose)cout << "Finished" << endl;cout.flush();
// //     // //Determine x and y axes
// //     //     double xmin = 1000000, xmax = 0, ymin = 1000000, ymax = 0;
// //     //     for(int i = 0; i < num_vertices((*I).allGraph);i++){
// //     //       if((*I).allGraph[i].x < xmin)xmin = (*I).allGraph[i].x;
// //     //       if((*I).allGraph[i].x > xmax)xmax = (*I).allGraph[i].x;
// //     //       if((*I).allGraph[i].y < ymin)ymin = (*I).allGraph[i].y;
// //     //       if((*I).allGraph[i].y > ymax)ymax = (*I).allGraph[i].y;
// //     //     }
    
// //     //Printout top matter
// //     pl << "set term post eps enhanced color \"Arial\" 26" << endl;
// //     pl << "set output '" << fprefix << ".LoopShadowRange.eps'" << endl;
// //     pl << "" << endl;
// //     pl << "set size 3,3" << endl;
// //     pl << "" << endl;
// //     pl << "set style line 1 lt 1 pt 7 lc rgb \"blue\" lw 4  ps 1 " << endl;
// //     pl << "set style line 2  lt 1 lc rgb \"red\" pt 7 ps 0.1" << endl;
// //     pl << "set style line 3  lt 1 lc rgb \"green\" pt 7 ps 1" << endl;
// //     pl << "set style line 4  lt 1 lc rgb \"brown\" pt 7 ps 1" << endl;
// //     pl << "set style line 5 lt 1 pt 7 lc rgb \"black\" lw 4  ps 1 " << endl;
// //     pl << "set style line 6 lt 1 pt 7 lc rgb \"black\" lw 8" << endl;
// //     pl << "set style line 7 lt 3 pt 7 lc rgb \"purple\" lw 5  ps 1 " << endl;
// //     pl << endl << endl;
    
// //     pl << "set origin 0,0" << endl;
   
// //     pl << "set style fill solid 1.0 noborder " << endl;
// //     // pl << "set xrange[" << xmin-10 << ":" << xmax+10 << "]" << endl;
// //     //     pl << "set yrange[" << ymin-10 << ":" << ymax+10 << "]" << endl;
// //     pl << "unset border" << endl;
// //     pl << "set title \"Islet " << (*I).IsletNum <<": " << (*I).nalpha+(*I).nbeta+ (*I).ndelta << " cells\"" << endl;
// //   }
// //   else pl.open((directory+"figures/"+fprefix+".LoopShadowRange.pl").c_str(), ios::app);
// //   cout << "Shadowrange.size is " << (*shadowrange).size() << endl;
// //   for(int i = 0; i < (*shadowrange).size();i=i+2){
// //     cout << "i = " << i << endl;
    									    
// //     pl << "set obj ";
// //     if(color == "blue") pl <<c*1000+ i/2+pathnum*10 +1000001;
// //     else pl <<c*1000+ i/2+pathnum*10 +1;
// //     pl << " circle   arc [" << (*shadowrange)[i]/(2*PI)*360.0 << " : " << (*shadowrange)[i+1]/(2*PI)*360.0 << "] fc rgb \"" << color << "\" at  " << (*Gp)[(*Ep)[c][pathnum]].x << "," << (*Gp)[(*Ep)[c][pathnum]].y << " size  8  front" << endl;
// //   }
// //   //if(pathnum == 41){
// //   int endpos = (*Ep)[c].size()-2;
// //   cout << "endpos = " << endpos << endl;
// //   while ((*Ep)[c][endpos] >99999)endpos--;
// //   if(pathnum == endpos){
// //     cout << "pathnum = endpos" << endl;
// //     pl << "plot '" << fprefix <<".Islet" << sIsletNum <<  ".edges' w l ls 5 notitle ,\\" << endl;
// //     pl << "    '" << fprefix <<  ".exterioredges' w l ls 6 notitle,\\" << endl;
    
// //     pl << "    \"<awk '{if($4 == \\\"a\\\")print}' " <<  fprefix <<".Islet" << sIsletNum <<  ".abd.vertices\" w p ls 2 notitle,\\" << endl;
// //     pl << "    \"<awk '{if($4 == \\\"b\\\")print}' " <<  fprefix <<".Islet" << sIsletNum <<  ".abd.vertices\" w p ls 3 notitle,\\" << endl;
// //     pl << "    \"<awk '{if($4 == \\\"d\\\")print}' " <<  fprefix <<".Islet" << sIsletNum <<  ".abd.vertices\" w p ls 4 notitle,\\" << endl;
// //     //pl << "    \"" <<  fprefix <<".Islet" << sIsletNum << ".abd.vertices\" u 1:2:($3) w labels font \"Arial,25\" offset 2  notitle,\\" << endl;
// //     //pl << "    '" << fprefix <<  ".mantleedges' w l ls 7 notitle,\\" << endl;
// //     //pl << "    \"<awk '{if($3 == \\\"a\\\")print}' " <<  fprefix << ".mantlevertices\" w p ls 5 notitle,\\" << endl;
// //     //pl << "    \"<awk '{if($3 == \\\"b\\\")print}' " <<  fprefix << ".mantlevertices\" w p ls 1 notitle,\\" << endl;
// //     //pl << "    \"<awk '{if($3 == \\\"d\\\")print}' " <<  fprefix << ".mantlevertices\" w p ls 5 notitle" << endl;
    


// //   }
 
// //   pl << endl << endl;
// //   pl.close();
  
// //   return;
// // }

//*********************************************************************************************//
//*********************************************************************************************//
//*********************************************************************************************//

void PrintoutSubjectGofr(vector<vector<double> > *Gofr, string tfile){
  bool verbose=0;
  //string group = getGroup(fprefix);

  //Printout Gr file
  //string ssubjno = static_cast<ostringstream*>( &(ostringstream() << subjno) )->str();
  ////if(verbose)cout << "Printing out Gofr for " << fprefix << "." << ssubjno << endl;
  
  ofstream fout(tfile.c_str());
  for(int i = 0; i < (*Gofr).size();i++){
    for(int j = 0; j < (*Gofr)[i].size();j++)
      fout << (*Gofr)[i][j] << "\t";
    fout << endl;
  }
  fout.close();
  return;
}
//*********************************************************************************************//
//*********************************************************************************************//
//*********************************************************************************************//

// // void PrintoutGnuplotSubjectGofr(vector<vector<double> > *Gofr, int subjno, string fprefix, string directory, string type){
// //   bool verbose=0;
// // string ssubjno = static_cast<ostringstream*>( &(ostringstream() << subjno) )->str();
// //   cout << "Check " << directory << fprefix << ".pl for a figure created." << endl;  
// //   string group = getGroup(fprefix);
// //   //Printout Gr file
// //   PrintoutSubjectGofr(Gofr, subjno,fprefix,directory,type);
    
// //   //create pl file
// //   int fignum = (subjno-1)/12, pos = (subjno-1)%12;
// //   string sfignum = static_cast<ostringstream*>( &(ostringstream() << fignum) )->str();
  
// //   cout << "Printing out pl file " << directory << "figures/" << group << "." << type<<".gr.fig" << sfignum  << ".pl... ";
// //   ofstream pl((directory+"figures/"+group+"."+type+".gr.fig" + sfignum +".pl").c_str(), ios::app);
  
// //   if(pos==0){
// //     //Printout top matter
// //     pl << "set term post eps enhanced color \"Arial\" 26" << endl;
// //     pl << "set output '" << group << "." << type << ".gr.fig" << sfignum << ".eps'" << endl;
// //     pl << "" << endl;
// //     pl << "set size 3,4" << endl;
// //     pl << "set multiplot" << endl;
// //     pl << "" << endl;
// //     pl << "set style line 1 lt 1 pt 7 lc rgb \"black\" lw 8" << endl;
// //     pl << endl << endl;
// //   }
  
// //   pl << "set origin " << pos%3 << "," << 3-pos/3 << endl;
// //   pl << "set size 1,1 " << endl;
// //   pl << "set arrow from 0,1 to 50,1 nohead lt -1 lw 2 #creates line at y=1" << endl;
// //   pl << "set xtics font \"Arial-Bold,30\"" << endl;
// //   pl << "set ytics font \"Arial-Bold,30\"" << endl;
// //   pl << "unset border" << endl;
// //   pl << "set xrange [0:50]" << endl;
// //   pl << "set title \"Subject # " << subjno << endl;
  
// //   pl << "plot '" << fprefix <<".subj" << ssubjno << "." << type << ".gr' u 1:2  w l ls 1 notitle" << endl; 
// //   pl << endl << endl;
// //   pl.close();
  
// //   return;
// // }



// // //*********************************************************************************************//
// // //*********************************************************************************************//
// // //*********************************************************************************************//
// // void PrintoutMantlePaths(IsletGraph *I, string fprefix,string directory, string type){
// //   bool verbose=0;

// //   Graph * Gp;
// //   vector<vector<int> > * Mp;
// //   vector<int> *Mpt;
// //   if(verbose)cout << "Inside Printout Mantle Paths ..." << endl;
// //   if(verbose){
// //     cout << "Inside PrintoutMantlePaths: The ip's of I are " ;
// //     for(int ip = 0; ip < (*I).IPs.size();ip++)
// //       cout << (*I).IPs[ip].CellNum << "(" << (*I).IPs[ip].x << "," << (*I).IPs[ip].y << ") ";
// //     cout << endl;
// //   }

// //   if(type == "b" || type =="abd"){
// //     if(verbose) cout << "Printing out beta mantle" << endl;
// //     Mp = &(*I).bMantle;
// //     Mpt = &(*I).bMantle_type;
// //     Gp =&(*I).adGraph;
// //     ofstream gnupl((directory+"figures/" + fprefix+".bmantleedges").c_str());
// //     ofstream gnupv((directory+"figures/" + fprefix+".bmantlevertices").c_str());
// //     for(int c = 0; c < (*Mp).size(); c++){
// //       if(verbose)cout << " c = " << c  << endl;
// //       for(int j = 0; j < (*Mp)[c].size();j++){
// // 	if(verbose) cout << "  j = " << j << " and vertex is " <<(*Gp)[(*Mp)[c][j]].IsletCellNum <<  endl;
// // 	if((*Mp)[c][j] < 100000)gnupv << (*Gp)[(*Mp)[c][j]].x << "\t" << (*Gp)[(*Mp)[c][j]].y << "\t" << (*Gp)[(*Mp)[c][j]].type <<  endl;
// // 	//Check if edge between neighboring mantle
// // 	if(j < (*Mp)[c].size()-1){
// // 	  if((*Mp)[c][j]> 99999 && (*Mpt)[c]==1){
// // 	    //Find IP in IPlist
// // 	    int ippos = 0;
// // 	    for(int ip = 0; ip < (*I).IPs.size();ip++){
// // 	      if((*Mp)[c][j] == (*I).IPs[ip].CellNum){
// // 		ippos = ip;
// // 		break;
// // 	      }
// // 	    }
// // 	    if(verbose)cout << "Printing IP " << (*I).IPs[ippos].CellNum << "(" << (*I).IPs[ippos].x  << "," << (*I).IPs[ippos].y << ")" << endl;
// // 	    gnupl << (*I).IPs[ippos].x << "\t" << (*I).IPs[ippos].y << endl;
// // 	    if((*Mp)[c][j+1] > 99999){
// // 	      int ippos = 0;
// // 	      for(int ip = 0; ip < (*I).IPs.size();ip++){
// // 		if((*Mp)[c][j+1] == (*I).IPs[ip].CellNum){
// // 		  ippos = ip;
// // 		  break;
// // 		}
// // 	      }
// // 	      if(verbose)cout << "Printing IP " << (*I).IPs[ippos].CellNum << "(" << (*I).IPs[ippos].x  << "," << (*I).IPs[ippos].y << ")" << endl;
// // 	      gnupl << (*I).IPs[ippos].x << "\t" << (*I).IPs[ippos].y << endl;
// // 	    }
// // 	    else gnupl << (*Gp)[(*Mp)[c][j+1]].x << "\t" << (*Gp)[(*Mp)[c][j+1]].y << endl;
// // 	  }
// // 	  else if((*Mp)[c][j+1] > 99999 && (*Mpt)[c]==1){
// // 	    int ippos = 0;
// // 	    for(int ip = 0; ip < (*I).IPs.size();ip++){
// // 	      if((*Mp)[c][j+1] == (*I).IPs[ip].CellNum){
// // 		ippos = ip;
// // 		break;
// // 	      }
// // 	    }
// // 	    if(verbose)cout << "Printing IP " << (*I).IPs[ippos].CellNum << "(" << (*I).IPs[ippos].x  << "," << (*I).IPs[ippos].y << ")" << endl;
// // 	    gnupl << (*Gp)[(*Mp)[c][j]].x << "\t" << (*Gp)[(*Mp)[c][j]].y << endl;
// // 	    gnupl << (*I).IPs[ippos].x << "\t" << (*I).IPs[ippos].y << endl;
// // 	  }
// // 	  else if ((*Mp)[c][j] < 100000 && (*Mp)[c][j+1] < 100000){
// // 	    edge_t e; bool b;
// // 	    tie(e,b) = edge((*Mp)[c][j], (*Mp)[c][j+1], (*Gp));
// // 	    if(b){
// // 	      if(verbose)cout << "   edges exists" << endl;
// // 	      gnupl << (*Gp)[(*Mp)[c][j]].x << "\t" << (*Gp)[(*Mp)[c][j]].y << endl;
// // 	      gnupl << (*Gp)[(*Mp)[c][j+1]].x << "\t" << (*Gp)[(*Mp)[c][j+1]].y << endl;
// // 	      gnupl << endl;
// // 	    }
// // 	  }
// // 	}
// //       }
// //     }
// //     gnupl.close();
// //     gnupv.close();
// //   }
// //   if(type == "ad" || type =="abd"){
// //     if(verbose) cout << "Printing out ad mantle" << endl;
// //     Mp = &(*I).adMantle;
// //     Mpt = &(*I).adMantle_type;
// //     Gp = &(*I).bGraph;
// //     ofstream gnupl((directory+"figures/" + fprefix+".admantleedges").c_str());
// //     ofstream gnupv((directory+"figures/" + fprefix+".admantlevertices").c_str());
// //     for(int c = 0; c < (*Mp).size(); c++){
// //       for(int j = 0; j < (*Mp)[c].size();j++){
// // 	if(verbose){
// // 	  if((*Mp)[c][j] < 100000) cout << "  j = " << j << " and vertex is " <<(*Gp)[(*Mp)[c][j]].IsletCellNum <<  endl;
// // 	  else cout << "  j = " << j << " and vertex is " <<(*Mp)[c][j] <<  endl;
// // 	}
// // 	if((*Mp)[c][j] < 100000)gnupv << (*Gp)[(*Mp)[c][j]].x << "\t" << (*Gp)[(*Mp)[c][j]].y << "\t" << (*Gp)[(*Mp)[c][j]].type <<  endl;

// // 	//Check if edge between neighboring mantle
// // 	if(j < (*Mp)[c].size()-1){
// // 	  if((*Mp)[c][j]> 99999 && (*Mpt)[c]==1){
// // 	    //Find IP in IPlist
// // 	    int ippos = 0;
// // 	    for(int ip = 0; ip < (*I).IPs.size();ip++){
// // 	      if((*Mp)[c][j] == (*I).IPs[ip].CellNum){
// // 		ippos = ip;
// // 		break;
// // 	      }
// // 	    }
// // 	    if(verbose)cout << "Printing IP " << (*I).IPs[ippos].CellNum << "(" << (*I).IPs[ippos].x  << "," << (*I).IPs[ippos].y << ")" << endl;
// // 	    gnupl << (*I).IPs[ippos].x << "\t" << (*I).IPs[ippos].y << endl;
// // 	    if((*Mp)[c][j+1] > 99999){
// // 	      int ippos = 0;
// // 	      for(int ip = 0; ip < (*I).IPs.size();ip++){
// // 		if((*Mp)[c][j+1] == (*I).IPs[ip].CellNum){
// // 		  ippos = ip;
// // 		  break;
// // 		}
// // 	      }
// // 	      if(verbose)cout << "Printing IP " << (*I).IPs[ippos].CellNum << "(" << (*I).IPs[ippos].x  << "," << (*I).IPs[ippos].y << ")" << endl;
// // 	      gnupl << (*I).IPs[ippos].x << "\t" << (*I).IPs[ippos].y << endl;
// // 	    }
// // 	    else gnupl << (*Gp)[(*Mp)[c][j+1]].x << "\t" << (*Gp)[(*Mp)[c][j+1]].y << endl;
// // 	  }
// // 	  else if((*Mp)[c][j+1] > 99999 && (*Mpt)[c]==1){
// // 	    int ippos = 0;
// // 	    for(int ip = 0; ip < (*I).IPs.size();ip++){
// // 	      if((*Mp)[c][j+1] == (*I).IPs[ip].CellNum){
// // 		ippos = ip;
// // 		break;
// // 	      }
// // 	    }
// // 	    if(verbose)cout << "Printing IP " << (*I).IPs[ippos].CellNum << "(" << (*I).IPs[ippos].x  << "," << (*I).IPs[ippos].y << ")" << endl;
// // 	    gnupl << (*Gp)[(*Mp)[c][j]].x << "\t" << (*Gp)[(*Mp)[c][j]].y << endl;
// // 	    gnupl << (*I).IPs[ippos].x << "\t" << (*I).IPs[ippos].y << endl;
// // 	  }
// // 	  else if ((*Mp)[c][j] < 100000 && (*Mp)[c][j+1] < 100000){
// // 	    edge_t e; bool b;
// // 	    tie(e,b) = edge((*Mp)[c][j], (*Mp)[c][j+1], (*Gp));
// // 	    if(b){
// // 	      if(verbose)cout << "   edges exists" << endl;
// // 	      gnupl << (*Gp)[(*Mp)[c][j]].x << "\t" << (*Gp)[(*Mp)[c][j]].y << endl;
// // 	      gnupl << (*Gp)[(*Mp)[c][j+1]].x << "\t" << (*Gp)[(*Mp)[c][j+1]].y << endl;
// // 	      gnupl << endl;
// // 	    }
// // 	  }
// // 	}
// //       }
// //     }
// //     gnupl.close();
// //     gnupv.close();
// //   }
// //   return;
// // }






// // void PrintoutMeanAndStdDev(vector<vector<double> >* msd, string type, string fileprefix,string directory ){
// //   //Create file
// //   ofstream fout((directory+fileprefix+".grpeaks").c_str(), ios::app);
  
// //   for(int i = 0; i < 2;i++)
// //     fout << (*msd)[i][0] << "\t" << (*msd)[i][1] << "\t" ;
// //   if(type == "adb")fout << endl;
// //   fout.close();
  
// //   return ;
// // }
  


// // //*********************************************************************************************//
// // //*********************************************************************************************//
// // //*********************************************************************************************//


// // void PrintoutVerticesWithIPs(Graph* Gp, string fprefix, string directory,string celltype,bool append){
// //    std::pair<vertex_iter, vertex_iter> vp;
   
// //    ofstream gnupb;
// //    if(append)gnupb.open((directory+"figures/" + fprefix+"." + celltype+".vertices").c_str(),ios::app);
// //    else gnupb.open((directory+"figures/" + fprefix+"." + celltype+".vertices").c_str());

// //    //Printout points
// //    for(vp=vertices(*Gp); vp.first!=vp.second; vp.first++){
// //      if((celltype == "abd") ||
// // 	(((*Gp)[*vp.first].type == "a") && (celltype == "aa" || celltype == "ab"|| celltype == "ad" || celltype == "adm" || celltype == "adb")) ||
// // 	(((*Gp)[*vp.first].type == "b") && (celltype == "ab" || celltype == "bb"|| celltype == "bd" || celltype == "adb")) ||
// // 	(((*Gp)[*vp.first].type == "d") && (celltype == "ad" || celltype == "bd"|| celltype == "dd" || celltype == "adm" || celltype == "adb"))){
// //        gnupb << (*Gp)[*vp.first].x << "\t" << (*Gp)[*vp.first].y  << "\t"<< (*Gp)[*vp.first].IsletCellNum << "\t" << (*Gp)[*vp.first].type << "\n";
// //      }
// //    }
   
// //    gnupb.close();
// //    return; 
// // }


