#include "./GraphNN_VersionForPaper.h"
#define PI 3.14159

using namespace boost;

int main(int argc, char **argv){
  bool verbose=0;//used for debugging
  int vectorcount = 0;
  //std::cout << "Number of parameters " <<argc << endl;
  //*****************************************************************************//
  //                          Commandline Printout                               //
  //*****************************************************************************//
  if(argc<4){
    std::cout <<"Error: The number of parameters are not matched." << endl;
    std::cout <<"The command line should have:  " << endl;
    std::cout <<"./Simplex filename.tsv  AVtype option optionals" << endl;
    std::cout << "where AVtype (alpha,beta,delta) is:" << endl;
    std::cout << "\t rgc - rfp,gfp,cy5" << endl;
    std::cout << "\t 121113 - 12,11,13" << endl;
    std::cout << "\t rgcArea - if rgc and Islet area is included " << endl;
    std::cout << "\t 321Area - if 123 (a,b,d) and Islet area is included " << endl;
    std::cout << "\t abdArea - if 123 abd Islet area is included " << endl;
    std::cout << "where option is: " << endl;
    std::cout << "\t G - Only compute G(r) and not loops and partial loops on graphs" <<endl;
    std::cout << "\t B - Calculate loops and partial loops on graphs" <<endl;
    std::cout << "\t S - Subject specific information" << endl;
    std::cout << "and optionals: " << endl;
    std::cout << "\t -fout [filename] - out filename prefix " << endl;
    std::cout << "\t -dir [directory] - directory for results" << endl;
    std::cout << "\t -sort [sorttype]- sort by Islet size " << endl;
    std::cout << "\t\t sorttypes are: descending(default; biggest islet first) ascending(smallest islet first)" << endl;
    std::cout << "\t -runGr y - runs G(r) per islet" << endl;
    std::cout << "\t -grdir [directory] - directory for G(r) results if already computed" << endl;
    std::cout << "\t -figdir [directory] - directory for figures" << endl;
    std::cout << "\t -singleSubj [subject #] - only runs given subject number" << endl;
    std::cout << "\t -IsletNum [islet #] - only runs given islet" << endl;
    std::cout << "\t -print [printtype] - prints out figures" << endl;
    std::cout << "\t\t printtypes are: gr, graph" << endl;
    std::cout << endl;
    exit(0);
  }

  //*****************************************************************************//
  //                               Start Timer                                   //
  //*****************************************************************************//
  time_t time0 = time(0);
  time_t timeStart, timeEnd;
  time (&timeStart);

   
  //*****************************************************************************//
  //                        Set Command line variables                           //
  //*****************************************************************************//

  string file, AVtype, opt;
  int arg_count = 0;

  while(arg_count < argc){

    string temp = argv[arg_count++];
    if (temp.length() < 23){
      continue;
    }
    string temp_sub = temp.substr(temp.length() - 23);
    if(temp_sub == "Simplex_VersionForPaper"){
      //cout << "Found!" << endl;
      //getchar();
      file=argv[arg_count++];
      AVtype = argv[arg_count++];
      opt = argv[arg_count++];
      break;
    }

  }
  
  //cout << "filename::" << file << ", AVtype::" << AVtype << ", opt::" << opt << endl; 

  //Set tags
  string fileprefix, sorttype,directory,gr_directory,fig_directory;
  bool IsletSort=1,singleSubj = 0, printgr = 0,runGr = 0,printGraph = 0;
  int subjno = 0,isletNum = 0;

  while (arg_count < argc){

    string argvi = argv[arg_count++];

    //cout << "argument::" << argvi << endl;
    if(argvi.at(0) != '-'){
      cout << "Continuing because not a - flag" << endl;
      continue;
    }

    string argv1 = argv[arg_count++];

    //cout << "found a flag" << argvi << " " << argv1;

    if (argv1.at(0)=='-'){
      cout << "ERROR IN PARS" << endl;
      exit(1);
    }

    if(argvi == "-fout") fileprefix = argv1;
    else if(argvi == "-dir") directory= argv1;
    else if(argvi == "-grdir") gr_directory= argv1;
    else if(argvi == "-figdir") fig_directory= argv1;
    else if(argvi == "-sort") {
        IsletSort = 1;
	      sorttype = argv1;
    }
    else if(argvi == "-runGr") runGr = 1;
    else if(argvi == "-print"){
	      if(argv1=="gr") printgr = 1;
	      if(argv1=="graph") printGraph = 1;
    }

  }

  //cout << "filename::" << file << " AVtype::" << AVtype << " opt::" << opt << " runGr::" << runGr << endl; 

  //*****************************************************************************//
  //                       Create outputfile header                              //
  //*****************************************************************************//
 
 
  //string subjprefix = getFileName(getFileheader(file));

  //cout << "fileprefix contains " << fileprefix << endl;

  //if(fileprefix.size()==0){
  //  fileprefix = subjprefix;
  //  if(directory.size()==0)
  //    directory = getDirectory(file);
  //}
  //else if(directory.size() == 0)directory = getDirectory(fileprefix);
  //std::cout << "Analyzing file " <<file << endl;
  //std::cout << "Results are located in directory " << directory << endl;
  //std::cout << " With header " << fileprefix << endl;
  

  //*****************************************************************************//
  //                  Create Islet Graphs and add vertices                       //
  //*****************************************************************************//
  //cout << "Adding vertices to allGraph..." ;cout.flush();
  vector<IsletGraph> Islets;  vectorcount++;
  AddVertices(&Islets,file, "abd", AVtype);
  //cout << " Finished " << endl;cout.flush();
  //cout << "  There are " << Islets.size() <<" islets "<< endl;cout.flush();

  //  //*****************************************************************************//
  // //                  Printout subject-specific information                      //
  // //*****************************************************************************//
  // if(opt == "S"){
  //   PrintoutSubjectSpecificInfo(&Islets,fileprefix,subjprefix);
  //   return 0;
  // }

  //*****************************************************************************//
  //              Sort Islet Graphs by # of vertices (optional)                  //
  //*****************************************************************************//
  // if(IsletSort){
  //   cout << "sorting Islets..." ; 
  //   if(sorttype == "ascending")sort(Islets.begin(), Islets.end(), compareByGraphSizeSmallestFirst);
  //   else sort(Islets.begin(), Islets.end(), compareByGraphSizeBiggestFirst);
  //   cout << " Finished" << endl;
  // }

  //*****************************************************************************//
  //                         Add Edges to allGraph                               //
  //*****************************************************************************//
  double shadowthresh = 4.0;
  bool AddEdgesFailed = 0;
  if(gr_directory.size() == 0)gr_directory = "./Thresholds";
  if(fig_directory.size() == 0)fig_directory = "./figures";

  bool setGrFailed = setGrThresholds(&Islets,runGr,printgr,gr_directory,fig_directory,fileprefix); //from Thresholds folder or calculates it
  if(setGrFailed){
      string tfile = gr_directory + "/" + fileprefix + ".smooth2.thresh.minBetPeak2And3.dat";
      ofstream throut;
      throut.open(tfile.c_str());
      throut << "Gofr calculation failed";
      throut.close();
      return 1;
  }

  //exit(1);

  for(int i = 0; i < Islets.size();i++){
    AddEdgesFailed = AddEdges(&Islets[i].allGraph,"bb",Islets[i].bbthresh,shadowthresh);
    AddEdgesFailed = AddEdges(&Islets[i].allGraph,"adm",Islets[i].admthresh,shadowthresh);
    AddEdgesFailed = AddEdges(&Islets[i].allGraph,"adb",Islets[i].adbthresh,shadowthresh);
  }
      
  //PrintoutGnuplot(&Islets[0].allGraph,fileprefix+".grEdges");
  
  for(int i = 0;i < Islets.size();i++){
     //In allGraph-clean out intersections between adm and bb edges

    cout << "Making graph planar by removing longer edges" << endl;
    MakeGraphPlanarByRemovingLongerEdges(&Islets[i].allGraph,1);
    
    //Add vertices from allGraph to b- and adGraphs
    Create_adAnd_bGraphs(&Islets[i]);
    cout << "  Finished add edges for islet i = " << i << " and Islet Number " << Islets[i].IsletNum << ": " << endl;cout.flush();
    cout << "allGraphs: There are " << num_vertices(Islets[i].allGraph) << " vertices and " << num_edges(Islets[i].allGraph) << " edges" <<  endl;
    cout << "adGraphs: There are " << num_vertices(Islets[i].adGraph) << " vertices and " << num_edges(Islets[i].adGraph) << " edges" <<  endl;
    cout << "bbGraphs: There are " << num_vertices(Islets[i].bGraph) << " vertices and " << num_edges(Islets[i].bGraph) << " edges" <<  endl;
    //if(i < 400) PrintoutGnuplotIslets4x5(&Islets[i],i ,fileprefix+".Islets");
  }
  // if(runGr)return 0;
  // if(printGraph){
  //   int IsletVectorNum = 0;
  //   if(isletNum != 0){
  //     for(int i = 0; i < Islets.size();i++){
  // 	if(Islets[i].IsletNum == isletNum){
  // 	  IsletVectorNum = i;
  // 	  break;
  // 	}
  //     }
  //   }
  //   PrintoutGnuplotIslet(&Islets,IsletVectorNum,fileprefix);
  // }

  if(opt == "C") {
    PrintoutGnuplotRandomlyChooseIsletsBySizeandFraction(&Islets,fileprefix);
    return 0;
  }
  
  
  ////////////////////////////////////////////
  // MANU
  // At this point, we have the graphs

   // Write all vertices
  std::pair<vertex_iter, vertex_iter> vp;
  ofstream gnupb;

  string gnupbfile = gr_directory + "/" + "vertices.csv";
  //cout << "Saving " << gnupbfile;
   //string gnupbfile = "../VerticesEdges/" + fileprefix + ".vertices";

  gnupb.open(gnupbfile.c_str());

  Graph* Gp = &(Islets[0].allGraph);

   for(vp=vertices(*Gp); vp.first!=vp.second; vp.first++){
     if((*Gp)[*vp.first].CellNum < 99999){
       	 gnupb << (*Gp)[*vp.first].x << "\t" << (*Gp)[*vp.first].y  << "\t"<< (*Gp)[*vp.first].IsletCellNum << "\t" << (*Gp)[*vp.first].type << "\n";
     }
   }
   gnupb.close();

   // Add ad-edges and bb-edges
   //gnupbfile = "../VerticesEdges/" + fileprefix + ".edges";

   gnupbfile = gr_directory + "/" + "edges.csv";
   //cout << "Saving " << gnupbfile;


   gnupb.open(gnupbfile.c_str());

   std::pair<edge_iter, edge_iter> ep;

   // Add ad-edges
   Gp = &(Islets[0].adGraph);
   for(ep = edges((*Gp));ep.first!=ep.second;ep.first++){
     if((*Gp)[source(*ep.first,(*Gp))].IsletCellNum <(*Gp)[target(*ep.first,(*Gp))].IsletCellNum){
	       vertex_t v1 = source(*ep.first,(*Gp));
	       vertex_t v2 = target(*ep.first,(*Gp)); 
	       gnupb << (*Gp)[v1].IsletCellNum << "," << (*Gp)[v2].IsletCellNum <<endl;
     }
   }

    // Add bb-edges
   Gp = &(Islets[0].bGraph);
   for(ep = edges((*Gp));ep.first!=ep.second;ep.first++){
     if((*Gp)[source(*ep.first,(*Gp))].IsletCellNum <(*Gp)[target(*ep.first,(*Gp))].IsletCellNum){
	       vertex_t v1 = source(*ep.first,(*Gp));
	       vertex_t v2 = target(*ep.first,(*Gp)); 
	       gnupb << (*Gp)[v1].IsletCellNum << "," << (*Gp)[v2].IsletCellNum <<endl;
     }
   }

   gnupb.close();

  ////////////////////////////////////////////

   exit(0);

  // //*****************************************************************************//
  // //                       Option B: Calculate Loops                             //
  // //*****************************************************************************//

  vector<vector<vertex_t> > * Pp;

  if(opt == "B"){
    cout << "opt B: Calculating loops" << endl;
    for(int i = 0;i < Islets.size();i++){
       //if(verbose){
    	//PrintoutGnuplot(&Islets[i].allGraph,fileprefix+".allGraph");
    	//PrintoutGnuplot(&Islets[i].bGraph,fileprefix+".bGraph");
    	//PrintoutGnuplot(&Islets[i].adGraph,fileprefix+".adGraph");
       //}
       
       //cout << "Press key to set components and label exterior paths done." << endl;
       //getchar();

       //Set Components for each islet and label cells as isolated, interior, or exterior
       SetComponentsAndLabelExteriorPaths(&Islets[i],"b",fileprefix);
       SetComponentsAndLabelExteriorPaths(&Islets[i],"ad",fileprefix);

       //cout << "Set components and label exterior paths done. Press key to continue." << endl;
       //getchar();


  //     //if(i < 240)PrintoutGnuplotExteriorPaths4x5(&Islets[i],i,fileprefix+".exteriorpaths");
  //     //if(i < 240) PrintoutGnuplotIslets4x5(&Islets[i],i ,fileprefix+".Islets");
      
       //cout << "Label all beta cells...Press key..." << endl;
       //getchar();
       //Labels all beta cells that are located in loops
       // Default sthresh = 'gr' and type = 'all'
       
       FindAD_BetaSets(&Islets[i],fileprefix);

       ////////////////////////////////////////////
       // MANU
       // Print ad-loops around b-sets
       
       Gp = &(Islets[i]).adGraph;
       Pp = &(Islets[i]).betasetLoops;

       gnupb.open(("../Loops/" + fileprefix + ".adLoops").c_str());
       for(int i = 0; i < (*Pp).size();i++){
         for(int j = 0; j < (*Pp)[i].size();j++){
	           //gnupb << (*Gp)[(*Pp)[i][j]].IsletCellNum << ","  ;
             
	           gnupb << (*Gp)[(*Pp)[i][j]].x << "," << (*Gp)[(*Pp)[i][j]].y << "," ;
	           //gnupb << (*Gp)[(*Pp)[i][j+1]].x << "\t" << (*Gp)[(*Pp)[i][j+1]].y << endl;
	           //gnupb << endl;
         }
	       gnupb << endl;
       }
       gnupb.close();


       Gp = &(Islets[i]).bGraph;
       Pp = &(Islets[i]).adsetLoops;

       gnupb.open(("../Loops/" + fileprefix + ".bLoops").c_str());
       for(int i = 0; i < (*Pp).size();i++){
         for(int j = 0; j < (*Pp)[i].size();j++){
	           //gnupb << (*Gp)[(*Pp)[i][j]].IsletCellNum << ","  ;

	           gnupb << (*Gp)[(*Pp)[i][j]].x << "," << (*Gp)[(*Pp)[i][j]].y << "," ;
	           //gnupb << (*Gp)[(*Pp)[i][j+1]].x << "\t" << (*Gp)[(*Pp)[i][j+1]].y << endl;
	           //gnupb << endl;
         }
	       gnupb << endl;
       }
       gnupb.close();



       ////////////////////////////////////////////


       //if(i < 240)PrintoutGnuplotSetLoops(&Islets[i],i ,fileprefix+".setLoops");
      
  //     if(verbose) cout << "Finished ad and betasets" << endl;
  //     FindPartialLoops(&Islets[i],shadowthresh,fileprefix);
  //     //if(i< PrintoutGnuplotPartialLoops(&Islets[i],fileprefix+".PartialLoops");
  //     if(verbose)cout << "Before printout3x4" << endl;cout.flush();
  //     if(i < 240)PrintoutGnuplotPartialLoops3x4(&Islets[i],i,fileprefix+".PartialLoops");
  //     if(verbose)cout << "After printout3x4" << endl;cout.flush();
  //     PrintoutPartialLoopResults(&Islets[i], fileprefix);
  //   }
  //   cout << "finished option B" << endl;
    
    }
   
  // //Remove Islet vector
  // if(verbose)cout << "Before removing" << endl;cout.flush();
  // RemoveVector(&Islets); vectorcount--;
  // if(verbose)cout << "After removing" << endl;cout.flush();
  // if(vectorcount != 0){
  //   cout << "***Memory problem: make sure all vectors are deleted***" << endl;
  //   return 1;
  }

  //*****************************************************************************//
  //                           Stop Timer and Print                              //
  //*****************************************************************************//
  time (&timeEnd);
  double dtime = difftime(timeEnd, timeStart);
  double dmins = dtime/60; 
  int mins = static_cast<int>(dmins);
  int secs = static_cast<int>((dmins - mins)*60 + 0.01);
  double dhours = static_cast<double>(mins)/60;
  int hours = static_cast<int>(dhours);
  mins = static_cast<int>((dhours - hours)*60 + 0.01);
  double ddays = static_cast<double>(hours)/24;
  int days = static_cast<int>(ddays);
  hours = static_cast<int>((ddays - days)*24 + 0.01);
  
  std::cout << "Running time = " << days << " days " << hours << " hours " 
	    << mins << " minutes " << secs << " seconds\n"; 

  return 0;
} 


  

