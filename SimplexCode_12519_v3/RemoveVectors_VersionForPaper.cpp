#include "./GraphNN_VersionForPaper.h"

void ClearIsletGraphVectors(IsletGraph *I){
 //Clear out components
  RemoveVector2d(&((*I).allComponents));
  RemoveVector2d(&((*I).bComponents));
  RemoveVector2d(&((*I).adComponents));

    //Clear out Exterior Paths
  RemoveVector2d(&((*I).allExteriorPaths));
  RemoveVector2d(&((*I).bExteriorPaths));
  RemoveVector2d(&((*I).adExteriorPaths));

//    //Clear out sets
//   RemoveVector2d(&((*I).betasets));
//   RemoveVector2d(&((*I).adsets));

//   //Clear out setloops
//   RemoveVector2d(&((*I).betasetLoops));
//   RemoveVector2d(&((*I).adsetLoops));

//    //Clean out IPs
//   RemoveVector(&((*I).IPs));
 
  //Remove mantles
  RemoveVector2d(&((*I).bMantle));
  RemoveVector2d(&((*I).adMantle));

  //Remove mantle types
  RemoveVector(&((*I).bMantle_type));
  RemoveVector(&((*I).adMantle_type));

//   //Remove loop comps
//   RemoveVector(&((*I).betasetloopcomp));
//   RemoveVector(&((*I).adsetloopcomp));

  //Remove mantle percents
  RemoveVector2d(&((*I).bMantlePercent));
  RemoveVector2d(&((*I).adMantlePercent));
  
//   //Remove mantle numcomps
//   RemoveVector(&((*I).bMantle_ncomp));
//   RemoveVector(&((*I).adMantle_ncomp));

   return;
}

//*********************************************************************************************//
//*********************************************************************************************//
//*********************************************************************************************//

void RemoveGraph(Graph * Graph){
  std::pair<edge_iter, edge_iter> ep; 
  std::pair<vertex_iter, vertex_iter> vp; 
  
  // for(ep = edges(*Graph);ep.first != ep.second; ep.first++){
  //   (*Graph)[*ep.first].adCompPathNums.clear();vector<int> ().swap((*Graph)[*ep.first].adCompPathNums);
  //   (*Graph)[*ep.first].bCompPathNums.clear();vector<int> ().swap((*Graph)[*ep.first].bCompPathNums);
  // }
  
  (*Graph).clear();
  
  return;
}

//*********************************************************************************************//
//*********************************************************************************************//
//*********************************************************************************************//

void RemoveIsletGraph(IsletGraph* I){
  ClearIsletGraphVectors(I);
  
  //Remove Graphs
  RemoveGraph(&(*I).allGraph);
  RemoveGraph(&(*I).adGraph);
  RemoveGraph(&(*I).bGraph);
   
  return;
}

//*********************************************************************************************//
//*********************************************************************************************//
//*********************************************************************************************//

void ClearEdge(Graph * g, edge_t e){
 
  return;
}

// void ClearVertex(Graph * g, vertex_t v){
//   std::pair<out_edge_iter, out_edge_iter> oep;
//   for(oep = out_edges(v,*g);oep.first!=oep.second;oep.first++)ClearEdge(g,*oep.first);
//   (*g)[v].IPparents.clear();vector<int> ().swap((*g)[v].IPparents);
 
//   return;
// }

// void RemoveCell(Cell * cell){
//   (*cell).IPparents.clear();vector<int> ().swap((*cell).IPparents);
//   return;
// }






// void RemovePath(Path *P){
//   (*P).path.clear();vector<vertex_t> ().swap((*P).path);
//   return;
// }

// void RemoveVector(vector<IsletGraph> * I){
//     for(int i = 0; i < (*I).size();i++)
//     RemoveIsletGraph(&(*I)[i]);
//     (*I).clear();vector<IsletGraph> ().swap(*I);
//     return;
// }

 void RemoveVector(vector<vertexListWeight>* vec){
  (*vec).clear();vector<vertexListWeight> ().swap(*vec);
  return;
}

void RemoveVector(vector<int>* vec){
  (*vec).clear();vector<int> ().swap(*vec);
  return;
}
 void RemoveVector(vector<string>* vec){
   (*vec).clear();vector<string> ().swap(*vec);
   return;
 }

// void RemoveVector(vector<Path>* vec){
//   for(int i = 0; i < (*vec).size();i++)
//     RemovePath(&(*vec)[i]);
//   (*vec).clear();vector<Path> ().swap(*vec);
//   return;
// }
void RemoveVector(vector<double>* vec){
  (*vec).clear();vector<double> ().swap(*vec);
  return;
}

void RemoveVector(vector<VertexDistance>* vec){
  (*vec).clear();vector<VertexDistance> ().swap(*vec);
  return;
}
void RemoveVector(vector<CellAngle>* vec){
  (*vec).clear();vector<CellAngle> ().swap(*vec);
  return;
}

void RemoveVector(vector<vertex_t>* vec){
  (*vec).clear();vector<vertex_t> ().swap(*vec);
  return;
}

void RemoveVector(vector<edge_t>* vec){
  (*vec).clear();vector<edge_t> ().swap(*vec);
  return;
}
// void RemoveVector(vector<Cell>* vec){
//   for(int i = 0; i < (*vec).size();i++)RemoveCell(&((*vec)[i]));
//   (*vec).clear();vector<Cell> ().swap(*vec);
//   return;
// }
void RemoveVector2d(vector<vector<int> > * vec){
   for(int i = 0; i < (*vec).size();i++){
     (*vec)[i].clear();    vector<int> ().swap((*vec)[i]);
   }
  (*vec).clear(); vector<vector<int> > ().swap(*vec);
  return;
}
void RemoveVector2d(vector<vector<double> > * vec){
  for(int i = 0; i < (*vec).size();i++){
    (*vec)[i].clear();vector<double> ().swap((*vec)[i]);
  }
  (*vec).clear(); vector<vector<double> > ().swap(*vec);
  return;
}
void RemoveVector2d(vector<vector<vertex_t> > * vec){
  for(int i = 0; i < (*vec).size();i++){
    (*vec)[i].clear();vector<vertex_t> ().swap((*vec)[i]);
  }
  (*vec).clear(); vector<vector<vertex_t> > ().swap(*vec);
  return;
}
// void RemoveVector3d(vector<vector<vector<double> > >* vec){
//   for(int i = 0; i < (*vec).size();i++){
//     for(int j = 0; j < (*vec)[i].size();j++){
//       (*vec)[i][j].clear(); vector<double> ().swap((*vec)[i][j]);
//     }
//     (*vec)[i].clear();vector<vector<double> > ().swap((*vec)[i]);
//   }
//   (*vec).clear();vector<vector<vector<double> > > ().swap(*vec);
//   return;
// }

// void RemoveVector3d(vector<vector<vector<int> > >* vec){
//   for(int i = 0; i < (*vec).size();i++){
//     for(int j = 0; j < (*vec)[i].size();j++){
//       (*vec)[i][j].clear(); vector<int> ().swap((*vec)[i][j]);
//     }
//     (*vec)[i].clear();vector<vector<int> > ().swap((*vec)[i]);
//   }
//   (*vec).clear();vector<vector<vector<int> > > ().swap(*vec);
//   return;
// }
