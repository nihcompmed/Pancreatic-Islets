#include "GraphNN_VersionForPaper.h"

bool compareByAngle(const CellAngle &a, const CellAngle &b){
    return a.Angle < b.Angle;
}

bool compareByDistance(const VertexDistance &a, const VertexDistance &b){
    return a.distance < b.distance;
}

bool compareByFirstEntry(const vector<double> &a, const vector<double> &b){
  return a[0] < b[0];
}

bool compareByGraphSizeBiggestFirst(const IsletGraph &a, const IsletGraph &b){
  return num_vertices(a.allGraph) > num_vertices(b.allGraph);
}

 bool compareByGraphSizeSmallestFirst(const IsletGraph &a, const IsletGraph &b){
  return num_vertices(a.allGraph) < num_vertices(b.allGraph);
}

bool compareByIsletNumSmallestFirst(const IsletGraph &a, const IsletGraph &b){
  return a.IsletNum < b.IsletNum;
}

bool compareByLastEntryBiggestFirst(const vector<double> &a, const vector<double> &b){
  return a[a.size()-1] > b[b.size()-1];
}

bool compareByWeight(const vertexListWeight &a, const vertexListWeight &b){
    return a.Weight < b.Weight;
}

bool compareByParent(const vertexListWeight &a, const vertexListWeight &b){
    return a.Parent < b.Parent;
}

// bool compareByAngleBiggestFirst(const CellAngle &a, const CellAngle &b){
//     return a.Angle > b.Angle;
// }
// bool compareByCellNum(const Cell &a, const Cell &b){
//   return a.CellNum < b.CellNum;
// }


// bool compareByPathDistance(const Path &a, const Path &b){
//     return a.distance < b.distance;
// }
// bool compareByPathtotalDistance(const Path &a, const Path &b){
//     return a.totalDistance < b.totalDistance;
// }

 bool compareBySecondEntry(const vector<int> &a, const vector<int> &b){
   return a[1] < b[1];
 }
 bool compareByLastEntry(const vector<int> &a, const vector<int> &b){
   return a[a.size()-1] < b[b.size()-1];
 }





// // inline bool vertex_t::operator==(const vertex_t& lhs, const vertex_t& rhs){ 
// //   return lhs.CellNum == rhs.CellNum;
// // }
// // inline bool vertex_t::operator!=(const vertex_t& lhs, const vertex_t& rhs){
// //   return !operator==(lhs,rhs);
// // }


// // // bool compareByEvalue(const eig &a, const eig &b)
// // // {
// // //     return a.evalue < b.evalue;
// // // }

// // bool compareByCompSizeBiggestFirst(const vector<int> &a, const vector<int> &b)
// // {
// //   return a.size() > b.size();
// // }
// // bool compareByLastEntryBiggestFirst(const vector<double> &a, const vector<double> &b)
// // {
// //   return a[1] > b[1];
// // }

