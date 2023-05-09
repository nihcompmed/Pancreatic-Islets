#include <algorithm>   //needed for sort (vectors)
#include <boost/config.hpp>  //basic boost header file
#include <boost/graph/adjacency_list.hpp> //used for creating graph
#include <boost/graph/connected_components.hpp> 
#include <boost/graph/dijkstra_shortest_paths.hpp> 
#include <boost/graph/graph_traits.hpp> //needed for creating vertex_t
#include "CommonMath_VersionForPaper.h"
#include "FileIO_VersionForPaper.h"  //needed for stripping directory and suffix from filename
#include <string> //to use *.at()
#include <time.h>   //needed for run time calculation
#include <vector> 
#include <iomanip>

 using namespace std; 


struct Cell{
  double x, y;
  string type; //cell type i.e. a  = alpha, b = beta, d = delta
  int  CellNum; //the number that is incremented when a new cell is added i.e. ranges from 0 to N where N is the number of vertices in graph
  int IsletCellNum; //the cellnumber found in the original islet file
  /* int ncellcomp; //the number of cells in the component  */
  int CompNum; //the component number ranging from 0 to #of components-1  
  int vother;  //if vertex is in abdGraph and is a beta, then vother is the CellNum in the bGraph
  /* int componentlabel; //0=isolated,1=on interior, 2=on exterior edge 3= on interior edge, 4= on interior and exterior edge  */
  /* vector<int> IPparents; //for IP cells only; these are the vertex endpts of the edges that intersect *\/ */
  
  bool operator==(const Cell&other)const {return (x ==other.x && y == other.y);}
  bool operator!=(const Cell&other)const {return CellNum !=other.CellNum;}

};

struct CellEdge{
  double distance;
  /* double edge_weight;//used for shortest path algorithms */
  /* int label; //0=original;1=added (used when saturating graph); 2 = one parent is an IP */
   vector<int> adCompPathNums;//used for determining Partial mantles 
   vector<int> bCompPathNums; 
};


//Defines graph, vertex_t and edge_t 
typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS,Cell,CellEdge,boost::property<boost::vertex_color_t, boost::default_color_type,        boost::property<boost::vertex_degree_t, int,
  boost::property<boost::vertex_in_degree_t, int,
  boost::property<boost::vertex_out_degree_t, int> > > > > Graph;
typedef boost::graph_traits<Graph>::edge_descriptor edge_t;
typedef boost::graph_traits<Graph>::edge_iterator edge_iter;
typedef boost::graph_traits<Graph>::out_edge_iterator out_edge_iter;
typedef boost::graph_traits<Graph>::vertex_descriptor vertex_t; 
typedef boost::graph_traits<Graph>::vertex_iterator vertex_iter;

struct IsletGraph{
  Graph allGraph;
  Graph bGraph;
  Graph adGraph;
  int IsletNum;
  int nalpha;
  int nbeta;
  int ndelta;
  double IsletArea; 
  double admthresh;
  double bbthresh;
  double adbthresh;
  vector<vector<vertex_t> > allComponents,bComponents,adComponents; //vertices are wrt individual Graphs 
  vector<vector<vertex_t> > allExteriorPaths,bExteriorPaths,adExteriorPaths;//vertices are wrt individual Graphs 
  vector<vector<vertex_t> > betasets,adsets;//vertices are wrt individual Graphs 
  vector<int>betasetloopcomp, adsetloopcomp; 
  vector<vector<vertex_t> > betasetLoops,adsetLoops; 
  vector<vector<vertex_t> >bMantle; //ad cells surrounding each b-component; vertices are wrt ad Graph 
  vector<vector<vertex_t> >adMantle; 
  /* vector<Cell> IPs; */
   vector<int> bMantle_type, adMantle_type;//0 = unassigned; 1 = in set based on gr threshold; 2 = in set based on 2*gr threshold; 3 = in partial Mantle  */
   vector<vector<double> > bMantlePercent, adMantlePercent; 
  /* vector<int> bMantle_ncomp, adMantle_ncomp; */
};
/* struct ParentEdge{ */
/*   vertex_t v,ev1, ev2; */
/*   double distance; */
/* }; */

/* struct Path{ */
/*   vertex_t cell1,cell2; */
/*   double distance; */
/*   double totalDistance; */
/*   vector<vertex_t> path; */
/* }; */
  

struct VertexDistance{ //changed dist to distance
  vertex_t vertex;
  double distance;
};

struct CellAngle{
  vertex_t Cell;
  double Angle;
};


struct vertexListWeight{
  vertex_t vertex;
  int Parent;
  double Weight;
};

/* //AddEdgesVertices */
bool AddEdges(Graph* , string , float , double); 
/* bool AddEdgesSaturated(Graph* , string , double, int label = 0); */
/* void AddIPEdges(Graph* , string); */
/* void AddIPVertices(Graph* , string ); */
void AddVertices(vector<IsletGraph>* , string , string, string );  
 void CopyComponentEdgesFromGraph(IsletGraph* , int ,Graph * , string); 
/* void CopyComponentEdgesFromGraphWithIPs(Graph* , int ,Graph * ); */
void CopyComponentVerticesFromGraph(IsletGraph* , int ,Graph * , string, string votherType = ""); 
/* void CopyComponentVerticesFromGraphWithIPs(Graph* , int ,Graph *); */
 void CopyGraph(Graph*,Graph* ); 
 void CopyGraphEdges(Graph*, Graph *, string type = "abd"); 
/* void CopyGraphMinusIP(Graph*,Graph*); */
void CopyGraphVertices(Graph*, Graph *, string type = "abd", bool setvother = 0); 
void CopyIsletGraph(IsletGraph* , IsletGraph* ); 
/* void CopyOppositeComponents(IsletGraph* , Graph * ,string ,int , int , Graph * ); */
 void Create_adAnd_bGraphs(IsletGraph *); 
void InitializeEdge(Graph * , edge_t , double); 
/* void InitializeIP(Graph * , vertex_t ,IsletGraph *); */
 void InitializeIslet(IsletGraph * , Graph* , int , double); 
 void InitializeVertex(Graph * , vertex_t ,double ,double , int , string, int ,int vother = 99999); 
void InitializeVertexFromVertex(Graph * , Graph* , vertex_t , vertex_t, int, int vother = 99999 ); 

/* void RemoveInteriorVertices( Graph * ); */

//CompareBy
bool compareByAngle(const CellAngle &, const CellAngle &);
bool compareByFirstEntry(const vector<double> &, const vector<double> &); 
bool compareByGraphSizeBiggestFirst(const IsletGraph &, const IsletGraph &);
bool compareByGraphSizeSmallestFirst(const IsletGraph &, const IsletGraph &);
bool compareByIsletNumSmallestFirst(const IsletGraph &, const IsletGraph &);
bool compareByLastEntryBiggestFirst(const vector<double> &, const vector<double> &);
bool compareByWeight(const vertexListWeight &, const vertexListWeight &);
/* bool compareByAngleBiggestFirst(const CellAngle &, const CellAngle &); */
/* bool compareByCellNum(const Cell &, const Cell &); */
bool compareByDistance(const VertexDistance &, const VertexDistance &); 
/* bool compareByPathDistance(const Path &, const Path &); */
/* bool compareByPathtotalDistance(const Path &, const Path &); */

 bool compareBySecondEntry(const vector<int> &, const vector<int> &); 
 bool compareByLastEntry(const vector<int> &, const vector<int> &); 



bool compareByParent(const vertexListWeight &, const vertexListWeight &);

/* //Measures */
void AddIPtoGraph(IsletGraph * , Graph * ,double , double ,edge_t ,edge_t ); 
void CalcGofr(vector<IsletGraph>* ,vector<vector<double> > *, string ,int Gno=99999); //returns Gofr without printing  
/* void CalcGofrUsingShadow(vector<IsletGraph>* ,vector<vector<double> > *, string ,int); */
/* bool CleanIntersectionPointsInPaths(Graph *, vector<vertex_t> *, vector<vertex_t> *, bool isSetLoop = 0, vector<vector<double> > * set = NULL); */
int connected_components_DS(Graph *, vector<int> *); 
/* void CreateMaskOfIslet(Graph* ,vector<vector<bool> > *); */
void DetermineThreshold(vector<IsletGraph>* , string,bool );  
/* void DeterminePlanarGraph(IsletGraph * , Graph * , Graph *, string, string ); */
 void FindCellShadows(IsletGraph * , vertex_t , vertex_t , vector<double>*, double ); 
void FindExteriorPath(IsletGraph*,Graph* ,vertex_t , vector<vertex_t>* , vertex_t * , int); 
void FindInteriorPath(IsletGraph*,Graph* ,vertex_t , vector<vertex_t>* , vertex_t*  , int ,double ,double );  
 void FindNextCellInPath(IsletGraph* ,Graph*, vector<vertex_t>*, double*,int,vertex_t* ); 
 void FindPartialLoops( IsletGraph *,double,string); 
 void FindShadowRange(IsletGraph * ,  int , int , string ,double,  vector<double> *); 
void MakeGraphLocallyPlanar(IsletGraph * ,Graph * , vertex_t ); 
/* void MakeGraphPlanar(IsletGraph *,Graph *);//Isletgraph is also sent in order to keep IP numbering consistant  */
void MakeGraphPlanarByRemovingLongerEdges(Graph * ,bool opposite=0); 
/* double MinBetween15and35(vector<vector<double> > * ); */
/* double MinBetweenPeaks2and3(vector<vector<double> > *, vector<vector<double> > *); */
 double MinAfter2ndHighestPeak(vector<vector<double> > * , vector<vector<double> > * ); 
 void PeakFinder(vector<vector<double> > * , vector<vector<double> >*);  
/* double PolygonArea(Graph *,vector<vertex_t> *); */
/* void RemoveVerticesAwayFromSets(Graph * ,string, double); */
/* void RotateIsletMaximalDistanceXaxis(IsletGraph *, string); */
void SetComponentsAndLabelExteriorPaths(IsletGraph*,string, string,string sthresh = "gr");  

//bool setGrThresholds(vector<IsletGraph>* , string ,bool,bool , string gr_directory = "./"); 
bool setGrThresholds(vector<IsletGraph>* , bool , bool , string , string , string );


bool PathBetweenTwoVertices(Graph * , vertex_t , vertex_t , vector<vertex_t> * ); 
/* void RandomizeOppositeEdges(Graph *, string ); */
/* void RandomizeOppositeEdgesTest(IsletGraph *, string ); */
bool ShortestWeightedPathBetweenTwoVertices(Graph * , vertex_t , vertex_t , vector<vertex_t> *);
void SmoothVector5(vector<vector<double> > * , vector<vector<double> > * );  

/* //Printout */
void PrintoutEdges(Graph *,string,string directory = "./",bool append = 0);  
void PrintoutExteriorPaths(IsletGraph *, string ,string , string type="adb" ); 
void PrintoutGnuplot(Graph *, string , string directory="./",string celltype="abd"); 
/* void PrintoutGnuplotWithIPs(Graph * , string , string directory="./",string celltype="abd"); */
void PrintoutGnuplotExteriorPaths(IsletGraph * , string , string directory="./" , string type = "abd"); 

//void PrintoutGnuplotGofr4x5(vector<vector<double> > *, IsletGraph*,int,int, string , string, string,double );

void PrintoutGnuplotGofr4x5(vector<vector<double> > *, IsletGraph* , string , string , string , double );


/* void PrintoutGnuplotHull3x4(IsletGraph* , int , string, string directory= "./"); */
/* void PrintoutGnuplotIslet(vector<IsletGraph> * , int , string, string directory = "./",string option = "" ); */
void PrintoutGnuplotIslet(IsletGraph *  , string, string directory = "./",string option = "" ); 
/* void PrintoutGnuplotIslets4x5(IsletGraph* , int , string, string directory = "./", string type = "adb"); */
/* void PrintoutGnuplotExteriorPaths4x5(IsletGraph* , int , string, string directory = "./", string type = "adb"); */
/* void PrintoutGnuplotPartialLoops(IsletGraph* , string , string directory = "./", string type = "adb"); */
/* void PrintoutGnuplotPartialLoopsShadowRange(IsletGraph* , int , int , vector<double> * , string , string directory ="./", string type = "b",string color = "red" ); */
/* void PrintoutGnuplotPartialLoops3x4(IsletGraph* , int, string , string directory = "./", string type = "abd"); */
void PrintoutGnuplotRandomlyChooseIsletsBySizeandFraction(vector<IsletGraph> *I,string fileprefix);
void PrintoutGnuplotSetLoops(IsletGraph* , int , string, string directory ="./", string type = "adb"); 
/* void PrintoutGnuplotSubjectGofr(vector<vector<double> > *, int , string , string , string); */
/* void PrintoutIPVertices(Graph* , string); */
/* void PrintoutIPEdges(Graph* , string); */
/* void PrintoutLoopResults(vector<IsletGraph> *,string , string directory = "./" ); */
/* void PrintoutMantlePaths(IsletGraph * , string ,string , string ); */
/* void PrintoutMeanAndStdDev(vector<vector<double> >* , string , string ,string ); */
 void PrintoutPartialLoopResults(IsletGraph * , string , string directory = "./"); 
void PrintoutSetLoops(IsletGraph *, string ,string directory = "./", string type = "adb" ); 

//void PrintoutSubjectGofr(vector<vector<double> > * , int , string , string , string); 
void PrintoutSubjectGofr(vector<vector<double> > *, string );

/* void PrintoutSubjectSpecificInfo(vector<IsletGraph>* ,string,string ); */
void PrintoutVertices(Graph* , string ,string directory = "./", string celltype = "abd", bool append = 0); 
/* void PrintoutVerticesWithIPs(Graph* , string ,string directory = "./", string celltype = "abd", bool append = 0); */

//RemoveVectors.cpp
void ClearEdge(Graph * , edge_t );
 void ClearIsletGraphVectors(IsletGraph *); 
/* void ClearVertex(Graph * , vertex_t ); */
/* void RemoveCell(Cell *); */
void RemoveGraph(Graph *); 
void RemoveIsletGraph(IsletGraph*); 
/* void RemovePath(Path *); */
void RemoveVector(vector<CellAngle>*); 
void RemoveVector(vector<double>* ); 
void RemoveVector(vector<int>* );
/* void RemoveVector(vector<IsletGraph>*); */
/* void RemoveVector(vector<Path>* ); */
void RemoveVector(vector<string>* ); 
void RemoveVector(vector<VertexDistance>* ); 
 void RemoveVector(vector<vertexListWeight>*); 
void RemoveVector(vector<vertex_t>*); 
void RemoveVector(vector<edge_t>*); 
/* void RemoveVector(vector<Cell>*); */
void RemoveVector2d(vector<vector<double> > *); 
void RemoveVector2d(vector<vector<int> > * ); 
void RemoveVector2d(vector<vector<vertex_t> > * );
/* void RemoveVector3d(vector<vector<vector<double> > >* ); */
/* void RemoveVector3d(vector<vector<vector<int> > >* ); */

/* //SpanningTree.cpp */
void FindAD_BetaSets(IsletGraph*,string, string sthresh = "gr", string type = "all"); 
void FindCellSets(Graph*,Graph*,Graph *, vector<vector<vertex_t> > * ,vector<int>*, double); 
/* void FindPrelimCellSets(Graph* ,Graph*, Graph* ,vector<vector<vertex_t> > *); */
bool FindSetsAndLoops(IsletGraph* , string ,string, string sthresh = "gr"); 
 int FindSetLoops(IsletGraph* I,Graph*,Graph *, vector<vector<vertex_t> > *,vector<vector<vertex_t> >* ,vector<int>*, string,string sthresh = "gr"); 
double isLeft( double ,double , double ,double , double , double  );  

// EDITED BY MANU
 void SpanningTree(Graph*, Graph * );  

 int WindingNumber(Graph* , vector<vertex_t>* ,double, double  );  

