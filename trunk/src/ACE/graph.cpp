#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/kruskal_min_spanning_tree.hpp>



typedef struct {
  int nodePos;
  int nodeNeg;
  int stamp;
} Edges;

typedef std::vector<dataEdge *> Edges;

static Edges sEdges;

void initGraph(int nbEdge,int NbNodes){
  sEdges.reserve(nbEdge);
  ;
}

void graphAddEdge(int n1,int n2,void *ptr){
}

int computeGraphMST(){

}
