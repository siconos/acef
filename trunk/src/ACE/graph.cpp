#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/kruskal_min_spanning_tree.hpp>
#include "graph.h"
#include <fstream>
#include <iostream>


using namespace boost;
  typedef adjacency_list < vecS, vecS, undirectedS,
    no_property, property < edge_weight_t, int > > Graph;
  typedef graph_traits < Graph >::edge_descriptor Edge;
  typedef graph_traits < Graph >::vertex_descriptor Vertex;
  typedef std::pair<int, int> E;
typedef std::list<E> listE;

typedef std::list<E>::iterator listInterE;


typedef std::map<E,componentCAP *> mapEdge;

static mapEdge sEdges;
static int sNbNodes=0;
static int sNbEdges=0;
static int sNb=0;
static E* spE=0;
static int *weights=0;
static listE sInE;
static listE sOutE;

static components sInComp;
static components sOutComp;
static componentsIt sInIt;
static componentsIt sOutIt;



int initGraph(int nbNodes, int nbEdges){
  sNbNodes=nbNodes;
  sNbEdges=nbEdges;
  sNb=0;
  spE=0;
  weights=0;
  if (sNbEdges==0){
    return 0;
  }
  spE = (E*) calloc(sNbEdges,sizeof(E));
  weights = (int*) calloc(sNbEdges,sizeof(int));
  return 1;
}

int nextEdgeInMST(int& np,int& nn,componentCAP **c){
  if (sNbEdges==0){
    return 0;
  }
  if (sInIt == sInComp.end())
    return 0;
  np = (*sInIt)->mNodePos;
  nn = (*sInIt)->mNodeNeg;
  (*c)=(componentCAP *)(*sInIt);
  sInIt++;
  return 1;
}
int nextEdgeOutMST(int& np,int& nn,componentCAP **c){
  if (sNbEdges==0){
    return 0;
  }
  if (sOutIt == sOutComp.end())
    return 0;
  np = (*sOutIt)->mNodePos;
  nn = (*sOutIt)->mNodeNeg;
  (*c)=(componentCAP *)(*sOutIt);
  sOutIt++;
  return 1;
}


void graphAddEdge(int n1,int n2,componentCAP *c){
  if (sNbEdges < sNb+1){
    ACE_INTERNAL_ERROR("graphAddEdge");
  }
  sEdges[E(n1,n2)]=c;
  spE[sNb]=E(n1,n2);
  sNb+=1;
}

int computeGraphMST(){
  if (!spE)
    return 0;
  Graph g(spE, spE + sNb, weights, sNbNodes);

  property_map < Graph, edge_weight_t >::type weight = get(edge_weight, g);
  std::vector < Edge > spanning_tree;

  kruskal_minimum_spanning_tree(g, std::back_inserter(spanning_tree));

  //edge in mst
    for (std::vector < Edge >::iterator ei = spanning_tree.begin();
       ei != spanning_tree.end(); ++ei) {
      sInE.push_back(E(source(*ei, g),target(*ei, g)));
  }

  graph_traits<Graph>::edge_iterator eiter, eiter_end;
  for (tie(eiter, eiter_end) = edges(g); eiter != eiter_end; ++eiter) {
    if (std::find(spanning_tree.begin(), spanning_tree.end(), *eiter)
        != spanning_tree.end())
      ;
    else
      sOutE.push_back(E(source(*eiter, g),target(*eiter, g)));
  }


  
  listInterE it;
  printf("print edge in:\n");
  for(it = sInE.begin(); it != sInE.end(); it++){
    
    if (sEdges.find(E(source(*it, g),target(*it, g))) != sEdges.end()){
      sInComp.push_back(sEdges[E(source(*it, g),target(*it, g))]);
      //std::cout<<"find edge "<<sEdges[E(source(*it, g),target(*it, g))]->mNodePos<<"\n";
    }
    std::cout<<source(*it, g)<<"  "<<target(*it, g)<<"\n";
  }
  std::cout<<"print edge out:\n";
  for(it = sOutE.begin();it != sOutE.end(); it++){
      if (sEdges.find(E(source(*it, g),target(*it, g))) != sEdges.end()){
      sOutComp.push_back(sEdges[E(source(*it, g),target(*it, g))]);
      //std::cout<<"find edge "<<sEdges[E(source(*it, g),target(*it, g))]->mNodePos<<"\n";
      }else{
	std::cout<<"cant find edge\n";
      }
      std::cout<<source(*it, g)<<"  "<<target(*it, g)<<"\n";
    }
  sInIt = sInComp.begin();
  sOutIt= sOutComp.begin();

  return 1;
  }
void stopGraph(){
  if (spE)
    free(spE);
  if(weights)
    free(weights);
  sInComp.clear();
  sOutComp.clear();
  sInE.clear();
  sOutE.clear();
  sEdges.clear();
}
void printEdges(){

  
}
