/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2019 INRIA.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
*/

/*! \file graph.cpp

*/
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/kruskal_min_spanning_tree.hpp>
#include <boost/graph/visitors.hpp>
#include <boost/graph/breadth_first_search.hpp>
#include <boost/graph/dijkstra_shortest_paths.hpp>
#include <boost/graph/graph_utility.hpp>
#include "graph.h"
#include <fstream>
#include <iostream>

using namespace std;
using namespace boost;

  typedef adjacency_list < vecS, vecS, undirectedS,
    no_property, property < edge_weight_t, int > > Graph;
  typedef graph_traits < Graph >::edge_descriptor Edge;
  typedef graph_traits < Graph >::vertex_descriptor Vertex;
  typedef std::pair<int, int> E;
typedef std::list<E> listE;

typedef std::list<E>::iterator listInterE;




typedef property<vertex_color_t, default_color_type, 
         property<vertex_distance_t,int> > VProperty;
typedef int weight_t;
typedef property<edge_weight_t,weight_t> EProperty;

typedef adjacency_list<vecS, vecS, directedS, VProperty, EProperty > Graph_;

//static int sss_i=0;

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
static components sInBFSComp;
static components sOutComp;
static componentsIt sInIt;
static componentsIt sInItEnd;

static componentsIt sOutIt;

template < class Tag>
struct edge_printer
 : public boost::base_visitor< edge_printer< Tag> >
{
  typedef Tag event_filter;

  edge_printer()  { }

  template <class T, class Graph_>
  void operator()(T x, Graph_& g) {
    //printf("\t%d-%d-%d",sss_i, source(x, g),target(x, g));
    sInBFSComp.push_back(sEdges[E(source(x, g),target(x, g))]);
    //sss_i++;
  }
  
};
template < class Tag>
edge_printer< Tag>
print_edge( Tag) {
  return edge_printer< Tag>();
}


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
  if (sInIt == sInItEnd)
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
  ACE_CHECK_IERROR(sNbEdges >= sNb+1,"graphAddEdge, sNbEdges >= sNb+1 ");
  ACE_CHECK_IERROR(sEdges.find(E(n1,n2)) == sEdges.end(),"graphAddEdge, edge already in the graphe ");
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
	ACE_INTERNAL_WARNING("computeGraphMST: cant find edge.");
      }
      std::cout<<source(*it, g)<<"  "<<target(*it, g)<<"\n";
    }

  /*
  int N = sNbNodes;
  Graph_ G(N);
  boost::property_map<Graph_, vertex_index_t>::type 
    vertex_id = get(vertex_index, G);

  componentsIt cIt;
  
  for (cIt = sInComp.begin();cIt != sInComp.end(); cIt++){
    componentCAP *c;
    c=(componentCAP *)(*cIt);
    add_edge(c->mNodePos, c->mNodeNeg, 1, G);
  }

   
  ACE_MESSAGE("breadth_first_search: Starting graph\n");

  boost::breadth_first_search
    (G, vertex(0, G), 
     visitor(make_bfs_visitor(
           print_edge( on_examine_edge()))));

  int naux = sInBFSComp.size();
  int nauxbis = sInComp.size();
  if (sInBFSComp.size() != sInComp.size()){
  ACE_INTERNAL_WARNING("computeGraphMST BFS size!!");*/
    sInItEnd = sInComp.end();
    sInIt =  sInComp.begin();
    /* }else{
    sInItEnd = sInBFSComp.end();
    sInIt =  sInBFSComp.begin();
    }*/
  
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
  sInBFSComp.clear();
  sOutE.clear();
  sEdges.clear();
}
void printEdges(){

  
}
