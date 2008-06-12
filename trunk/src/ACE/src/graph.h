#include "componentcap.h"

void stopGraph();
void printEdges();
int initGraph(int nbNodes, int nbEdges);
void graphAddEdge(int np,int nn,componentCAP *c);
int computeGraphMST();
int nextEdgeOutMST(int& np,int& nn,componentCAP **c);
int nextEdgeInMST(int& np,int& nn,componentCAP **c);
