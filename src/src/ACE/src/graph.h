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

/*! \file graph.h

*/
#include "componentcap.h"

void stopGraph();
void printEdges();
int initGraph(int nbNodes, int nbEdges);
void graphAddEdge(int np,int nn,componentCAP *c);
int computeGraphMST();
int nextEdgeOutMST(int& np,int& nn,componentCAP **c);
int nextEdgeInMST(int& np,int& nn,componentCAP **c);
