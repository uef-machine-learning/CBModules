

#ifndef GRAPH_H
#define GRAPH_H

Graph *GraphRead_v2(char* FileName);
GraphVector *ReadGraphVectorFromFile_v2(FILE *f, Graph *g, int index);
void GraphWrite_v2(char* FileName, Graph *g,  int AllowOverWrite);
Graph* bruteForcekNNGraph(TRAININGSET* TS, int K, int mutual);
Graph* sampledkNNGraph(TRAININGSET* TS, int K, float sample);

double* dimBasedDensity(TRAININGSET* TS, int K, int nearestType);
double* knnGraphDensity(TRAININGSET* TS, int knnK, float sample);

#endif
