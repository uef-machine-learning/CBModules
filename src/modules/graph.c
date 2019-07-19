#ifndef _GRAPH_H_
#define _GRAPH_H_


#include <math.h>
#include <stdio.h>
#include <limits.h>
#include <float.h>

#include "cb.h"
#include "file.h"
#include "graph.h"
#include "knngraph.h"
#include "memctrl.h"
#include "random.h"

//Graph *ReadGraphFileHeader(FILE *f)
//{
  //int nvec, k, dim, i;
  //Graph *g;

  //fscanf(f, "Graph v. %i", &i);  // Read version
  //fscanf(f, "%i\n", &nvec);
  //fscanf(f, "%i\n", &k);
  //fscanf(f, "%i\n", &dim);

  //g = AllocateMemoryForGraph(nvec, k, dim);

  //while (getc(f) == '-' )             // Read dashes and '\n'
    //;

  //return g;
//}

/* ----------------------------------------------------------------- */

//GraphVector *AllocateMemoryForGraphVectorK(Graph *g, int index)
//{
  //GraphVector *v;

  //v = (GraphVector *) allocate(sizeof(GraphVector));
  //v->data = (int *) allocate(sizeof(int)*GraphGetDim(g));
  //v->kindices = (int *) allocate(sizeof(int)*GraphGetK(g));
  //v->needed = (int *) allocate(sizeof(int)*GraphGetK(g));
  //v->distances = (double *) allocate(sizeof(double)*GraphGetK(g));
  //v->index = index;
  //v->max_k =GraphGetK(g);
  //v->k = 0;

  //if (v == NULL || v->data == NULL || v->kindices == NULL || v->needed == NULL)
    //{
    //ErrorMessage("ERROR: Allocate error: AllocateMemoryForGraphVector()\n");
    //return NULL;
    //}

  //return v;
//}



Graph *GraphRead_v2(char* FileName)
{

    int nvec;
  int   i, j, min = INT_MAX, max = 0;
  FILE *f;
  Graph *g;
  int *data;

  f = FileOpen(FileName, INPUT, NO);
  fscanf(f, "%i\n", &nvec);
  printf("Reading %d graph nodes\n",nvec);
  //g = AllocateMemoryForGraph(nvec, 10 /*k*/, 10 /*dim*/);
  //TODO: allocate nodes dynamically.
  g = AllocateMemoryForGraph(nvec, 300 /*k*/, 10 /*dim*/);

  //g = ReadGraphFileHeader(f);
  if (g == NULL) {return NULL;}

  g->totalEdges=0;

  for (i=0; i< nvec; i++)
  {
      g->vectors[i] = ReadGraphVectorFromFile_v2(f, g, i);
  }

  fclose(f);

  //for (i=0; i< GraphGetNumberVectors(g); i++)
    //{
    //data = GraphGetVectorCoord(g->vectors[i]);
    //for (j = 0; j < GraphGetDim(g); j++)
      //{
      //if (data[j] < min)  min = data[j];
      //if (data[j] > max)  max = data[j];
      //}
    //}

  //g->maxcoord = max;
  //g->mincoord = min;

  return g;
}


int WriteGraphk_v2(FILE *f, Graph *g)
{

  fprintf(f, "Graph v. %d\n", 2);
  fprintf(f, "%d\n", g->nvec);
  //fprintf(f, "%d\n", g->k);
  fprintf(f, "%d\n", g->dim);

  fprintf(f, "-------------------------------------");
  putc(10, f);

  return 0;
}

void GraphWrite_v2(char* FileName, Graph *g,  int AllowOverWrite)
{
  int   i;
  FILE* f;

  f = FileOpen(FileName, OUTPUT, AllowOverWrite);
  //WriteGraphkFileHeader(f, g);

  fprintf(f, "%d\n", g->nvec);

  for(i=0; i<g->nvec; i++)
    {
    WriteGraphVectorToFile_v2(f, g, GraphGetVector(g, i));
    }

  fclose(f);
}


#define LINEBUF 10000
GraphVector *ReadGraphVectorFromFile_v2(FILE *f, Graph *g, int index)
{
  int i = 0, j = 0;
  int id,k,edgeTo;
  double dist;
  char *buffer, *temp;
  GraphVector *v;

  buffer = (char *) allocate(sizeof(char)*LINEBUF);

  temp = buffer;

  if (readline(f, buffer)) {
    free(temp);
    return 0;
  }
  //printf("ReadGraphVectorFromFile_v2\n");

  id = (int) strtol(buffer, &buffer, 10);
  k = (int) strtol(buffer, &buffer, 10);
  //printf("ID:%d ",id);
  //printf("%d ",k);
  v = GraphAddNode(g, NULL);
  for(i=0; i<k; i++) {
      edgeTo = (int) strtol(buffer, &buffer, 10);
      GraphAddEdge(g, index , edgeTo);
      g->totalEdges++;
      //printf ("%d ",edgeTo);
  }
  for(i=0; i<k; i++) {
      dist = (double) strtold(buffer, NULL);
      v->distances[i] = dist;
      //GraphAddEdge(g, index , edgeTo);
      //printf ("%f ",dist);
  }

  //printf("\n");



#ifdef DISABLED00
  v = AllocateMemoryForGraphVector(g, index);


	v->kindices[j++] = (int) strtol(buffer, &buffer, 10);

  while (*buffer) {
    while (*buffer && isspace(*buffer)) buffer++;
    if (*buffer) {
      if (i < GraphGetDim(g))
        {
	v->data[i++] = (int) strtol(buffer, &buffer, 10);
        }
      else
	v->kindices[j++] = (int) strtol(buffer, &buffer, 10);
      buffer++;
    }
  }
#endif

  free(temp);
  return v;
}



int WriteGraphVectorToFile_v2(FILE *f, Graph *g, GraphVector *v)
{
  int i;



  fprintf(f, "%d ",v->index);
  fprintf(f, "%d ",v->k);

  //for (i = 0; i < g->dim; i++)
    //{
    //fprintf(f, "%d ",v->data[i]);
    //}
  //for (i = 0; i < g->k; i++)
  for (i = 0; i < v->k; i++)
    {
    fprintf(f, "%d ",v->kindices[i]);
    }

  //TODO: if weighted
  for (i = 0; i < v->k; i++)
    {
        if(v->distances) {
            fprintf(f, "%f ",(float) (v->distances[i]));
        }
        else {
            fprintf(f, "0 ");
        }
    }

  fprintf(f, "\n");
  return 0;
}


static int cmpDoubleDesc(const void* a, const void* b, const void* info)
{
    int id1=*((int*)a);
    int id2=*((int*)b);
    return (((double*)info)[id1] > ((double*)info)[id2] ? 1 : 0);
}

static int cmpDoubleAsc(const void* a, const void* b, const void* info)
{
    int id1=*((int*)a);
    int id2=*((int*)b);
    return (((double*)info)[id1] < ((double*)info)[id2] ? 1 : 0);
}

double* kNNDensityForDim(TRAININGSET* TS, int knnK, int dim ) {

    kNNGraph* knng;
    Graph* g;
    int N;
    int edgeTo=0;
    N = BookSize(TS);
    int i,j;
    double eps = 1e-20;
    //int L=K-1;
    llong Dist=MAXLLONG;

    //knng = init_kNNGraph(N, knnK, knnK);
    double* density = calloc(N,sizeof(double)); // Init to zero
    double* dimvals = calloc(N,sizeof(double)); // Init to zero
    int* idx = malloc(sizeof(int)*N);
    //float* density = malloc(sizeof(float)*N);

    //printf("dim:%d knnK:%d\n",dim,knnK);

    for(i=0; i<N; i++)
    {
        idx[i]=i;
        dimvals[i] = (double) VectorScalar(TS,i,dim);
    }

    QuickSort(idx, N, sizeof(int), dimvals, cmpDoubleAsc);

    for( i=0; i<N; i++) {
        printf("dim:%d i:%d idx:%d dimval:%f\n",dim,i,idx[i], dimvals[idx[i]]);
    }

    for( i=0; i<N; i++)
    {
        int curid=idx[i];
        int idxLeft=i-1;
        int idxRight=i+1;
        int nearestId=0;
        double nearestVal=0;
        double nearestLeftVal=0;
        double nearestRightVal=0;
        double curNeighborVal=0;
        double curDens=0;
        for( j=0; j < knnK; j++)
        {
            nearestLeftVal=DBL_MAX;
            nearestRightVal=DBL_MAX;
            if(idxLeft >= 0) {
                nearestLeftVal=abs(dimvals[curid]-dimvals[idx[idxLeft]]);
            }
            if(idxRight < N) {
                nearestRightVal=abs(dimvals[curid]-dimvals[idx[idxRight]]);
            }
            if(nearestLeftVal < nearestRightVal) {
                curNeighborVal = nearestLeftVal;
                idxLeft--;
                //printf(" select left ");
            }
            else {
                idxRight++;
                //printf(" select right ");
                curNeighborVal = nearestRightVal;
            }
            curDens += curNeighborVal;
            //printf("dim:%d i:%d idRight:%d idLeft:%d nearestLeftVal:%f nearestRightVal:%f\n",dim,i,idxRight,idxLeft,nearestLeftVal,nearestRightVal);
        }

        //printf("XXdim:%d i:%d dist/density: %0.0f\n",dim,i,curDens);
        curDens = knnK/(curDens+eps);
        density[curid] = curDens;

        //if(i < 100 ) { printf("i=%d %f %f \n",i,dimvals[curid],density[curid]); }
    }

    return density;
    //printf("END Brute forcing\n");
}

double* kSlidingWindowDensityForDim(TRAININGSET* TS, int knnK, int dim ) {

    kNNGraph* knng;
    Graph* g;
    int N;
    int edgeTo=0;
    N = BookSize(TS);
    int i,j,windowStart,windowEnd;
    //int L=K-1;
    llong Dist=MAXLLONG;
    double eps = 1e-20;

    //knng = init_kNNGraph(N, knnK, knnK);
    double* density = calloc(N,sizeof(double)); // Init to zero
    double* dimvals = calloc(N,sizeof(double)); // Init to zero
    int* idx = malloc(sizeof(int)*N);
    int curid;
    double curval,curDens;
    //float* density = malloc(sizeof(float)*N);

    //printf("dim:%d knnK:%d\n",dim,knnK);

    for(i=0; i<N; i++)
    {
        idx[i]=i;
        dimvals[i] = (double) VectorScalar(TS,i,dim);
    }

    QuickSort(idx, N, sizeof(int), dimvals, cmpDoubleAsc);

    //for( i=0; i<N; i++) {
    //printf("dim:%d i:%d idx:%d dimval:%f\n",dim,i,idx[i], dimvals[idx[i]]);
    //}



    int kd = knnK/2;
    int leftRemove=0; //Index to remove from sum
    int rightRemove=0; //Index to remove from sum
    int leftAdd=0; //Index to remove from sum
    int rightAdd=0; //Index to remove from sum

    double sumLeft = 0;
    double sumRight = 0;
    for( i=1; i<=knnK; i++) {
        curid=idx[i];
        sumRight += dimvals[idx[curid]];
    }

    int windowLeftLength = 0;
    int windowRightLength = knnK;

    //TODO: Remove
    for( i=0; i<N; i++) {
        density[i] = 1.0;
    }
    for( i=0; i<N; i++)
    //for( i=kd; i<(N-kd); i++)
    {
        curid=idx[i];
        curval = dimvals[curid];
        curDens= abs(sumLeft - curval*((double)windowLeftLength)) + abs(sumRight - curval*((double)windowRightLength));
        curDens = knnK/(curDens+eps);
        density[curid] = curDens;

        leftRemove  = i - kd;
        rightRemove  = i+1;
        leftAdd = i;
        rightAdd = i+kd+1;
        //printf("i=%d curDens=%0.12f sumLeft=%f sumRight=%f leftAdd=%d leftRemove=%d rightAdd=%d rightRemove=%d windowLeftLength:%d windowRightLength:%d\n",i,curDens,sumLeft,sumRight,leftAdd,leftRemove,rightAdd,rightRemove,windowLeftLength,windowRightLength);
        //if(i >= kd && i < N-kd-1) {
            //printf("OK\n");
            //sumLeft = sumLeft - dimvals[idx[leftRemove]] + dimvals[idx[leftAdd]];
            //sumRight = sumRight - dimvals[idx[rightRemove]] + dimvals[idx[rightAdd]];
        //}
        if(i == N-1) {break;}

        if(leftRemove < 0 || rightAdd >= N) {
            sumLeft = sumLeft + dimvals[idx[leftAdd]];
            sumRight = sumRight - dimvals[idx[rightRemove]] ;
            windowLeftLength++;
            windowRightLength--;
            //printf("endpoint\n");
        } else {
            //printf("OK\n");
            sumLeft = sumLeft - dimvals[idx[leftRemove]] + dimvals[idx[leftAdd]];
            sumRight = sumRight - dimvals[idx[rightRemove]] + dimvals[idx[rightAdd]];

        }

    }

    return density;
    //printf("END Brute forcing\n");
}


double* kSlidingWindowDensityForDimOld(TRAININGSET* TS, int knnK, int dim ) {

    kNNGraph* knng;
    Graph* g;
    int N;
    int edgeTo=0;
    N = BookSize(TS);
    int i,j,windowStart,windowEnd;
    //int L=K-1;
    llong Dist=MAXLLONG;
    double eps = 1e-20;

    //knng = init_kNNGraph(N, knnK, knnK);
    double* density = calloc(N,sizeof(double)); // Init to zero
    double* dimvals = calloc(N,sizeof(double)); // Init to zero
    int* idx = malloc(sizeof(int)*N);
    //float* density = malloc(sizeof(float)*N);

    //printf("dim:%d knnK:%d\n",dim,knnK);

    for(i=0; i<N; i++)
    {
        idx[i]=i;
        dimvals[i] = (double) VectorScalar(TS,i,dim);
    }

    QuickSort(idx, N, sizeof(int), dimvals, cmpDoubleAsc);

    //for( i=0; i<N; i++) {
    //printf("dim:%d i:%d idx:%d dimval:%f\n",dim,i,idx[i], dimvals[idx[i]]);
    //}

    for( i=0; i<N; i++)
    {
        int curid=idx[i];
        double curDens=0;
        windowStart=i - knnK/2;
        if(windowStart < 0) { windowStart = 0; }
        windowEnd=windowStart+knnK;
        if(windowEnd >= N) { windowEnd = N-1; windowStart = windowEnd-(knnK);}
        //printf("windowStart:%d windowEnd:%d\n",windowStart,windowEnd);

        for( j=windowStart; j <= windowEnd; j++)
        {
            curDens+=abs(dimvals[curid]-dimvals[idx[j]]);
            //printf("dim:%d i:%d idRight:%d idLeft:%d nearestLeftVal:%f nearestRightVal:%f\n",dim,i,idxRight,idxLeft,nearestLeftVal,nearestRightVal);
        }

        //printf("XXdim:%d i:%d curid:%d dist/density: %0.0f\n",dim,i,curid,curDens);
        curDens = knnK/(curDens+eps);
        density[curid] = curDens;

        //if(i < 100 ) { printf("i=%d %f %f \n",i,dimvals[curid],density[curid]); }
    }

    return density;
    //printf("END Brute forcing\n");
}



double* kNNDensityForDimOld(TRAININGSET* TS, int knnK, int dim ) {

    kNNGraph* knng;
    Graph* g;
    int N;
    int edgeTo=0;
    N = BookSize(TS);
    int i,j;
    //int L=K-1;
    llong Dist=MAXLLONG;

    knng = init_kNNGraph(N, knnK, knnK);
    double* density = calloc(N,sizeof(double)); // Init to zero

    //printf("Estimating density for dim=%d\n",dim);
    for( i=0; i<N - 1; i++)
    {
        //if(i % 100 == 0 ) { printf("i=%d\n",i); }
        for( j=i+1; j<N; j++)
        {
            //Dist = abs((llong)v1[dim] - (llong)v2[dim]);
            Dist = abs(VectorScalar(TS,i,dim)-VectorScalar(TS,j,dim));

            updatekNN(knng, i, j, (float) Dist);
            updatekNN(knng, j, i, (float) Dist);
        }
    }

    for( i=0; i<N - 1; i++)
    {
        for ( j=0; j<knnK; j++)
        {
            density[i] += knng->list[i].items[j].dist/knnK;
        }
        density[i] = 1/density[i];
        //density[i] = g->vectors[i]->distances[knnK-1];
        //density[i] = (double) 1/knng->list[i].items[knnK-1].dist;
    }
    return density;
    //printf("END Brute forcing\n");
}

double* dimBasedDensity(TRAININGSET* TS, int K, int nearestType) {
    int i,j;
    int D = VectorSize(TS);
    int N = BookSize(TS);
    //double* density = malloc(sizeof(double)*N);
    double* density = calloc(N,sizeof(double)); // Init to zero
    double* densityForDim;
    for (i=0;i<D;i++)
    {
        if(nearestType==0)
        {
            densityForDim=kNNDensityForDim( TS,  K, i);
        }
        else if(nearestType==1) {
            densityForDim=kSlidingWindowDensityForDim( TS,  K, i);
        }
        else { assert(0);}

        for( j=0; j<N; j++)
        {
            density[j] += densityForDim[j];
        }
        free(densityForDim);
    }
    return density;
}

double* knnGraphDensity(TRAININGSET* TS, int knnK, float sample) {
    int i,j,N;
    N=BookSize(TS);
    //double* density = malloc(sizeof(double)*N);
    double* density = calloc(N,sizeof(double)); // Init to zero
    double eps = 1e-20;

    Graph* g;
    //printf("Initial centroids based on density\n");
    //printf("knnK=%d\n",knnK);
    if(sample > 0 && sample < 1) {
        g = sampledkNNGraph(TS, knnK,sample);
    }
    else {
        g = bruteForcekNNGraph(TS, knnK,0);
    }
    for(i=0; i<N; i++)
    {
        //density[i] = 1/(((double)(g->vectors[i]->distances[knnK-1]))+eps);
        density[i] = 1/(((double)(g->vectors[i]->distances[knnK-1]))+eps);
        for(j=0; j<knnK; j++) {
            density[i] += ((double)(g->vectors[i]->distances[j]));
        }
        density[i] = density[i]/knnK;
        if (density[i] < eps) { density[i] = eps; }
        //printf("i:%d density:%f\n",i,density[i]);
        density[i] = 1/density[i];

    }
    return density;
}



//TODO: mutual or not as parameter
Graph* bruteForcekNNGraph(TRAININGSET* TS, int K, int mutual) {
    kNNGraph* knng;
    Graph* g;
    int N;
    int edgeTo=0;
    N = BookSize(TS);
    int i,j;
    int L=K-1;
    llong Dist=MAXLLONG;

    knng = init_kNNGraph(N, K, K);

    printf("START Brute forcing\n");
    for( i=0; i<N - 1; i++)
    {
        //if(i % 100 == 0 ) { printf("i=%d\n",i); }
        for( j=i+1; j<N; j++)
        {
            Dist = sqrt(VectorDistance(Vector(TS, i),
                                       Vector(TS, j),
                                       VectorSize(TS),
                                       MAXLLONG,
                                       EUCLIDEANSQ));
            updatekNN(knng, i, j, (float) Dist);
            updatekNN(knng, j, i, (float) Dist);
        }
    }
    printf("END Brute forcing\n");

    //debug_graph(knng);

    //Maximum K*2 neighbors because kNN grpah may be converted to undirected graph.
    g = AllocateMemoryForGraph(N, K*2, VectorSize(TS));
    g->TS=TS;
    for( i=0; i<N; i++)
    {
        GraphAddNode(g, NULL);
        // Set data vector values for graph
        for(int i_dim=0; i_dim<VectorSize(TS); i_dim++)
        {
            g->vectors[i]->data[i_dim] = VectorScalar(TS,i,i_dim);
        }
    }

    for( i=0; i<N; i++)
    {
        for( j=0; j<K; j++)
        {
            double df = (double) knng->list[i].items[j].dist;
            //if(i==0) {printf("df:%f\n",df);}
            edgeTo = get_kNN_item_id(knng, i, j);
            if(mutual==1)
            {
                GraphAddMutualEdge(g,i,edgeTo);
            }
            else {
                GraphAddEdgeDist(g, i, edgeTo, df);
            }
        }
    }

    free_kNNGraph(knng);
    return g;
}

Graph* sampledkNNGraph(TRAININGSET* TS, int K, float sample) {
    kNNGraph* knng;
    Graph* g;
    int N;
    int edgeTo=0;
    N = BookSize(TS);
    int i,j,refPoint;
    int L=K-1;
    llong Dist=MAXLLONG;

    knng = init_kNNGraph(N, K, K);

    int numSamples = (int) (sample*N);

    //if (numSamples < K) { numSamples = K;}
    if (numSamples/2 < K) { K=numSamples/2;}
    printf("numSamples:%d K:%d\n",numSamples,K);
    int* samples = getRandomSampleInts(N-1, N-1);


    printf("START sampled kNN graph calculation sample:%f numSamples:%d\n",sample,numSamples);
    for( i=0; i<N - 1; i++)
    {
        //printf("%d ZZ\n",samples[i]);
        //if(i % 100 == 0 ) { printf("i=%d\n",i); }
        //for( j=i+1; j<N; j++)
        for( j=0; j<numSamples; j++)
        {

            refPoint = IRZ(N-1);
            //printf("%d AA\n",(i+j) % (N-1));
            //refPoint = samples[(i+j) % (N-1)];
            //printf("%d %d BB\n",(i+j) % N,refPoint);
            //int refPoint = j;
            // To prevent i == refPoint case:
            if(refPoint >= i) {refPoint++;}
            if(refPoint >= N) {refPoint=0;}
            Dist = sqrt(VectorDistance(Vector(TS, i),
                                       Vector(TS, refPoint),
                                       VectorSize(TS),
                                       MAXLLONG,
                                       EUCLIDEANSQ));
            updatekNN(knng, i, refPoint, (float) Dist);
            updatekNN(knng, refPoint, i, (float) Dist);
        }
    }
    printf("END Brute forcing\n");

    //debug_graph(knng);

    //Maximum K*2 neighbors because kNN grpah may be converted to undirected graph.
    g = AllocateMemoryForGraph(N, K*2, VectorSize(TS));
    g->TS=TS;
    for( i=0; i<N; i++)
    {
        GraphAddNode(g, NULL);
        // Set data vector values for graph
        for(int i_dim=0; i_dim<VectorSize(TS); i_dim++)
        {
            g->vectors[i]->data[i_dim] = VectorScalar(TS,i,i_dim);
        }
    }

    for( i=0; i<N; i++)
    {
        for( j=0; j<K; j++)
        {
            double df = (double) knng->list[i].items[j].dist;
            //if(i==0) {printf("df:%f\n",df);}
            edgeTo = get_kNN_item_id(knng, i, j);
            GraphAddEdgeDist(g, i, edgeTo, df);
        }
    }

    free_kNNGraph(knng);
    return g;
}


#endif

