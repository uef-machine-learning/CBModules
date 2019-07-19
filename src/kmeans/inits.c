/*-------------------------------------------------------------------*/
/* inits.c         Sami Sieranoja                                    */
/*                                                                   */
/* Different K-means initializations                                 */
/*-------------------------------------------------------------------*/

#include <float.h>
#include <stddef.h> //NULL definition
#include <stdio.h>
#include <math.h>

#include "cb.h"
#include "graph.h"
#include "sort.h"

#include "inits.h"
#include "kmeans.h"


#include "random.h"
//#include "interfc.h"
//#include "reporting.h"
//#include "limits.h"

#include "memctrl.h"
#include "reporting.h"




void SelectRandomRepresentatives(TRAININGSET *pTS, CODEBOOK *pCB)
{

    int k, n, x, Unique;

    for (k = 0; k < BookSize(pCB); k++)
    {
        do
        {
            Unique = 1;
            x = irand(0, BookSize(pTS) - 1);
            for (n = 0; (n < k) && Unique; n++)
                Unique = !EqualVectors(Vector(pTS, x), Vector(pCB, n), VectorSize(pCB));
        }
        while (!Unique);

        CopyVector(Vector(pTS, x), Vector(pCB, k), VectorSize(pCB));
        VectorFreq(pCB, k) = 0;
    }

}


/*-------------------------------------------------------------------*/
/* Initializes codebook with randomly selected vectors from dataset. */
/* Duplicate vectors are allowed.                                    */


void SelectRandomRepresentatives2(TRAININGSET *pTS, CODEBOOK *pCB) /* mm */
{
    int i,j;

    TRAININGSET tempTS;
    CreateNewCodebook(&tempTS,BookSize(pTS),pTS);
    CopyCodebook(pTS,&tempTS);
    for (i=0; i<BookSize(pCB); i++)
    {
        j = IRI(i,BookSize(&tempTS));
        CopyVector(Vector(&tempTS,j),Vector(pCB,i),VectorSize(&tempTS));
        VectorFreq(pCB,i)=1;
        CopyVector(Vector(&tempTS,i),Vector(&tempTS,j),VectorSize(&tempTS));
    }
    FreeCodebook(&tempTS);
}


/*-------------------------------------------------------------------*/
/* Pasi's solution 20.9.2016                                         */
/* (1) Shuffle training set. (2) Select first k vectors.             */
/* Someone else please test this...                                  */
/*-------------------------------------------------------------------*/


void RandomCodebook(TRAININGSET *pTS, CODEBOOK *pCB)
{
    int i;

    ShuffleTS(pTS);
    for(i=0; i<BookSize(pCB); i++)
    {
        CopyNode( &Node(pTS,i), &Node(pCB,i), VectorSize(pCB));
    }
}


/*-------------------------------------------------------------------*/
/* Marko's variant: uniform sampling. Not random if TS sorted!       */
/*-------------------------------------------------------------------*/


void SelectRandomRepresentativesbyMarko(TRAININGSET *pTS, CODEBOOK *pCB)
{
    int i, j, k;

    k = BookSize(pTS) / BookSize(pCB);

    for (i = 0; i < BookSize(pCB); i++)
    {
        /* interval [0,BookSize(pTS)) is divided M subintervals and
           random number is chosen from every subinterval */
        if (i == (BookSize(pCB) - 1))
        {
            /* random number generator must be initialized! */
            j = IRI(i*k, BookSize(pTS));
        }
        else
        {
            j = IRI(i*k, (i+1)*k);
        }

        CopyVector(Vector(pTS, j), Vector(pCB, i), VectorSize(pTS));
        VectorFreq(pCB, i) = 1;
    }

}


/*-------------------------------------------------------------------*/
/*  Simple vector by Ville  */
/*-------------------------------------------------------------------*/


llong* GetLongVector(int dim)
{
    llong *ptr;

    ptr = (llong *) malloc(sizeof(llong)*dim);

    if (ptr == NULL)
    {
        printf("Out of memory: in GetLongVector()\n");
        exit(-1);
    }

    return ptr;
}


/*-------------------------------------------------------------------*/
/* Random selection of type k-means++:                               */
/* The Advantages of Careful Seeding by Arthur and Vassilvitskii     */
/* Code by Ville Hautamäki (following original C++ implementation)   */
/*-------------------------------------------------------------------*/


void SelectRandomWeightedRepresentatives(TRAININGSET *pTS, CODEBOOK *pCB, int numLocalTries)
{
    int i, index, n, numCenters, centerCount, localTrial;
    llong *vector, currentPot, randval, newPot,bestNewPot; /* temprand */

    n = BookSize(pTS);
    numCenters = BookSize(pCB);
    index = irand(0, n - 1);
    CopyVector(Vector(pTS, index), Vector(pCB, 0), VectorSize(pCB));
    vector = GetLongVector(n);

    currentPot = 0;

    for (i=0; i<n-1; i++)
    {
        vector[i] = VectorDist(Vector(pTS, i), Vector(pCB, 0), VectorSize(pCB));
        currentPot += vector[i];
    }

    for (centerCount = 1; centerCount < numCenters; centerCount++)
    {
        bestNewPot = -1;
        // [SS]: Run for numLocalTries and select the clustering with best MSE
        for (localTrial = 0; localTrial < numLocalTries; localTrial++)
        {

            // [SS]: Select a random point as new centroid (=index) with probability based on its distance to already selected centroids
            randval = (llong) (frand() * (float) currentPot);
            /*  temprand = randval; */
            for (index = 0; index < n; index++)
            {
                if (randval <= vector[index]) break;
                else randval -= vector[index];
            }
            newPot = 0;

            // [SS]: Essentially calculate MSE for clustering containing already selected centroids plus new centroid candidate (=index)
            for (i= 0; i < n; i++)
                newPot += min(VectorDist(Vector(pTS, i),
                                         Vector(pTS, index), VectorSize(pCB)), vector[i]);
            // Store the best result
            if ((bestNewPot < 0) || (newPot < bestNewPot))
            {
                bestNewPot = newPot;
                /* bestNewIndex = index; */
            }
        }
        currentPot = bestNewPot;
        for (i= 0; i < n; i++)
            vector[i] = min(VectorDist(Vector(pTS, i),
                                       Vector(pTS, index), VectorSize(pCB)), vector[i]);
        CopyVector(Vector(pTS, index), Vector(pCB, centerCount), VectorSize(pCB));
    }

    free(vector);
}



/* Create codebook (targetCB) by in each iteration selecting the vector from
   sourceCB which has maximum of minimum distances to vectors in targetCB.
 */
void  CreateMaxMinCodebook(TRAININGSET *pTS, CODEBOOK *sourceCB, CODEBOOK *targetCB, int k)
{
    int i,j,i_source,i_target,N;
    int i_dist_max_min;
    //llong dist, dist_min, dist_max_min;
    llong dist_max_min, dist;
    llong* dist_min;
    int* i_dist_min;
    llong* dists;
    N = BookSize(pTS);
    dists = malloc(sizeof(llong)*BookSize(sourceCB));
    i_dist_min = malloc(sizeof(int)*BookSize(sourceCB));

    for(i_source=0; i_source<BookSize(sourceCB); i_source++)
    {
        i_target=0;
        dists[i_source] = MAXLLONG; // Should be 9223372036854775807
        i_dist_min[i_source] = 0;
    }


    //PrintMessage("Ready to loop...\n");
    // BookSize(targetCB) should be >= 1
    for( i=BookSize(targetCB); i<k; i++ )
    {
        dist_max_min=0;

        // Find Maximum of minimum distances
        for(i_source=0; i_source<BookSize(sourceCB); i_source++)
        {
            i_target=i-1;
            dist = VectorDist(Vector(targetCB,i_target), Vector(sourceCB,i_source), VectorSize(pTS));
            if (dist < dists[i_source])
            {
                i_dist_min[i_source] = i_target;
                dists[i_source] = dist;
            }
            dist_min = dists[i_source];

            if(dist_max_min < dist_min)
            {
                dist_max_min = dist_min;
                i_dist_max_min = i_source;
            }
        }

        //printf("i_dist_max_min: %d
        // dist_max_min:%llu\n",i_dist_max_min,dist_max_min);

        ChangeCodebookSize(targetCB, i+1);
        CopyVector( Vector(sourceCB,i_dist_max_min), Vector(targetCB,i), VectorSize(pTS));
    }


    free(dists);
    free(i_dist_min);
}


/* Gonzalez, Teofilo F. "Clustering to minimize the maximum intercluster
   distance." Theoretical Computer Science 38 (1985): 293-306.   "The
   construction of the new cluster is accomplished by first finding a node, v_i,
   in one of the first j clusters (B1, ... , Bj) whose distance to the head of
   the cluster it belongs is maximal." */

void  MaxMinInitialCentroids(TRAININGSET *pTS, CODEBOOK *pCB, int InitialPoint)
{
    int k;
    //int InitialPoint=1;
    VECTORTYPE mean_v;
    k = BookSize(pCB);
    mean_v = CreateEmptyVector(VectorSize(pCB));

    if(!(InitialPoint == 1 || InitialPoint == 2))
    {
        printf("ERROR:InitialPoint\n"); exit(1);
    }

    ChangeCodebookSize(pCB, 1);

    //Choosing first vector for codebook
    if(InitialPoint==1)
    {
        //Random vector from pTS
        CopyVector(Vector(pTS,IRZ(BookSize(pTS))), Vector(pCB,0), VectorSize(pTS));
    }
    else if(InitialPoint==2)
    {
        // Select mean as first centroid
        CodebookCentroid(pTS,mean_v);
        //printf("Mean vector:\n");
        //PrintVector(mean_v,VectorSize(pCB),1);
        //printf("\n");
        CopyVector(mean_v, Vector(pCB,0), VectorSize(pTS));
    }
    else {
        exit(0);
    }



    // TODO: Or select point closest to mean as first centroid

    CreateMaxMinCodebook(pTS, pTS, pCB,  k);

    free(mean_v);
    //exit(0);
}


static int cmpDist(const void* a, const void* b, const void* info)
/* Ascending order assumed */
{

    //llong d1= ((llong*)info)[a];
    int id1=*((int*)a);
    int id2=*((int*)b);
    return (((llong*)info)[id1] < ((llong*)info)[id2] ? 1 : 0);
}

/* J. A. Hartigan and M. A. Wong, “Algorithm AS 136: A k-means clustering algorithm,” Journal of the Royal Statistical Society. Series C (Applied Statistics), vol. 28, no. 1, p. 100–108, 1979.
 * (Appears in additional comments section of article)
 * */

void  SortingHeuristic(TRAININGSET *pTS, CODEBOOK *pCB, int J)
{
    int i,j;
    int randId;
    int N=BookSize(pTS);
    int K=BookSize(pCB);
    float fac = ((float) N)/((float) K);
    VECTORTYPE referencePoint;
    referencePoint = CreateEmptyVector(VectorSize(pCB));

    randId = IRZ(BookSize(pTS));
    CopyVector(Vector(pTS,randId), referencePoint, VectorSize(pTS));
    // Alternative deterministic (original?) way, select mean vector
    // as reference: CodebookCentroid(pTS,referencePoint);


    llong* dists;
    dists = malloc(sizeof(llong)*BookSize(pTS));
    int* idx = malloc(sizeof(int)*BookSize(pTS));

    for(i=0; i<N; i++)
    {
        idx[i] = i;
        dists[i] = VectorDistance(referencePoint, Vector(pTS, i), VectorSize(pTS), MAXLLONG, EUCLIDEANSQ);
    }


    // printf("Reference point:%d\n",randId);
    QuickSort(idx, N, sizeof(int), dists, cmpDist);


    for(i=0; i<K; i++)
    {
        //j = i*(N/K); // Original version, choose at beginning of each bin

        //j = i*(N/K) +(N/K)/2; // Choose at half way of each bin
        j = (int) (i*fac +fac/2); // Choose at half way of each bin
        assert(j < N);
        CopyVector(Vector(pTS, idx[j]), Vector(pCB, i), VectorSize(pTS));
    }

}

/* Bradley and U. Fayyad, Refining initial points for k-means clustering, Int. Conf. on Machine Learning, 91-99, San Francisco, 1998.*/

void  BradleyInitialCentroids(TRAININGSET *pTS, CODEBOOK *pCB, int J)
{
    int k,i;
    PARTITIONING Pnew,P;
    CODEBOOK*  CBsample;
    int N=BookSize(pTS);
    k = BookSize(pCB);
    // Recommended values from article: J=10, sample = 10%
    int numSample=N/10;
    int dsize=BookSize(pTS);
    int* isPointSelected = calloc(BookSize(pTS),sizeof(int)); // Init with zeros
    int centroidBookSize=k*J;
    if (dsize < centroidBookSize) { dsize = centroidBookSize; }

    printf("BradleyInitialCentroids\n  N:%d numSample:%d k:%d J:%d dsize:%d centroidBookSize:%d \n",N,numSample,k,J,dsize, centroidBookSize);
    CODEBOOK**  codebooks = malloc(sizeof(CODEBOOK*)*J);

    llong* distance = malloc(sizeof(llong)*dsize);
    llong* distanceInit = malloc(sizeof(llong)*dsize);
    TRAININGSET* TSsample;
    TSsample = malloc(sizeof(TRAININGSET));


    TRAININGSET* CentroidClusters;
    CentroidClusters = malloc(sizeof(TRAININGSET));


    double totalTime, error, currError, minError;
    int better, iter, totalIter;
    int MaxIter=0;
    int quietLevel=0;
    int useInitial=1;
    int i_k=0;
    minError=-1;

    CreateNewTrainingSet(CentroidClusters,
                         k*J,
                         pTS->BlockSizeX,
                         pTS->BlockSizeY,
                         pTS->BytesPerElement,
                         pTS->MinValue,
                         pTS->MaxValue,
                         pTS->Preprocessing,
                         pTS->GenerationMethod);
    TotalFreq(CentroidClusters)=k*J; //TODO: Should this be set automatically in
                                     // CreateNewTrainingSet

    CreateNewTrainingSet(TSsample,
                         numSample,
                         pTS->BlockSizeX,
                         pTS->BlockSizeY,
                         pTS->BytesPerElement,
                         pTS->MinValue,
                         pTS->MaxValue,
                         pTS->Preprocessing,
                         pTS->GenerationMethod);

    int randId;
    int cid=0;
    int bestCB=0;
    for(cid=0; cid<J; cid++)
    {
        // Select numSample different random points.
        for(i=0; i<numSample; i++)
        {
            randId = IRZ(BookSize(pTS));
            if(isPointSelected[randId] == (cid+1)) { i--; continue; }
            CopyVector(Vector(pTS,randId), Vector(TSsample,i), VectorSize(pTS));
            VectorFreq(TSsample, i) = 1;
            isPointSelected[randId] = cid+1;
        }
        TotalFreq(TSsample)=numSample;
        CBsample = malloc(sizeof(CODEBOOK));
        CreateNewCodebook(CBsample, k, TSsample);
        SelectRandomRepresentatives(TSsample, CBsample);

        CreateNewPartitioning(&Pnew, TSsample, k);
        GenerateOptimalPartitioningGeneral(TSsample, CBsample, &Pnew, MSE);
        CalculateDistances(TSsample, CBsample, &Pnew, distance);

        totalIter = 0;
        currError = error = 0;
        better = iter = 0;

        KMeansIterate(TSsample, CBsample, &Pnew, distance, quietLevel, 0, &iter,
                      totalTime, &error, useInitial, MaxIter);
        codebooks[cid] = CBsample;
        if ( error < minError || cid==0)
        {
            bestCB=cid;
            minError = error;
        }

        for(i_k=0; i_k<k; i_k++)
        {
            CopyVector(Vector(CBsample,i_k), Vector(CentroidClusters,cid*k+i_k), VectorSize(pTS));
            VectorFreq(CentroidClusters, cid*k+i_k) = 1;
        }

    }

    cid=0;

    CODEBOOK* best = malloc(sizeof(CODEBOOK));
    CODEBOOK* cbc = malloc(sizeof(CODEBOOK));
    CreateNewCodebook(best, k, CentroidClusters);
    CreateNewCodebook(cbc, k, CentroidClusters);
    for(cid=0; cid<J; cid++)
    {
        CBsample=codebooks[cid];
        for (i=0; i< k; i++)
        {
            CopyVector(Vector(CBsample,i), Vector(cbc,i), VectorSize(CBsample));
        }
        CreateNewPartitioning(&P, CentroidClusters, k);
        GenerateOptimalPartitioningGeneral(CentroidClusters, cbc, &P, MSE);
        CalculateDistances(CentroidClusters, cbc, &P, distance);

        KMeansIterate(CentroidClusters, cbc, &P, distance, quietLevel=0, 0, &iter,
                      totalTime, &error, useInitial, MaxIter);
        if ( error < minError || cid==0)
        {
            CopyCodebook(cbc,best);
            //FreeCodebook(best);
            bestCB=cid;
            minError = error;
            //best = cbc;
        }
        else {
            //FreeCodebook(cbc);
        }
        FreePartitioning(&P);
    }

    for (i=0; i< k; i++)
    {
        CopyVector(Vector(best,i), Vector(pCB,i), VectorSize(CBsample));
    }

    for(cid=0; cid<J; cid++)
    {
        //TODO: Should free codebooks?
        //FreeCodebook(codebooks[cid]);
    }

    printf("Best CB:%d %f\n",bestCB,minError);

    //TODO: This part of method not yet implemented:
    // Re assignment of empty clusters: If any of the K clusters have no
    // membership (which often happens when clustering over small subsamples),
    // the corresponding initial estimates of the empty cluster centroids are
    // set to data elements which are farthest from their assigned cluster
    // center, and classic K-Means is called again from these new initial
    // centriods

    free(isPointSelected);
    free(codebooks);
    free(distance);
    free(distanceInit);
    free(TSsample);
    free(CentroidClusters);
}


void  RandomPartitionInitialCentroids(TRAININGSET *pTS, CODEBOOK *pCB, int NumTries)
{
    int k,i,i_try;
    double error;
    double min_error = -1;
    CODEBOOK* CB2;

    PARTITIONING P;
    VECTORTYPE mean_v;
    k = BookSize(pCB);
    CB2=malloc(sizeof(CODEBOOK));

    CreateNewCodebook(CB2, k, pTS);

    for(i_try=0; i_try < NumTries; i_try++)
    {
        CreateNewPartitioning(&P, pTS, k);

        //PutAllInFirstPartition(pTS, &P);
        for(i=0; i<BookSize(pTS); i++)
        {
            AddToPartition(pTS, &P, IRZ(k), i);
        }

        GenerateOptimalCodebook(pTS, CB2, &P);
        error = AverageErrorForSolution(pTS, CB2, &P, MSE);
        //printf("MSE:%f",error);
        if ( error < min_error || i_try==0)
        {
            min_error = error;
            CopyCodebook(CB2,pCB);
            //printf(" (improved) ");
        }
        //printf("\n");

        FreePartitioning(&P);
    }

    printf("Performed %d random partitions\n",i_try);

    FreeCodebook(CB2);
    free(CB2);
}



// Vector minus a-b
int* vecminus(int* a, int* b, int size) {
    int* r = malloc(sizeof(int)*size);
    int i;
    for (i=0; i<size; i++)
    {
        r[i] = a[i]-b[i];
    }
    return r;
}

// Dot product a*b
float dotprod(int *a,int *b, int size)
{
    int i;
    float A,B;
    float dotp = 0.0f;
    for (i=0; i<size; i++)
    {
        A = (float)(a[i]);
        B = (float)(b[i]);
        dotp += A*B;
    }
    return dotp;
}

static int cmpFloat(const void* a, const void* b, const void* info)
{
    int id1=*((int*)a);
    int id2=*((int*)b);
    return (((float*)info)[id1] < ((float*)info)[id2] ? 1 : 0);
}

// Sort ints based on float val in descending order
static int cmpFloatDesc(const void* a, const void* b, const void* info)
{
    int id1=*((int*)a);
    int id2=*((int*)b);
    return (((float*)info)[id1] > ((float*)info)[id2] ? 1 : 0);
}

static int cmpDoubleDesc(const void* a, const void* b, const void* info)
{
    int id1=*((int*)a);
    int id2=*((int*)b);
    return (((double*)info)[id1] > ((double*)info)[id2] ? 1 : 0);
}


int SplitCluster(TRAININGSET* TS, CODEBOOK* CBleft, CODEBOOK* CBright, TRAININGSET*  TSleft, TRAININGSET*  TSright, int kmIters) {

    PARTITIONING*  pP;
    PARTITIONING*  Pleft;
    PARTITIONING*  Pright;
    CODEBOOK* CBnew;

    CODEBOOK* CBparent = allocate(sizeof(CODEBOOK));

    pP=allocate(sizeof(PARTITIONING));
    Pleft=allocate(sizeof(PARTITIONING));
    Pright=allocate(sizeof(PARTITIONING));

    CBnew=allocate(sizeof(CODEBOOK));

    int useInitial=0;
    int InitMethod=1; //1=random;
    int max_iter=10;
    int iter=0;
    int numc=2;

    double totalTime, error;

    llong distanceInit[BookSize(TS)];
    llong distance[BookSize(TS)];

    CreateNewCodebook(CBparent, 2, TS);
    CreateNewPartitioning(pP, TS, 2);
    SelectRandomRepresentatives(TS, CBparent);
    GenerateOptimalPartitioningGeneral(TS, CBparent, pP, MSE);


    if(kmIters > 0)
    {
        // Default: Ten iterations of kmeans.
        KMeansIterate(TS, CBparent, pP, distance, 0, 0, &iter,
                      totalTime, &error, useInitial,kmIters);
    }
    // Else: just divide by random pair

    GetSubspaceForCentroid(TSleft,pP,Pleft,CBparent,CBleft,TS,0 /*label*/,1);
    GetSubspaceForCentroid(TSright,pP,Pright,CBparent,CBright,TS,1 /*label*/,1);
    TotalFreq(TSleft) = BookSize(TSleft); //TODO: should be done automatically
                                          // in GetSubspaceForCentroid ?
    TotalFreq(TSright) = BookSize(TSright);

    FreeCodebookFull(CBparent);
    FreePartitioningFull(pP);
    FreePartitioningFull(Pleft);
    FreePartitioningFull(Pright);

    return 0;
}

void hierarchicalSplittingCodebook(TRAININGSET* TS, CODEBOOK* pCBnew, int kmIters) {
    VECTORTYPE u;
    CODEBOOK CB;
    CODEBOOK*      pCB =&CB;
    CODEBOOK*      CBleft;
    CODEBOOK*      CBright;
    PARTITIONING P;
    PARTITIONING*  pP =&P;
    TRAININGSET*  curTS;

    TRAININGSET* TSleft;
    TRAININGSET* TSright;

    llong distanceInit[BookSize(TS)];
    int targetClusters = BookSize(pCBnew);

    TRAININGSET*    pTS =TS;
    PARTITIONING Pnew, Pinit;
    CODEBOOK CBnew, CBinit;

    TRAININGSET** ts_list;
    VECTORTYPE** centroids;
    int N = BookSize(TS);
    if (targetClusters > N) {targetClusters=N; }

    int i=0;
    int numc=1;
    int useInitial=0;
    int nn_id;
    double dleft;
    double dright;
    int num_ts=0;

    VECTORTYPE leftC;
    VECTORTYPE rightC;

    ts_list = allocate(sizeof(TRAININGSET*)*BookSize(TS));
    centroids = allocate(sizeof(VECTORTYPE*)*BookSize(TS));

    CreateNewCodebook(pCB, numc, pTS);
    CreateNewPartitioning(pP, pTS, numc);

    InitializeSolutions(TS, pCB, pP, &CBnew, &Pnew, &CBinit, &Pinit,
                        distanceInit, numc, useInitial);
    FreePartitioning(&Pnew);
    FreePartitioning(&Pinit);

    CBleft = allocate(sizeof(CODEBOOK));
    CBright = allocate(sizeof(CODEBOOK));
    TSleft = allocate(sizeof(TRAININGSET));
    TSright = allocate(sizeof(TRAININGSET));

    ts_list[0] =malloc(sizeof(CODEBOOK));
    CreateNewCodebook(ts_list[0], BookSize(TS), TS);
    CopyCodebook(TS, ts_list[0]);
    num_ts=1;
    centroids[0] = NULL;


    int part_id=0;
    for(i=0; i<targetClusters-1; i++)
    {
        int new_part_id=i+1;
        int max_size=0;

        if(0)
        {
            // Experimental. Choose random cluster instead of largest.
            part_id = IRZ(i+1);
            while ( BookSize(ts_list[part_id]) < 2) {
                part_id = IRZ(part_id+1);
            }
        }
        else {
            //Find Cluster with largest size
            //TODO: this step could be speeded up by keeping a sorted list.
            for(int i_ts=0; i_ts<num_ts; i_ts++)
            {
                if(BookSize(ts_list[i_ts]) > max_size)
                {
                    part_id=i_ts; max_size=BookSize(ts_list[i_ts]);
                }
            }
        }
        curTS=ts_list[part_id];
        /*printf("Split iter:%d part_id:%d
           size:%d\n",i,part_id,BookSize(curTS));*/
        /*dump_TS(curTS);*/
        TSleft = allocate(sizeof(TRAININGSET));
        TSright = allocate(sizeof(TRAININGSET));
        CBleft = allocate(sizeof(CODEBOOK));
        CBright = allocate(sizeof(CODEBOOK));
        SplitCluster(curTS,CBleft,CBright,TSleft,TSright,kmIters);
        //printf("L:%d, R:%d\n",TSleft->CodebookSize, TSright->CodebookSize);

        leftC=allocate(VectorSize(TS)*sizeof(VECTORELEMENT));
        CopyVector(Vector(CBleft,0), leftC, VectorSize(CBleft));

        rightC=allocate(VectorSize(TS)*sizeof(VECTORELEMENT));
        CopyVector(Vector(CBright,0), rightC, VectorSize(CBright));

        FreeCodebookFull(ts_list[part_id]);
        ts_list[part_id] = TSleft;
        ts_list[new_part_id] = TSright;
        if(centroids[part_id]!=NULL) {deallocate(centroids[part_id]); }
        centroids[part_id] = leftC; //TODO: use graph vector->data instead of
                                    // centroids[];
        centroids[new_part_id] = rightC;

        FreeCodebook(CBleft);
        FreeCodebook(CBright);
        deallocate(CBleft);
        deallocate(CBright);

        num_ts++;
    }

    for(i=0; i<targetClusters; i++)
    {
        CopyVector(centroids[i], Vector(pCBnew,i), VectorSize(TS));
        deallocate(centroids[i]);
        FreeCodebookFull(ts_list[i]);
    }
    deallocate(centroids);
    deallocate(ts_list);

}


void  DensitySampledCentroids(TRAININGSET *pTS, CODEBOOK *pCB, int variant) {
    //int KMeansDensityPeaks(TRAININGSET *pTS, CODEBOOK *pCB, PARTITIONING *pP,
        //int clus, int repeats, int quietLevel) {
    printf("DensitySampledCentroids\n");
    int i,K;
    int N = BookSize(pTS);
    //int numSamples = 5000;
    //int numSamples = sqrt(log(N))*sqrt(N);
    K=BookSize(pCB);
    int numSamples = variant*K;
    //int numSamples = sqrt(N);
    //int numSamples = N/10;
    //int numSamples = N;
    if(numSamples > N ) {numSamples=N;}
    printf("numSamples=%d\n",numSamples);
    int* samples = getRandomSampleInts(N-1, numSamples);
    TRAININGSET *sampleTS;
    CODEBOOK *sampleCB=malloc(sizeof(CODEBOOK));

    //sqrt(log(N))*sqrt(N)

    printf("KMPEAKS\n");
    sampleTS = createSampleDataset(pTS, numSamples);
    //PrintCodebook(sampleTS);

    //CreateNewCodebook(pCB, K, pTS);
    //CreateNewCodebook(sampleCB, clus, sampleTS);
    //CreateNewCodebook(sampleCB, clus, pTS);

    DensityInitialCentroids(sampleTS, pCB, 0);
    return;
}

static int cmpDoubleAsc(const void* a, const void* b, const void* info)
{
    int id1=*((int*)a);
    int id2=*((int*)b);
    return (((double*)info)[id1] < ((double*)info)[id2] ? 1 : 0);
}

double* calcDeltasFromWindow( TRAININGSET *pTS, double* density, int* densityOrder) {

    int i,j,K,N,D,higherDenPoint,curPoint;
    N=BookSize(pTS);
    D=VectorSize(pTS);
    double dist,minDist;
    double* dists = calloc(N,sizeof(double));


    double* dimvals = calloc(N,sizeof(double)); // Init to zero
    int* evidence = calloc(N,sizeof(int)); // Init to zero
    int* idx = malloc(sizeof(int)*N);

    int dim = 0;

    for(dim=0;dim<D;dim++) {
        for(i=0; i<N; i++)
        {
            idx[i]=i;
            dimvals[i] = (double) VectorScalar(pTS,i,dim);
        }

        QuickSort(idx, N, sizeof(int), dimvals, cmpDoubleAsc);

        K=25;
        for(i=0; i<N; i++)
        {
            int startw = i - K;
            if(startw < 0) { startw = 0;}
            int endw = startw + K*2;
            if(endw >= N) { endw = N-1; startw=endw-2*K;}

            for(j=startw; j<endw; j++)
            {
                if (density[idx[j]] > density[idx[i]]) {
                    evidence[idx[i]]++;
                }
            }
            if(dim == D -1) {
                // printf("i=%d evidence=%d\n",i,evidence[idx[i]]);
            }
        }
    }


    // printf("Candidates:\n");
    int candidateCount=0;
    for(i=0; i<N; i++)
    {
        curPoint=densityOrder[i];
        dists[curPoint] = 0.0001;
        //if(evidence[curPoint] < K) {
        if(1) {
            // printf("i=%d evidence=%d\n",i,evidence[i]);
            candidateCount++;

            /********/
            for(j=i-1; j>=0; j--)
            {
                higherDenPoint=densityOrder[j];
                dist=(double) VectorDistance(Vector(pTS,curPoint), Vector(pTS,higherDenPoint), D, MAXLLONG, EUCLIDEAN);
                if(j==i-1) { minDist=dist; }
                if(dist < minDist) { minDist=dist; }
            }
            dists[curPoint] = minDist;
            /********/
        }
    }
    // printf("candidateCount=%d\n",candidateCount);



    return dists;

}



double* calcDeltasPerDim( TRAININGSET *pTS, double* density, int* densityOrder ) {

    int i,j,K,N,D,higherDenPoint,curPoint;
    N=BookSize(pTS);
    D=VectorSize(pTS);
    double dist,minDist;
    double* dists = calloc(N,sizeof(double));
    //double** distsPerDim = malloc(sizeof(double*)*D);
    //int* idx = malloc(sizeof(int)*N);
    //int* deltaOrder = malloc(sizeof(int)*N);


    // printf("DP caclDeltasPerDim\n");

    for(i=0; i<D; i++)
    {
        //distsPerDim[i] = malloc(sizeof(double)*N);
    }

    dists[densityOrder[0]] = FLT_MAX;


   //QuickSort(idx2, N, sizeof(int), dists, cmpDoubleDesc);

    // "d_i is measured by computing the minimum distance between the point i
    // and any other point with higher density:"
    //Calculate distances to point with higher density

    //Should be max(d_ij) according to article, but does not affect the sorting

    for(int dim=0; dim<D; dim++)
    {
        for(i=1; i<N; i++)
        {
            curPoint=densityOrder[i];
            for(j=i-1; j>=0 && kmOpt.densityPeaksFilterMethod!=5; j--)
            {
                higherDenPoint=densityOrder[j];
                //dist=(double) VectorDistance(Vector(pTS,curPoint), Vector(pTS,higherDenPoint), D, MAXLLONG, EUCLIDEAN);
                //dist=(double) VectorDistance(Vector(pTS,curPoint), Vector(pTS,higherDenPoint), D, MAXLLONG, EUCLIDEAN);
                dist = abs(VectorScalar(pTS,curPoint,dim) - VectorScalar(pTS,higherDenPoint,dim));
                // If first time
                if(j==i-1) { minDist=dist; }

                if(dist < minDist) { minDist=dist; }
            }
            //distsPerDim[D][curPoint] = minDist;
            dists[curPoint] += minDist;
        }
    }

    for(i=0; i<N; i++)
    {
        dists[i] = density[i]*dists[i];
    }


    return dists;

}



double* calcDeltas( TRAININGSET *pTS, double* density, int* densityOrder ) {

    int i,j,K,N,D,higherDenPoint,curPoint;
    N=BookSize(pTS);
    D=VectorSize(pTS);
    double dist,minDist;
    double* dists = malloc(sizeof(double)*N);


    printf("DP filter method:%d\n",kmOpt.densityPeaksFilterMethod);


    // "d_i is measured by computing the minimum distance between the point i
    // and any other point with higher density:"
    //Calculate distances to point with higher density

    //Should be max(d_ij) according to article, but does not affect the sorting
    dists[densityOrder[0]] = FLT_MAX;

    for(i=1; i<N; i++)
    {
        curPoint=densityOrder[i];
        for(j=i-1; j>=0 && kmOpt.densityPeaksFilterMethod!=5; j--)
        {
            higherDenPoint=densityOrder[j];
            dist=(double) VectorDistance(Vector(pTS,curPoint), Vector(pTS,higherDenPoint), D, MAXLLONG, EUCLIDEAN);
            if(j==i-1) { minDist=dist; }
            if(dist < minDist) { minDist=dist; }
        }
        if (kmOpt.densityPeaksFilterMethod==0)
        {
            dists[curPoint] = minDist;
        }
        else if ( kmOpt.densityPeaksFilterMethod==1 )
        {
            dists[curPoint] = minDist*density[curPoint];
        }
        else if ( kmOpt.densityPeaksFilterMethod==5 )
        {
            // TODO: remove?
            dists[curPoint] = density[curPoint];
        }
        else if ( kmOpt.densityPeaksFilterMethod==99 )
        {
            dists[curPoint] = minDist;
        }


    }

    return dists;

}

/* Density peaks algorithm:
 * A. Rodriguez and A. Laio. Clustering by fast search and find of density
 **peaks. Science, 344(6191):1492-1496 2014.

 * For density estimation, use a method from Mitra 2002, but disregard factors
 **affecting the scale since only order of the densities are of interest here.
 *
 * Mitra, Pabitra, C. A. Murthy, and Sankar K. Pal. "Density-based multiscale
 **data condensation." IEEE Transactions on pattern analysis and machine
 **intelligence 24.6 (2002): 734-747.
 */
void  DensityInitialCentroids(TRAININGSET *pTS, CODEBOOK *pCB, int variant) {
    int i,j,K,N,D,higherDenPoint,curPoint;
    double dist,minDist;
    N=BookSize(pTS);
    D=VectorSize(pTS);
    int* idx = malloc(sizeof(int)*N);
    int* deltaOrder = malloc(sizeof(int)*N);
    K=BookSize(pCB);
    double* density = malloc(sizeof(double)*N);
    //double* dists = malloc(sizeof(double)*N);
    double* deltas = malloc(sizeof(double)*N);
    double densityTime;
    //int knnK=8;
    //int knnK=20;
    //int knnK=10;
    int knnK=(N/K)/2;
    if(knnK > 50) {knnK=50;}
    if(knnK < 2) {knnK=2;}
    if(variant > 2) {knnK=variant;}
    printf("N=%d knnK:%d filterMethod=%d\n",N,knnK,kmOpt.densityPeaksFilterMethod);

    //Graph* g;
    //printf("Initial centroids based on density\n");
    //printf("knnK=%d\n",knnK);
    //g = bruteForcekNNGraph(pTS, knnK,0);
    for(i=0; i<N; i++)
    {
        idx[i]=i;
        deltaOrder[i]=i;
        //density[i] = 1/((float)(g->vectors[i]->distances[knnK-1]));
    }

    SetClock(&densityTime);
    if ( kmOpt.densityMethod==0) {
        printf("DensityMethod DDDE\n");
        density = dimBasedDensity( pTS, knnK,0);
    }
    else if ( kmOpt.densityMethod==1) {
        printf("DensityMethod kNN graph\n");
        density = knnGraphDensity(pTS, knnK, kmOpt.sample);
    }
    else if ( kmOpt.densityMethod==2) {
        printf("DensityMethod DDDE (sliding window)\n");
        density = dimBasedDensity( pTS, knnK,1);
    }

    else if ( kmOpt.densityMethod==3) {
        printf("Random density\n");
        for (i=0;i<N;i++) {
            density[i] = drand();
        }
    }


    PrintMessage("density_time=%f\n", GetClock(densityTime));


    QuickSort(idx, N, sizeof(int), density, cmpDoubleDesc);

    QuickSort(deltaOrder, N, sizeof(int), deltas, cmpDoubleDesc);


    if(kmOpt.densityPeaksDeltaMethod==0) {
        deltas = calcDeltas( pTS, density, idx );
    } else if (kmOpt.densityPeaksDeltaMethod==1) {
        deltas = calcDeltasFromWindow( pTS, density, idx );
    }
    else {
        printf("ERROR:Unspecified delta method\n"); exit(1);
    }

    //Debug, print deltas and dists
    if ( kmOpt.densityPeaksFilterMethod==99 || kmOpt.QuietLevel >= 4 ) {
        printf("=== START delta,rho ===\n");
        for(i=0; i<N; i++)
        {
            printf("%f %0.10f\n",deltas[i],density[i]);
        }
        printf("=== END delta,rho ===\n");
    }


    // Sort by minimum distance to higher density point
    QuickSort(deltaOrder, N, sizeof(int), deltas, cmpDoubleDesc);


    // Take first K points
    for(i=0; i<K; i++)
    {
        CopyVector(Vector(pTS, deltaOrder[i]), Vector(pCB, i), D);
    }

    free(idx);
    free(deltaOrder);
    free(density);
    free(deltas);

    return;
}

void  ProjectionInitialCentroids(TRAININGSET *pTS, CODEBOOK *pCB, int variant) {
    int pointA,pointB,i,j,K,N,D,maxInd;
    llong dist,maxDist;
    N=BookSize(pTS);
    D=VectorSize(pTS);
    K=BookSize(pCB);
    float* scores=malloc(sizeof(float)*N);
    int* idx = malloc(sizeof(int)*N);
    VECTORTYPE projAxis;
    printf("Initial centroids with projection ");
    if(variant == 1) { printf("(two random points)\n"); }
    else if(variant == 2)
    {
        printf("(furthest of a random point)\n");
    }

    for(i=0; i<N; i++)
    {
        idx[i] = i;
    }
    pointA = irand(0,N-1);
    if(variant==1)
    {
        pointB = pointA;
        while(pointB == pointA) {
            pointB = irand(0,N-1);
        }
    }
    else if (variant==2)
    {

        maxDist=0;
        for(i=0; i<N; i++)
        {
            if(i==pointA) { continue; }
            dist = VectorDist(Vector(pTS,pointA), Vector(pTS,i), D);
            if (dist > maxDist)
            {
                maxDist = dist;
                maxInd = i;
            }
        }
        pointB=maxInd;

    }
    else {assert(0); }
    if(kmOpt.QuietLevel >=3) {

        printf("==== START Projection Points ====\n");
        //printf("A:%d B:%d\n",pointA,pointB);
        printf("%d %d\n",pointA,pointB);
        PrintVector(Vector(pTS,pointA),VectorSize(pTS),1); printf("\n");
        PrintVector(Vector(pTS,pointB),VectorSize(pTS),1); printf("\n");
        printf("==== END Projection Points ====\n");

    }

    projAxis=vecminus(Vector(pTS,pointA),Vector(pTS,pointB),D);

    for( i=0; i<N; i++)
    {
        scores[i] = dotprod(projAxis,Vector(pTS,i),D);
        //printf("s:%f\n",scores[i]);
    }

    QuickSort(idx, N, sizeof(int), scores, cmpFloat);

    for(i=0; i<K; i++)
    {
        j = i*(N/K); // Choose at start of each bin
        //j = i*(N/K) +(N/K)/2; // Choose at half way of each bin
        assert(j < N);
        CopyVector(Vector(pTS, idx[j]), Vector(pCB, i), VectorSize(pTS));
    }


    free(projAxis);
    free(scores);
    free(idx);

    return;
}

/*-------------------------------------------------------------------*/
/* Random selection by Luxburg:                                      */
/* von Luxburg, Clustering stability: an overview                    */
/* Foundations and Trends in Machine Learning, 2010                  */
/*-------------------------------------------------------------------*/


void  LuxburgInitialCentroids(TRAININGSET *pTS, CODEBOOK *pCB)
{
    int i,j, new, k, L, size;
    CODEBOOK CB2;
    CODEBOOK CBreduced;
    TRAININGSET TStmp;
    PARTITIONING P;
    int i_red,i_dist_min, i_dist_max_min;
    llong dist, dist_min, dist_max_min;
    int remove_count = 0;

    k    = BookSize(pCB);
    L    = (int) (k * log(k)/log(2));
    size = BookSize(pTS)/L;
    if ( L > BookSize(pTS)) { L = BookSize(pTS); }


    // Step 1: L = k*log(k) initial centroids.

    CreateNewCodebook(&CB2, L, pTS);
    CreateNewCodebook(&CBreduced, L, pTS);

    //Random codebook shuffles training set, so make a temporary for that.

    CreateNewCodebook(&TStmp, BookSize(pTS), pTS);
    CopyCodebook(pTS, &TStmp);
    RandomCodebook(&TStmp, &CB2);
    FreeCodebook(&TStmp);
    CreateNewPartitioning(&P, pTS, L);

    // Step 2: One iteration of k-means
    GenerateOptimalPartitioning(pTS, &CB2, &P);
    GenerateOptimalCodebook(pTS, &CB2, &P);

    // Step 3: Remove small centroids (size<N/L)
    int newcbsize=0;
    for( i=0; i<BookSize(&CB2); i++)
    {
        // Keep at least k centroids.
        if( CCFreq(&P,i) < size && BookSize(&CB2) - remove_count > k)
        {
            remove_count++;
        }
        else {

            CopyVector( Vector(&CB2,i), Vector(&CBreduced,newcbsize), VectorSize(pTS));
            //AddToCodebook(&CBreduced, Vector(&CB2,i)); // Produces segfault
            // (??)
            newcbsize++;
        }
    }
    if(remove_count > 0)
    {
        DecreaseCodebookSize(&CBreduced, newcbsize);
    }

    // Step 4: Heuristic selection (MaxMin) aiming at even distribution

    // First vector randomly
    ChangeCodebookSize(pCB, 1);
    new = IRZ(BookSize(&CBreduced));
    CopyVector( Vector(&CBreduced,new), Vector(pCB,0), VectorSize(pTS));
    FreePartitioning(&P);
    CreateNewPartitioning(&P, pTS, BookSize(&CBreduced));

    //TODO: is this needed?:
    GenerateOptimalPartitioning(&CBreduced, pCB, &P);

    // From article: "Repeat until K centers are selected: Select the next
    // center as the one that maximizes the minimum distance to the centers
    // already selected."
    CreateMaxMinCodebook(pTS, &CBreduced, pCB,  k);

    return;
}

void decreaseOverlapInit(TRAININGSET *TS, CODEBOOK *pCB, int variant) {
    int i,j,pointA,pointB;
    float p;
    printf("decreaseOverlapInit\n");

    TRAININGSET * oiTS; // Overlap increased TS
    TRAININGSET * oiCB; // Overlap increased CB
    PARTITIONING P;
    int N = BookSize(TS);
    int Nnew = N*variant;
    int k = BookSize(pCB);
    oiTS = malloc(sizeof(TRAININGSET));
    oiCB = malloc(sizeof(CODEBOOK));
    //VECTORTYPE

    //CreateNewCodebook(oiTS, Nnew, TS);
    CreateNewTrainingSet(oiTS,
                         Nnew,
                         TS->BlockSizeX,
                         TS->BlockSizeY,
                         TS->BytesPerElement,
                         TS->MinValue,
                         TS->MaxValue,
                         TS->Preprocessing,
                         TS->GenerationMethod);
    TotalFreq(oiTS)=Nnew;


    for (i=0; i<N; i++)
    {
        CopyNode(&Node(TS,i), &Node(oiTS,i), VectorSize(TS));
        //CopyVector( Vector(TS, i), Vector(oiTS, i), VectorSize(TS));
    }

    for (i=N; i<Nnew; i++)
    {

        p=frand();

        // Take mean of random points A and B
        pointA = irand(0, N - 1);
        pointB = irand(0, N - 1);
        //CopyNode(&Node(TS,i), &Node(oiTS,i), VectorSize(TS));
        for( j=0; j<VectorSize(TS); j++ )
        {
            VectorScalar(oiTS,i,j) = (llong) (p*VectorScalar(TS,pointA,j) + (1-p)*VectorScalar(TS,pointB,j));
            //VectorScalar(oiTS,i,j) = (VectorScalar(TS,pointA,j) + VectorScalar(TS,pointB,j))/2;
            VectorFreq(oiTS, i) = 1;
        }
    }

    CreateNewCodebook(oiCB, k, oiTS);

    CreateNewPartitioning(&P, oiTS, k);

    PerformKMeans(oiTS, oiCB, &P, k, 1 /*repeats*/,  1 /*InitMethod*/,  0 /*InitMethodVar*/,
            0 /*quietLevel*/, 0 /*useInitial*/, 0 /*MaxIter*/);

    for (i=0; i<k; i++)
    {
        //CopyNode(&Node(oiCB,i), &Node(pCB,i), VectorSize(TS));
        CopyVector(Vector(oiCB,i), Vector(pCB,i), VectorSize(TS));
    }


    FreeCodebookFull(oiTS);
    FreeCodebookFull(oiCB);
    FreePartitioning(&P);
}



