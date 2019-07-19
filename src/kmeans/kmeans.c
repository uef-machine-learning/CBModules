/*------------------------------------------------------------------------*/
/* KMEANS.C     Marko Tuononen                                            */
/*              Pasi Fränti                                               */
/*                                                                        */
/* Fast repeated K-means implementation. Uses activity detection from     */
/*                                                                        */
/*    "A fast exact GLA based on code vector activity detection"          */
/*    IEEE Trans. on Image Processing, 9 (8), 1337-1342, August 2000.     */
/*    Timo Kaukoranta, Pasi Fränti and Olli Nevalainen                    */
/*                                                                        */
/* Naming conventions used in the code                                    */
/*                                                                        */
/*    TS        training set (data objects)                               */
/*    CB        codebook (cluster representatives, centroids)             */
/*    P         partitioning (pointing from TS to CB)                     */
/*    p-prefix  pointer, e.g. pTS is pointer to the training set TS       */
/*                                                                        */
/* Changelog:                                                             */
/*                                                                        */
/* 0.65 (SS): Luxburg, MaxMin and Random partition init. (4.4.17)         */
/* 0.64 (PF): Luxburg initialization method (20.9.16)                     */
/* 0.63 (MM): Changed -R parameter max to 1.000.000                       */
/* 0.62 (MM): Something                                                   */
/*                                                                        */
/*------------------------------------------------------------------------*/

#define ProgName       "KMEANS"
#define VersionNumber  "Version 0.65" /* SS */
#define LastUpdated    "20.9.2016"

/*-------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "cb.h"
#include "random.h"
#include "interfc.h"
#include "reporting.h"
//#include "rs.h"
#include "memctrl.h"

#include "kmeans.h"
#include "inits.h"
#include "graph.h"


// Global variable
KmeansOpt kmOpt;


/* =================== PUBLIC FUNCTIONS ============================= */


// Hybrid method that calculates initial solution to sampled data set
// with density peaks and improves with K-means.

TRAININGSET*  createSampleDataset(TRAININGSET *pTS ,int numSamples) {
    TRAININGSET* TSsample;
    TSsample = malloc(sizeof(TRAININGSET));
    int N = BookSize(pTS);
    int i,randId;

    CreateNewTrainingSet(TSsample,
            numSamples,
            pTS->BlockSizeX,
            pTS->BlockSizeY,
            pTS->BytesPerElement,
            pTS->MinValue,
            pTS->MaxValue,
            pTS->Preprocessing,
            pTS->GenerationMethod);

    int* samples = getRandomSampleInts(N-1, numSamples);

    // Select numSample different random points.
    for(i=0; i<numSamples; i++)
    {
        randId=samples[i];
        CopyVector(Vector(pTS,randId), Vector(TSsample,i), VectorSize(pTS));
        VectorFreq(TSsample, i) = 1;
    }

    return TSsample;
}


int KMeansDensityPeaks(TRAININGSET *pTS, CODEBOOK *pCB, PARTITIONING *pP,
        int clus, int repeats, int quietLevel) {

    int i;
    int N = BookSize(pTS);
    //int numSamples = 500;
    int numSamples = sqrt(log(N))*sqrt(N);
    //int numSamples = sqrt(N);
    printf("numSamples=%d\n",numSamples);
    int* samples = getRandomSampleInts(N-1, numSamples);
    TRAININGSET *sampleTS;
    CODEBOOK *sampleCB=malloc(sizeof(CODEBOOK));

    //sqrt(log(N))*sqrt(N)

    printf("KMPEAKS\n");
    sampleTS = createSampleDataset(pTS, numSamples);
    PrintCodebook(sampleTS);

    CreateNewCodebook(pCB, clus, pTS);
    //CreateNewCodebook(sampleCB, clus, sampleTS);
    CreateNewCodebook(sampleCB, clus, pTS);

    DensityInitialCentroids(sampleTS, pCB, 20);
    PrintCodebook(pCB);

    //KMeansIterate(pTS, pCB, &Pnew, distance, quietLevel, i, &iter,
    //totalTime, &error, 1 /*Not used?*/, 0 /*MaxIter*/);


    return 0;
}



/* Gets training set pTS (and optionally initial codebook pCB or
   partitioning pP) as a parameter, generates solution (codebook pCB +
   partitioning pP) and returns 0 if clustering completed successfully.
   N.B. Random number generator (in random.c) must be initialized! */

int PerformKMeans(TRAININGSET *pTS, CODEBOOK *pCB, PARTITIONING *pP,
		  int clus, int repeats, int InitMethod, int InitMethodVar,
		  int quietLevel, int useInitial, int MaxIter)
{
    // If MaxIter > 0, limit k-means iterations to MaxIter.
  PARTITIONING  Pnew, Pinit;
  CODEBOOK      CBnew, CBinit;
  llong         distance[BookSize(pTS)];
  llong         distanceInit[BookSize(pTS)];
  double        totalTime, error, currError, initTime;
  int           i, better, iter, totalIter;

  SetClock(&totalTime);
  totalIter = 0;
  currError = error = 0;

  if ((clus < 1) || (BookSize(pTS) < clus) || (repeats < 0))
    {
    return 1;   /* clustering failed */
    }

  //printf("start_time=%f\n", GetClock(totalTime));
  InitializeSolutions(pTS, pCB, pP, &CBnew, &Pnew, &CBinit, &Pinit,
                      distanceInit, clus, useInitial);


  PrintHeader(quietLevel);
  //printf("InitMethodVar:%d\n",InitMethodVar);

  /* perform repeats time full K-means */
  for (i = 0; i < repeats || (i==0 && repeats==0); i++)
    {
    better = iter = 0;


    SetClock(&initTime);
    GenerateSolution(pTS, &CBnew, &Pnew, &CBinit, &Pinit, distance,
		     distanceInit, InitMethod, InitMethodVar, useInitial);
    printf("init_time=%f\n", GetClock(initTime));
    if(repeats == 0)
    {
        printf("Initialize codebook only. No k-means.\n");
        error = AverageErrorForSolution(pTS, &CBnew, &Pnew, MSE);
    }
    else
    {
        KMeansIterate(pTS, &CBnew, &Pnew, distance, quietLevel, i, &iter,
                totalTime, &error, useInitial, MaxIter);
    }

    totalIter += iter;

    /* got better result */
    if ((i == 0) || (error < currError))
      {
      CopyCodebook(&CBnew, pCB);
      CopyPartitioning(&Pnew, pP);
      currError = error;
      better = 1;
      }

    PrintRepeat(quietLevel, repeats, i, iter, error, GetClock(totalTime), better);
    }

  PrintFooterKM(quietLevel, currError, repeats, GetClock(totalTime), totalIter);

  FreeCodebook(&CBnew);
  FreePartitioning(&Pnew);
  FreeCodebook(&CBinit);
  FreePartitioning(&Pinit);

  return 0;
}


/* ------------------------------------------------------------------ */


char* KMeansInfo(void)
{
  char* p;
  int len;

  len = strlen(ProgName)+strlen(VersionNumber)+strlen(LastUpdated)+4;
  p   = (char*) malloc(len*sizeof(char));

  if (!p)
  {
    ErrorMessage("ERROR: Allocating memory failed!\n");
    ExitProcessing(FATAL_ERROR);
  }

  sprintf(p, "%s\t%s\t%s", ProgName, VersionNumber, LastUpdated);

  return p;
}



/* ====================== PRIVATE FUNCTIONS ========================= */


void InitializeSolutions(TRAININGSET *pTS, CODEBOOK *pCB,
       PARTITIONING *pP, CODEBOOK *pCBnew, PARTITIONING *pPnew, CODEBOOK *pCBinit,
       PARTITIONING *pPinit, llong *distanceInit, int clus, int initial)
{
  CreateNewCodebook(pCBnew, clus, pTS);
  CreateNewPartitioning(pPnew, pTS, clus);
  CreateNewCodebook(pCBinit, clus, pTS);
  CreateNewPartitioning(pPinit, pTS, clus);

  if (initial)
    {
    if (initial == 1)
      {
      GenerateOptimalPartitioningGeneral(pTS, pCB, pP, MSE);
      }
    else if (initial == 2)
      {
      GenerateOptimalCodebookGeneral(pTS, pCB, pP, MSE);
      }

    CopyCodebook(pCB, pCBinit);
    CopyPartitioning(pP, pPinit);
    CalculateDistances(pTS, pCBinit, pPinit, distanceInit);
  }
}


/* ------------------------------------------------------------------ */


void GenerateSolution(TRAININGSET *pTS, CODEBOOK *pCBnew,
			     PARTITIONING *pPnew, CODEBOOK *pCBinit,
			     PARTITIONING *pPinit, llong *distance,
			     llong *distanceInit, int InitMethod, int InitMethodVar, int initial)
{
  int i;


  if (initial)
    {
    CopyCodebook(pCBinit, pCBnew);
    CopyPartitioning(pPinit, pPnew);
    for (i = 0; i < BookSize(pTS); i++)
      {
      distance[i] = distanceInit[i];
      }
    }
  else
    {
    switch( InitMethod )
      {
          case 0:  DensityInitialCentroids(pTS, pCBnew,InitMethodVar);
          //case 0:  decreaseOverlapInit(pTS, pCBnew,InitMethodVar);
               break;
      case 1:  SelectRandomRepresentatives(pTS, pCBnew);
               break;
      case 2:  SelectRandomWeightedRepresentatives(pTS, pCBnew, InitMethodVar);
               break;
      case 3:  BradleyInitialCentroids(pTS, pCBnew,InitMethodVar);
               break;
      case 4:  ProjectionInitialCentroids(pTS, pCBnew,InitMethodVar);
               break;
      case 5:  LuxburgInitialCentroids(pTS, pCBnew);
               break;
      case 6:  MaxMinInitialCentroids(pTS, pCBnew, InitMethodVar);
               break;
      case 7:
               RandomPartitionInitialCentroids(pTS, pCBnew, InitMethodVar);
               break;
      case 8:
               hierarchicalSplittingCodebook(pTS, pCBnew,InitMethodVar);
               break;
      case 9:
               SortingHeuristic(pTS, pCBnew,InitMethodVar);
               break;

      case 101:  DensitySampledCentroids(pTS, pCBnew,InitMethodVar);
                 break;

      default: break;
      }

    GenerateOptimalPartitioningGeneral(pTS, pCBnew, pPnew, MSE);
    CalculateDistances(pTS, pCBnew, pPnew, distance);
    }
}


/* ------------------------------------------------------------------ */


//#define VISUALIZER 1

void KMeansIterate(TRAININGSET *pTS, CODEBOOK *pCBnew,
PARTITIONING *pPnew, llong *distance, int quietLevel, int i, int *iter,
double time, double *error, int initial, int MaxIter)
{
  double starttime = GetClock(time);

  int      activeCount;
  int      active[BookSize(pCBnew)];
  llong    cdist[BookSize(pCBnew)];
  double   oldError;

  *iter  = 0;
  *error = AverageErrorForSolution(pTS, pCBnew, pPnew, MSE);

#ifdef VISUALIZER
  PrintIterationKM(quietLevel, i, *iter, *error, GetClock(time),pCBnew);
#else
  PrintIterationKM(quietLevel, i, *iter, *error, GetClock(time));
#endif

  double inittime = GetClock(time) - starttime;

  /* Perform K-means iterations */
  do
    {
    (*iter)++;
    oldError = *error;

    OptimalRepresentatives(pPnew, pTS, pCBnew, active, cdist, &activeCount, NO);
    OptimalPartition(pCBnew, pTS, pPnew, active, cdist, activeCount, distance, NO);

    *error = AverageErrorForSolution(pTS, pCBnew, pPnew, MSE);
#ifdef VISUALIZER
    PrintIterationKM(quietLevel, i, *iter, *error, GetClock(time),pCBnew);
#else
    PrintIterationKM(quietLevel, i, *iter, *error, GetClock(time));
#endif

    }
  /* until no improvement or MaxIter reached */
  while (*error < oldError && !(MaxIter > 0 && *iter >= MaxIter));

#ifndef VISUALIZER

  /* Printing debug information */
  if ((quietLevel >= 4) && (*iter) > 0)
     {
     PrintIterationKMSummary(GetClock(time)-starttime, inittime);
     }
#endif

}


/* -------------------------------------------------------------------- */
/* Calculates data objects current distances to their cluster centroids */
/* -------------------------------------------------------------------- */


void CalculateDistances(TRAININGSET *pTS, CODEBOOK *pCB, PARTITIONING *pP, llong *distance)
{
  int i, j;

  for (i = 0; i < BookSize(pTS); i++)
    {
    j = Map(pP, i);
    distance[i] = VectorDistance(Vector(pTS, i), Vector(pCB, j), VectorSize(pTS), MAXLLONG, EUCLIDEANSQ);
    }
}


/*-------------------------------------------------------------------*/
// AKTIVITEETIN PÄIVITTÄMINEN TULEE TÄNNE

/* generates optimal codebook with respect to a given partitioning */
void OptimalRepresentatives(PARTITIONING *pP, TRAININGSET *pTS, CODEBOOK *pCB,
int *active, llong *cdist, int *activeCount, int travellerSearch)
{
  int i, j;
  VECTORTYPE v;

  j = 0;
  v = CreateEmptyVector(VectorSize(pCB));

  for(i = 0; i < BookSize(pCB); i++)
    {
    if (CCFreq(pP, i) > 0)
      {
      CopyVector(Vector(pCB, i), v, VectorSize(pCB));
      /* calculate mean values for centroid */
      PartitionCentroid(pP, i, &Node(pCB, i));
      /* if centroid changed, cluster is active */
      if (CompareVectors(Vector(pCB, i), v, VectorSize(pCB)) != 0)
        {
	if(travellerSearch)
          {
	  /* calculate the distance centroid moved */
          cdist[i] = VectorDistance(v, Vector(pCB, i), VectorSize(pTS), MAXLLONG, EUCLIDEANSQ);
          }
        active[j] = i;
        j++;
        }
      }
    else
      {
      VectorFreq(pCB, i) = 0;
      }
    }

  FreeVector(v);
  (*activeCount) = j;
}

/*-------------------------------------------------------------------*/
/* generates optimal partitioning with respect to a given codebook */
// AKTIIVINEN-PASIIVINEN VEKTORI MUUTOS


void OptimalPartition(CODEBOOK *pCB, TRAININGSET *pTS, PARTITIONING *pP,
int *active, llong *cdist, int activeCount, llong *distance, int quietLevel)
{
  int i, j, k;
  int nearest;
  llong error, dist;
  CODEBOOK CBact;

  if (quietLevel >= 5)  PrintMessage("\n Optimal Partition starts. ActiveCount=%i..\n", activeCount);

  /* all vectors are static; there is nothing to do! */
  if (activeCount < 1) return;

  /* creating subcodebook (active clusters) */
  if (quietLevel >= 5)  PrintMessage("Creating subcodebook...");
  CreateNewCodebook(&CBact, activeCount, pTS);
  for (i = 0; i < activeCount; i++)
    {
    CopyVector(Vector(pCB, active[i]), Vector(&CBact, i), VectorSize(pCB));
    }
  if (quietLevel >= 5)  PrintMessage("Done.\n");

  if (quietLevel >= 5)  PrintMessage("Looping ... ");
  for(i = 0; i < BookSize(pTS); i++)
     {
     if (quietLevel >= 5)  PrintMessage(" %i ", i);
     j     = Map(pP, i);
     k     = BinarySearch(active, activeCount, j);
     dist  = VectorDistance(Vector(pTS, i), Vector(pCB, j), VectorSize(pTS), MAXLLONG, EUCLIDEANSQ);

     // static vector - search subcodebook
     if (k < 0)
       {
       nearest = FindNearestVector(&Node(pTS,i), &CBact, &error, 0, EUCLIDEANSQ);
       nearest = (error < dist) ? active[nearest] : j;
       }
     // active vector, centroid moved closer - search subcodebook
     else if (dist < distance[i])
       {
       nearest = FindNearestVector(&Node(pTS,i), &CBact, &error, k, EUCLIDEANSQ);
       nearest = active[nearest];
       }
     // active vector, centroid moved farther - FULL search
     else
       {
       nearest = FindNearestVector(&Node(pTS,i), pCB, &error, j, EUCLIDEANSQ);
       }

     if (nearest != j)
       {
       /* closer cluster was found */
       ChangePartition(pTS, pP, nearest, i);
       distance[i] = error;
       }
     else
       {
       distance[i] = dist;
       }
    }

  FreeCodebook(&CBact);

  if (quietLevel >= 5)  PrintMessage("Optimal Partition ended.\n");
}

/*-------------------------------------------------------------------*/
/* arr must be sorted ascending order! */


int BinarySearch(int *arr, int size, int key)
{
  int top, bottom, middle;

  top = 0;
  bottom = size - 1;
  middle = (top + bottom) / 2;

  do
    {
    if (arr[middle] < key)     top    = middle + 1;
    else                       bottom = middle;
    middle = (top + bottom) / 2;
    }
  while (top < bottom);

  if (arr[middle] == key)    return middle;
  else                       return -1;
}


