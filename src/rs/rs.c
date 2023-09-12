/*-------------------------------------------------------------------*/
/* RS.C         Marko Tuononen, Pasi Franti                          */
/*                                                                   */
/* Advanced model implementation of Random Swap (RS)                 */
/*                                                                   */
/*   "Randomized local search algorithm for the clustering problem"  */
/*   Pattern Analysis and Applications, 3 (4), 358-369, 2000.        */
/*   Pasi Franti and Juha Kivijarvi                                  */
/*                                                                   */
/* This model also includes possibility to select deterministicly,   */
/* which cluster centroid to swap. In case of deterministic choice   */
/* in every iteration cluster that increases objective function      */
/* (MSE) least, if removed, is selected.                             */
/*                                                                   */
/* K-means -operation uses activity detection presented in           */
/*                                                                   */
/*   "A fast exact GLA based on code vector activity detection"      */
/*   IEEE Trans. on Image Processing, 9 (8), 1337-1342, August 2000. */
/*   Timo Kaukoranta, Pasi Franti and Olli Nevalainen                */
/*                                                                   */
/* Naming conventions used in the code                               */
/*                                                                   */
/*    TS        training set (data objects)                          */
/*    CB        codebook (cluster representatives, centroids)        */
/*    P         partitioning (pointing from TS to CB)                */
/*                                                                   */
/*    p-prefix  pointer, e.g. pTS is pointer to the training set TS  */
/*                                                                   */
/* ----------------------------------------------------------------- */
/*                                                                   */
/* Traveller search mode implementation is based on idea presented   */
/* in following paper:                                               */
/*                                                                   */
/*  "Faster and more robust point symmetry-based K-means algorithm"  */
/*  Pattern Recognition, 40, 410-422, 2007.                          */
/*  Kuo-Liang Chung, Jhin-Sian Lin                                   */
/*                                                                   */
/* The main idea of the mode is that the optimal centroid search for */
/* the closer data points is filtered to include only the current    */
/* centroid and all active centroids that have moved greater distance*/
/* and the current centroid.					     */
/*                                                                   */
/* ----------------------------------------------------------------- */
/*                                                                   */
/* HISTORY:                                                          */
/*                                                                   */
/* 0.34 PF  Added: RandomCodebook, LuxburgInitialization (20.9.16)   */
/* 0.33 PF  MonitorProgress mode + Stochastic variant (8.7.16)       */
/* 0.31 PF  Renamed RLS->RS; Refactoring code; CI-index (25.6.16)    */
/* 0.25 AH  Traveller search and new Q-level 4 features (28.2.10)    */
/* 0.24 MM  Correct random initialization (15.7.09)                  */
/* 0.22 VH  Correct Random initialization (15.4.09)                  */
/* 0.21 MT  Modified RSInfo() to print less information              */
/* 0.20 MT  Fixed SelectRandomDataObject, added automatic iter.count */
/*-------------------------------------------------------------------*/


#define ProgName       "RS"
#define VersionNumber  "Version 0.34"
#define LastUpdated    "20.9.2016"  /* PF */

/* converts ObjectiveFunction values to MSE values */
#define CALC_MSE(val) (double) (val) / (TotalFreq(pTS) * VectorSize(pTS))

#define AUTOMATIC_MAX_ITER  50000
#define AUTOMATIC_MIN_SPEED 1e-5
#define min(a,b) ((a) < (b) ? (a) : (b))

/*-------------------------------------------------------------------*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <float.h>
#include <string.h>
#include <assert.h>

#include "cb.h"
#include "random.h"
#include "interfc.h"
#include "reporting.h"
#include "limits.h"

#include "memctrl.h"
//#include "graph.h"
//#include "fastxnn.h"


#include "kmeans.h"


/* ========================== PROTOTYPES ============================= */

int PerformRS(TRAININGSET *pTS, CODEBOOK *pCB, PARTITIONING *pP, int iter,
    int kmIter, int deterministic, int travellerSearch, int quietLevel,
    int useInitialCB, int monitoring);
void InitializeSolution(PARTITIONING *pP, CODEBOOK *pCB, TRAININGSET *pTS,
    int clus);
void FreeSolution(PARTITIONING *pP, CODEBOOK *pCB);
YESNO StopCondition(double currError, double newError, int iter);
double GenerateInitialSolution(PARTITIONING *pP, CODEBOOK *pCB,
    TRAININGSET *pTS, int useInitialCB);
void SelectRandomRepresentatives(TRAININGSET *pTS, CODEBOOK *pCB);
void SelectRandomRepresentatives2(TRAININGSET *pTS, CODEBOOK *pCB); /* mm */
int SelectRandomDataObject(CODEBOOK *pCB, TRAININGSET *pTS);
void RandomCodebook(TRAININGSET *pTS, CODEBOOK *pCB);
void RandomSwap(CODEBOOK *pCB, TRAININGSET *pTS, int *j, int deterministic,
    int quietLevel);
void LocalRepartition(PARTITIONING *pP, CODEBOOK *pCB, TRAININGSET *pTS,
    int j, double time, int quietLevel);
void OptimalPartition(CODEBOOK *pCB, TRAININGSET *pTS, PARTITIONING *pP,
    int *active, llong *cdist, int activeCount, llong *distance, int quietLevel);
void KMeans(PARTITIONING *pP, CODEBOOK *pCB, TRAININGSET *pTS,
    llong *distance, int iter, int travellerSearch, int quietLevel, double time);
double ObjectiveFunction(PARTITIONING *pP, CODEBOOK *pCB, TRAININGSET *pTS);
int FindSecondNearestVector(BOOKNODE *node, CODEBOOK *pCB, int firstIndex,
    llong *secondError);
int SelectClusterToBeSwapped(TRAININGSET *pTS, CODEBOOK *pCB,
    PARTITIONING *pP, llong *distance);
char* RSInfo(void);


/* ========================== FUNCTIONS ============================== */

/* Gets training set pTS (and optionally initial codebook pCB or
   partitioning pP) as a parameter, generates solution (codebook pCB +
   partitioning pP) and returns 0 if clustering completed successfully.
   N.B. Random number generator (in random.c) must be initialized! */

int PerformRS(TRAININGSET *pTS, CODEBOOK *pCB, PARTITIONING *pP, int iter,
int kmIter, int deterministic, int travellerSearch, int quietLevel,
int useInitial, int monitoring)
{
  PARTITIONING  Pnew;
  CODEBOOK      CBnew, CBref;
  int           i, j, better;
  int           ci=0, ciPrev=0, ciZero=0, ciMax=0, PrevSuccess=0;
  int           CIHistogram[111];
  double        currError, newError;
  llong         distance[BookSize(pTS)];
  double        c, error;
  int           stop=NO, automatic=((iter==0) ? YES : NO);

  /* Error checking for invalid parameters */
  if ((iter < 0) || (kmIter < 0) || (BookSize(pTS) < BookSize(pCB)))
    {
    return 1;  // Error: clustering failed
    }

  /* Progress monitor uses input codebook as reference */
  if (monitoring)
    {
    CreateNewCodebook(&CBref, BookSize(pCB), pTS);
    CopyCodebook(pCB, &CBref);
    for( ci=0; ci<=100; ci++ ) CIHistogram[ci]=0;
    useInitial *= 100;  /* Special code: 0->0, 1->100, 2->200 */
    }

  InitializeSolution(&Pnew, &CBnew, pTS, BookSize(pCB));
  SetClock(&c);
  currError = GenerateInitialSolution(pP, pCB, pTS, useInitial);
  error = CALC_MSE(currError);
  if(useInitial) ciPrev = CentroidIndex(&CBref, pCB);
  else           ciPrev = 100;

  /* use automatic iteration count */
  if (automatic)  iter = AUTOMATIC_MAX_ITER;

  PrintHeader(quietLevel);
#ifndef VISUALIZER
  PrintIterationRS(quietLevel, 0, error, 0, GetClock(c), 1);
#endif

  /* Deterministic variant initialization */
  if (deterministic)
    {
    CalculateDistances(pTS, pCB, pP, distance);
    j = SelectClusterToBeSwapped(pTS, pCB, pP, distance);
    }

  /* - - - - -  Random Swap iterations - - - - - */

  for (i=1; (i<=iter) && (!stop); i++)
    {
    better = NO;

    /* generate new solution */
    CopyCodebook(pCB, &CBnew);
    CopyPartitioning(pP, &Pnew);
    RandomSwap(&CBnew, pTS, &j, deterministic, quietLevel);

    /* tuning new solution */
    LocalRepartition(&Pnew, &CBnew, pTS, j, c, quietLevel);
    KMeans(&Pnew, &CBnew, pTS, distance, kmIter, travellerSearch, quietLevel, c);

    newError = ObjectiveFunction(&Pnew, &CBnew, pTS);
    error    = CALC_MSE(newError);

    /* Found better solution */
    // printf("newError=%f currError=%f\n",newError,currError);
    if (newError < currError)
      {
      /* Check stopping criterion */
      if (automatic)
        {
        stop = StopCondition(currError, newError, i);
        }
      /* Monitoring outputs CI-value: relative to Prev or Reference */
      if(monitoring)
         {
         if(useInitial) ci = CentroidIndex(&CBnew, &CBref);
         else           ci = CentroidIndex(&CBnew, pCB);
         /* Update Success histogram */
         if( (ci>=0) && (ci<ciPrev) && (ci<100) )
           {
           /* printf("XXXX Prev=%d  Curr=%d  CI=%d  CIPrev=%i  Iter=%d \n", PrevSuccess, i, ci, ciPrev, (i-PrevSuccess)); */
           CIHistogram[ci] += (i-PrevSuccess);
           if(ci>ciMax)  ciMax = ci;
           if(ci==0)     ciZero = i;
           PrevSuccess = i;
           ciPrev      = ci;
           }
         if( (ci>ciPrev) && quietLevel ) printf("!!! CI increased %i to %i at iteration %d\n", ciPrev, ci, i);
         }

      CopyCodebook(&CBnew, pCB);
      CopyPartitioning(&Pnew, pP);

      currError = newError;
      better = YES;

      if (deterministic) /* Alterantive ro Random. But why here?  */
        {
        j = SelectClusterToBeSwapped(pTS, pCB, pP, distance);
        }
      }

#ifndef VISUALIZER
    PrintIterationRS(quietLevel, i, error, ci, GetClock(c), better);
#endif
    }

  /* - - - - -  Random Swap iterations - - - - - */

  error = CALC_MSE(currError);
#ifndef VISUALIZER
  PrintFooterRS(quietLevel, i-1, error, GetClock(c));
#endif

  if(monitoring)
     {
     PrintMessage("Total: %-7d   Swaps: ", ciZero);
     for( ci=0; ci<=ciMax; ci++ )
       {
       PrintMessage("%3d  ", CIHistogram[ci]);
       }
     PrintMessage("\n", ciZero);
     }

  FreeSolution(&Pnew, &CBnew);
  return 0;
}


/*-------------------------------------------------------------------*/


void InitializeSolution(PARTITIONING *pP, CODEBOOK *pCB, TRAININGSET *pTS,
int clus)
{
  CreateNewCodebook(pCB, clus, pTS);
  CreateNewPartitioning(pP, pTS, clus);
}


/*-------------------------------------------------------------------*/


void FreeSolution(PARTITIONING *pP, CODEBOOK *pCB)
{
  FreeCodebook(pCB);
  FreePartitioning(pP);
}


/*-------------------------------------------------------------------*/


YESNO  StopCondition(double currError, double newError, int iter)
{
  static double   currImpr=DBL_MAX, prevImpr=DBL_MAX;
  static int      prevIter=1;

  currImpr  = (double)(currError - newError) / (double)currError;
  currImpr /= (double) (iter - prevIter);
  if (AUTOMATIC_MIN_SPEED < currImpr + prevImpr)
     {
     prevImpr = currImpr;
     prevIter = iter;
     return(NO);
     }
  else  /* too slow speed, better to stop.. */
     {
     return(YES);
     }
}


/*-------------------------------------------------------------------*/


double GenerateInitialSolution(PARTITIONING *pP, CODEBOOK *pCB,
TRAININGSET *pTS, int useInitial)
{
  if (useInitial == 1)
  {
    GenerateOptimalPartitioningGeneral(pTS, pCB, pP, MSE);
  }
  else if (useInitial == 2)
  {
    GenerateOptimalCodebookGeneral(pTS, pCB, pP, MSE);
  }
  else
  {
    SelectRandomRepresentatives(pTS, pCB);
    GenerateOptimalPartitioningGeneral(pTS, pCB, pP, MSE);
  }

  return ObjectiveFunction(pP, pCB, pTS);
}



/*-------------------------------------------------------------------*/


int SelectRandomDataObject(CODEBOOK *pCB, TRAININGSET *pTS)
{
  int i, j, count = 0;
  int ok;

  do
    {
    count++;

    /* random number generator must be initialized! */
    j = IRZ(BookSize(pTS));

    /* eliminate duplicates */
    ok = 1;
    for (i = 0; i < BookSize(pCB); i++)
      {
      if (EqualVectors(Vector(pCB, i), Vector(pTS, j), VectorSize(pTS)))
        {
        ok = 0;
        }
      }
  }
  while (!ok && (count <= BookSize(pTS)));   /* fixed 25.01.2005 */

  return j;
}


/*-------------------------------------------------------------------*/
/* random number generator must be initialized! */


void RandomSwap(CODEBOOK *pCB, TRAININGSET *pTS, int *j, int deterministic,
                int quietLevel)
{
  int i;

  if (!deterministic)
    {
    *j = IRZ(BookSize(pCB));
    }

  i = SelectRandomDataObject(pCB, pTS);

  CopyVector(Vector(pTS, i), Vector(pCB, *j), VectorSize(pTS));
  if (quietLevel >= 5)  PrintMessage("Random Swap done: x=%i  c=%i \n", i, *j);
}



/*-------------------------------------------------------------------*/


void LocalRepartition(PARTITIONING *pP, CODEBOOK *pCB, TRAININGSET *pTS, int j,
double time, int quietLevel)
{
  if (quietLevel >= 5)  PrintMessage("Local repartition of vector %i \n", j);

  /* object rejection; maps points from a cluster to their nearest cluster */
  LocalRepartitioningGeneral(pTS, pCB, pP, j, EUCLIDEANSQ);

  /* object attraction; moves vectors from their old partitions to
     a the cluster j if its centroid is closer */
  RepartitionDueToNewVectorGeneral(pTS, pCB, pP, j, EUCLIDEANSQ);

  if (quietLevel >= 3)  PrintMessage("RepartitionTime= %f   ", GetClock(time));
}





/*-------------------------------------------------------------------*/
/* generates optimal partitioning with respect to a given codebook */
// AKTIIVINEN-PASIIVINEN VEKTORI MUUTOS


void OptimalPartitionTraveller(CODEBOOK *pCB, TRAININGSET *pTS, PARTITIONING *pP,
int *active, llong *cdist, int activeCount, llong *distance, int quietLevel)
{
  int i, j, k, l;
  int nearest;
  llong error, dist;
  CODEBOOK CBact;

  if (quietLevel >= 5)  PrintMessage("\n Optimal Partition (Traveller variant) starts. ActiveCount=%i..\n", activeCount);

  /* all vectors are static; there is nothing to do */
  if (activeCount < 1)  return;

  /* traveller codebook variables */
  CODEBOOK CBactTrvs[activeCount];
  int trvCount[activeCount];
  int trvArray[activeCount][BookSize(pCB)];
  int ptrv[BookSize(pCB)];

  if (quietLevel >= 5)  PrintMessage("Traveller arrays created.\n");

  /* creating subcodebook (active clusters) */
  if (quietLevel >= 5)  PrintMessage("Creating subcodebook...");
  CreateNewCodebook(&CBact, activeCount, pTS);
  for (i = 0; i < activeCount; i++) {
    CopyVector(Vector(pCB, active[i]), Vector(&CBact, i), VectorSize(pCB));
  if (quietLevel >= 5)  PrintMessage("Done.\n");

  /* Creating traveller codebooks */
  for (i = 0; i < activeCount; i++)
      {
      /* find centroids among active centroieds whose move distance is greater than current centroids */
      trvCount[i] = 0;
      for (l = 0; l < activeCount; l++)
          {
          if(active[i] == active[l] || cdist[active[i]] < cdist[active[l]])
              {
              trvArray[i][trvCount[i]++] = active[l];
              }
          }
      /* create traveller codebook */
      CreateNewCodebook(&CBactTrvs[i], trvCount[i], pTS);
      for (l = 0; l < trvCount[i]; l++)
          {
          CopyVector(Vector(pCB, trvArray[i][l]), Vector(&CBactTrvs[i], l), VectorSize(pCB));
          }
      ptrv[active[i]] = i;
      }

  if (quietLevel >= 5)  PrintMessage("Looping ... ");
  for(i = 0; i < BookSize(pTS); i++)
     {
     if (quietLevel >= 5)  PrintMessage(" %i ", i);
     j     = Map(pP, i);
     k     = BinarySearch(active, activeCount, j);
     dist  = VectorDistance(Vector(pTS, i), Vector(pCB, j), VectorSize(pTS), MAXLLONG, EUCLIDEANSQ);

     if (k < 0)  /* static vector */
       {
       // search subcodebook
       nearest = FindNearestVector(&Node(pTS,i), &CBact, &error, 0, EUCLIDEANSQ);
       nearest = (error < dist) ? active[nearest] : j;
       }
     else if (dist < distance[i])  /* active vector, centroid moved closer */
       {
       // search traveller codebook
       k = BinarySearch(trvArray[ptrv[j]], trvCount[ptrv[j]], j);
       nearest = FindNearestVector(&Node(pTS,i), &CBactTrvs[ptrv[j]], &error, k, EUCLIDEANSQ);
       nearest = trvArray[ptrv[j]][nearest];
       }

     else  /* active vector, centroid moved farther */
       {
       // search full codebook
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

  /* free traveller codebooks */
  for(i = 0; i < activeCount; i++)
     {
     FreeCodebook(&CBactTrvs[i]);
     }
  }

  if (quietLevel >= 5)  PrintMessage("Optimal Partition Traveller variant ended.\n");
}


/*-------------------------------------------------------------------*/
/* fast K-means implementation (uses activity detection method) */


void KMeans(PARTITIONING *pP, CODEBOOK *pCB, TRAININGSET *pTS, llong *distance,
int iter, int travellerSearch, int quietLevel, double time)
{

  double starttime = GetClock(time);

  int     i, activeCount;
  int     active[BookSize(pCB)];
  llong   cdist[BookSize(pCB)];

  CalculateDistances(pTS, pCB, pP, distance);

  double inittime = GetClock(time) - starttime;

  /* performs iter K-means iterations */
  for (i = 0; i < iter; i++)
    {
    /* OptimalRepresentatives-operation should be before
       OptimalPartition-operation, because we have previously tuned
       partition with LocalRepartition-operation */
    OptimalRepresentatives(pP, pTS, pCB, active, cdist, &activeCount, travellerSearch);
    if(travellerSearch)
       {
       OptimalPartitionTraveller(pCB, pTS, pP, active, cdist, activeCount, distance, quietLevel);
       }
    else
       {
       OptimalPartition(pCB, pTS, pP, active, cdist, activeCount, distance, quietLevel);
       }

    if (quietLevel >= 3)
      {
#ifndef VISUALIZER
      PrintIterationActivity(GetClock(time), i, activeCount, BookSize(pCB), quietLevel);
#endif
      }
    }

#ifndef VISUALIZER
  if ((quietLevel >= 4) && iter > 0)
     {
     PrintIterationKMSummary(GetClock(time)-starttime, inittime);
     }
#endif
}


/*-------------------------------------------------------------------*/


double ObjectiveFunction(PARTITIONING *pP, CODEBOOK *pCB, TRAININGSET *pTS)
{
  double sum = 0.0;
  int i, j;

  /* sum of squared distances of the data object to their
     cluster representatives */
  for (i = 0; i < BookSize(pTS); i++)
    {
    j = Map(pP, i);
    sum += VectorDistance(Vector(pTS, i), Vector(pCB, j),
           VectorSize(pTS), MAXLLONG, EUCLIDEANSQ) * VectorFreq(pTS, i);
    }

// printf("SUM=%f\n",sum);
  return sum;
}



/*-------------------------------------------------------------------*/


int FindSecondNearestVector(BOOKNODE *node, CODEBOOK *pCB,
                            int firstIndex, llong *secondError)
{
  int   i;
  int   secondIndex;
  llong e;

  secondIndex = -1;
  *secondError = MAXLLONG;

  for(i = 0; i < BookSize(pCB); i++)
    {
    e = VectorDistance(Vector(pCB,i), node->vector, VectorSize(pCB),
        *secondError, EUCLIDEANSQ);

      if ((e < *secondError) && (i != firstIndex))
    {
      *secondError = e;
      secondIndex  = i;
      }
    }
  return secondIndex;
}


/*-------------------------------------------------------------------*/
/* selects deterministicly, which cluster centroid to swap. one that
   increases objective function (MSE) least, if removed, is selected. */

int SelectClusterToBeSwapped(TRAININGSET *pTS, CODEBOOK *pCB,
                             PARTITIONING *pP, llong *distance)
{
  int i, j, k, min;
  llong error;
  llong priError[BookSize(pCB)];  /* current error; data objects are in
                                     their primary (closest) cluster) */
  llong secError[BookSize(pCB)];  /* error after partition is removed and
                                     data objects are repartitioned; data
                                     objects are in their secondary
                                     (second closest) cluster */

  /* initializing */
  for (i = 0; i < BookSize(pCB); i++)
    {
    priError[i] = 0;
    secError[i] = 0;
    }

  /* calculating primary and secondary cluster errors */
  for (i = 0; i < BookSize(pTS); i++)
    {
    j = Map(pP, i);
    k = FindSecondNearestVector(&Node(pTS,i), pCB, j, &error);
    /* k will not be stored, only the return error value is used */

    priError[j] += distance[i] * VectorFreq(pTS, i);
    secError[j] += error * VectorFreq(pTS, i);
    }

  /* finding cluster that increases objective function least */
  min = -1;
  error = MAXLLONG;
  for (j = 0; j < BookSize(pCB); j++)
    {
    if ((secError[j] - priError[j]) < error)
      {
      min = j;
      error = secError[j] - priError[j];
      }
    }

  return min;
}


/*-------------------------------------------------------------------*/


char* RSInfo(void)
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


/*-------------------------------------------------------------------*/
