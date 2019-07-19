/*-------------------------------------------------------------------*/
/* SPLIT.C        Timo Kaukoranta                                    */
/*                Juha Kivijärvi (GenerateSearchTree)                */
/*                                                                   */
/* - Module of Splitting method                                      */
/*                                                                   */
/*-------------------------------------------------------------------*/

#define ProgName       "SPLIT"
#define VersionNumber  "V.0.38"
#define LastUpdated    "28.04.2009" /* PF */

/* ----------------------------------------------------------------- */
/*      Changelog:                                                   */
/*                                                                   */
/* PF   0.38  Combined updates from different versions               */
/* JS   0.37  Local repartition speedup                              */
/*            Compute centroids as the new code vectors after split. */
/*            Remapping: Integrate faster partition switching.       */
/*            Bugfix in FurthestVectorInPartition: long -> llong     */
/*            Changed: MAXDOUBLE -> DBL_MAX, MAXLONG -> LONG_MAX     */
/* PF   0.36  added: FasterVariant using ChangePartitionFast         */
/*                                                                   */
/* ----------------------------------------------------------------- */

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <float.h>
#include <limits.h>

#include "bintree.h"
#include "cb.h"
#include "memctrl.h"
#include "owntypes.h"
#include "random.h"
#include "sa.h"
#include "sort.h"
#include "split.h"



/* ------------------- Parameters for the Split module -------------- */

static struct { int ShowProgress;
                int PartitionErrorType;
                int NewVectors;
                int NewVectorsHeuristically;
                int Hyperplane;
                int HyperplanePivot;
                int HyperplaneLine;
                int PartitionRemapping;
                int LocalGLA;
                int LocalGLAIterations;
                int GLAIterations;
              } ModuleParameters;

#define Value(x)      (ModuleParameters.x)


/* These enum definitions are from the *.fac file. */
enum { SplitAllPartitions, SquareError, SizeOfPartition, FurthestTwoVectors, ThirdMoment, SplitRandom, OptimalSelection };
enum { CurrentAndSigma, CurrentPlusMinusSigma, CurrentAndFurthest, nv4, CurrentAndRandom, TwoFurthest, TwoRandom, MeansOfCurrentAndFurthest, FasterVariant };
enum { HPPCentroid, HPPOptimal, HPPOptimalEquitz, HPPWeightedCentroid, HPPGoeddelBass, HPPMomentPreserving, HPPLloydScalarQuantization, HPPOptimalWu };
enum { HPLPCA, HPLRegression, HPLFurthest };


/*-----------------------------  M i s c  ----------------------------*/


typedef struct { int        whoami;
                 llong      DistortionValue;
                 llong      SquareError;
                 llong      FurthestDistance;
                 llong      ThirdMoment;
                 llong      DistortionChange;
                 int        TrVector1;
                 int        TrVector2;
                 VECTORTYPE Centroid1;
                 VECTORTYPE Centroid2;
               } CLUSTERINFO;

typedef int    PREPARTITIONING;

typedef struct { int    vindex;
                 double distance;
               } PROJECTION;

#define round(a)      ((a) + 0.5)
#define fpositive(a)  ((a) > FLT_EPSILON)
#define printfe       printf("File %s line %i ", __FILE__, __LINE__);printf
#define sqr(a)        ((a) * (a))
#define absd(a)       ((a) > 0 ? (a) : -(a))

#define MAXERROR      ((LONG_MAX) >> 2)


/* =================== Parameters for the Split module ============= */


void SetSplitParameters(int ShowProgress,
                        int PartitionErrorType,
                        int NewVectors,
                        int NewVectorsHeuristically,
                        int Hyperplane,
                        int HyperplanePivot,
                        int HyperplaneLine,
                        int PartitionRemapping,
                        int LocalGLA,
                        int LocalGLAIterations,
                        int GLAIterations)
{
  Value(ShowProgress)            = ShowProgress;
  Value(PartitionErrorType)      = PartitionErrorType;
  Value(NewVectors)              = NewVectors;
  Value(NewVectorsHeuristically) = NewVectorsHeuristically;
  Value(Hyperplane)              = Hyperplane;
  Value(HyperplanePivot)         = HyperplanePivot;
  Value(HyperplaneLine)          = HyperplaneLine;
  Value(PartitionRemapping)      = PartitionRemapping;
  Value(LocalGLA)                = LocalGLA;
  Value(LocalGLAIterations)      = LocalGLAIterations;
  Value(GLAIterations)           = GLAIterations;
}


/*===================  PRINTING FOR DEBUGGING  ======================*/


static void PrintPartition(TRAININGSET* TS, PARTITIONING* P, int Pindex)
{
  int j;

  printf("PP%5i: Freq=%5i Unique=%5i\n",
         Pindex, CCFreq(P,Pindex), UniqueVectors(P, Pindex));
  for( j = FirstVector(P,Pindex); ! EndOfPartition(j); j = NextVector(P,j) )
    {
    printf("j=%5i ", j);
    PrintVector(Vector(TS, j), VectorSize(TS), 1);
    }
  printf("\n");
}


/* ----------------------------------------------------------------- */


static void PrintBintree(BINTREE* Tree)
{
  STACK s;
  CLUSTERINFO* c;
  llong prevdist = -LONG_MAX;

  printf("Treesize=%i\n", Tree->nnodes);
  InitInOrderBintree(Tree, &s);
  while( (c = (CLUSTERINFO*)InOrderBintree(&s)) != NULL )
    {
    printf("whoami=%4i DistValue=%10.1f\n",
           c->whoami, (double)c->DistortionValue);
    if( prevdist > c->DistortionValue )
      {
      printf("order broken");
      exit( -1 );
      }
    prevdist = c->DistortionValue;
    }
}


/* =========================  COUNTERS  ============================ */


static llong* AllocateCounterVector(int Vsize)
{
  int k;
  llong* counter = (llong*)allocate( Vsize * sizeof(llong) );

  for( k = 0; k < Vsize; k++ )
    {
    counter[k] = 0LL;
    }
  return( counter );
}


/* ----------------------------------------------------------------- */


static void AddCounterToVector(llong*     counter,
                               VECTORTYPE vector,
                               int        Vsize,
                               int        maxelement)
{
  int k;

  for( k = 0; k < Vsize; k++ )
    {
    vector[k] += (VECTORELEMENT)counter[k];
    if( vector[k] > maxelement )
      {
      vector[k] = maxelement;
      }
    }
}


/* ----------------------------------------------------------------- */


static void SubtractCounterFromVector(llong*     counter,
                                      VECTORTYPE vector,
                                      int        Vsize,
                                      int        minelement)
{
  int k;

  for( k = 0; k < Vsize; k++ )
    {
    vector[k] -= (VECTORELEMENT)counter[k];
    if( vector[k] < minelement )
      {
      vector[k] = minelement;
      }
    }
}


/* ----------------------------------------------------------------- */



static void MeanOfTwoVectors(VECTORTYPE source1,
                             VECTORTYPE source2,
                             VECTORTYPE dest,
                             int        Vsize)
{
  int  k;

  for( k = 0; k < Vsize; k++ )
    {
    dest[k] = (VECTORELEMENT) ((llong)source1[k] + (llong)source2[k]) / 2;
    }
}


/* ================================================================= */


static long MSEDifference(BOOKNODE *v1,
                          BOOKNODE *v2,
                          int       Vsize,
                          llong     MinError)
{
  int  k;
  llong e;
  llong Error = 0;

  for( k = 0; k < Vsize; k++ )
    {
    e      = ((llong)v1->vector[k] - (llong)v2->vector[k]);
    Error += sqr(e);
    if( Error >= MinError )
      {
      break;
      }
    }

  return( Error );
}


/* ----------------------------------------------------------------- */


static llong ThirdMomentDifference(llong cent, llong x)
{
  long e = cent - x;

  if( e >= 0 )
    {
    return( sqr(e) );
    }
  else
    {
    return( - sqr(e) );
    }
}


/*-------------------------------------------------------------------*/


static void NearestVectorInPartition(TRAININGSET*  TS,
                                     PARTITIONING* P,
                                     int           Pindex,
                                     BOOKNODE*     V,
                                     int*          NearestV)
{
  int   j;
  llong error;
  llong MinError = MAXLLONG;

  for( j = FirstVector(P,Pindex); ! EndOfPartition(j); j = NextVector(P,j) )
    {
    error = VectorDistance(Vector(TS,j), V->vector, VectorSize(TS),
                           MinError, EUCLIDEANSQ);
    if( error < MinError )
      {
      MinError  = error;
      *NearestV = j;
      }
    }
}


/* ----------------------------------------------------------------- */


static void FurthestVectorInPartition(TRAININGSET*  TS,
                                      PARTITIONING* P,
                                      int           Pindex,
                                      BOOKNODE*     V,
                                      int*          FurthestV)
{
  int  j;
  llong error;
  llong MaxError = -1L;

  for( j = FirstVector(P,Pindex); ! EndOfPartition(j); j = NextVector(P,j) )
    {
    error = VectorDist(Vector(TS, j), V->vector, VectorSize(TS));
    if( error > MaxError )
      {
      MaxError = error;
      *FurthestV = j;
      }

    }
}


/* ----------------------------------------------------------------- */


static llong TwoFurthestVectorsInPartition(TRAININGSET*  TS,
                                           CODEBOOK*     CB,
                                           PARTITIONING* P,
                                           int           Pindex,
                                           int*          F1,
                                           int*          F2)
/* Returns the distance. */
{
  int   j, j2;
  int   count = 0;
  llong dist, distance = -1LL;

  if( Value(NewVectorsHeuristically) == YES )
    {
    /* From current codevector... */
    FurthestVectorInPartition(TS, P, Pindex, &Node(CB, Pindex), F1);
    /* ...and from that again. */
    FurthestVectorInPartition(TS, P, Pindex, &Node(TS, *F1), F2);
    distance = VectorDist(Vector(TS, *F1), Vector(TS, *F2), VectorSize(TS));
    }
  else
    {
    *F1 = ENDPARTITION;
    *F2 = ENDPARTITION;
    dist = -1LL;
    for( j = FirstVector(P,Pindex); ! IsLast(P,j); j = NextVector(P,j) )
      {
      for( j2 = NextVector(P,j); ! EndOfPartition(j2); j2 = NextVector(P,j2) )
        {
        dist = VectorDist(Vector(TS, j), Vector(TS, j2), VectorSize(TS));
        if( dist > distance )
          {
          distance = dist;
          *F1 = j;
          *F2 = j2;
          }
        }
      if( Value(ShowProgress) >= 5 )
        {
        printf("F2V:%4i j=%4i (%4i,%4i) d=%10.1f\n",
               ++count, j, *F1, *F2, (double)distance);
        fflush(stdout);
        }
      }
    if( Value(ShowProgress) >= 5 )
      {
      printf(" done.\n");
      }
    }

  if( Value(ShowProgress) >= 4 )
    {
    printf("F2V:(%4i,%4i) d=%.0f\n", *F1, *F2, (double)distance);
    PrintVector(Vector(TS, *F1), VectorSize(TS), 1);
    PrintVector(Vector(TS, *F2), VectorSize(TS), 1);
    fflush(stdout);
    }
  return( distance );
}


/*===================================================================*/


static llong PartitionError(TRAININGSET*  TS,
                            PARTITIONING* P,
                            int           Pindex,
                            BOOKNODE*     vec,
                            llong         MinDist)
{
  int   j;
  llong error = 0L;

  for( j = FirstVector(P, Pindex); ! EndOfPartition(j); j = NextVector(P, j) )
    {
    error += VectorDistance(vec->vector, Vector(TS, j), VectorSize(TS),
                            MinDist, EUCLIDEANSQ)
             * VectorFreq(TS, j);
    if( Value(ShowProgress) >= 5 )
      {
      printf("PartitionError: Pindex=%i j=%i error=%.1f freq=%i MD=%.1f\n",
             Pindex, j, (double)error, VectorFreq(TS, j), (double)MinDist);
      }
    if( error >= MinDist )
      {
      break;
      }
    }
  return( error );
}


/* ----------------------------------------------------------------- */


static double SE2MSE(TRAININGSET* TS, CODEBOOK* CB, CLUSTERINFO CI[])
{
  int   i;
  llong totalerror = 0LL;

  if( TS->TotalFreq > 0 )
    {
    for( i = 0; i < BookSize(CB); i++ )
      {
      totalerror += CI[i].DistortionValue;
      }

    return( totalerror / ((double)TS->TotalFreq * (double)VectorSize(TS)) );
    }
  else
    {
    printfe("SE2MSE: TS->TotalFreq = 0. \n");
    exit(-1);
    }
}


/* ----------------------------------------------------------------- */


static int RandomVectorFromPartition(TRAININGSET*  TS,
                                     PARTITIONING* P,
                                     int           Pindex,
                                     long          pfreq,
                                     int           discard)
{
  int  j;
  long picked;

  if( pfreq <= 0 )
    {
    printf("ERROR: RandomVectorFromPartition pfreq=%li\n", pfreq);
    }

  if( discard != ENDPARTITION )
    {
    pfreq -= VectorFreq(TS, discard);
    }

  picked = irand(1, pfreq);
  j = FirstVector(P, Pindex);
  if( discard != ENDPARTITION && j == discard )
    {
    j = NextVector(P, j);
    }
  picked -= VectorFreq(TS, j);
  while( picked > 0 )
    {
    j = NextVector(P, j);
    if( discard != ENDPARTITION && j == discard )
      {
      j = NextVector(P, j);
      }
    picked -= VectorFreq(TS, j);
    }

  return( j );
}


/* ----------------------------------------------------------------- */


static void PartitionSD(TRAININGSET*  TS,
                        PARTITIONING* P,
                        int           Pindex,
                        llong*        SD)
{
  int    j, k;
  llong* sqrsum;
  double mean, sqrmean;

  /* Initialize counters. */
  sqrsum = AllocateCounterVector(VectorSize(TS));
  for( k = 0; k < VectorSize(TS); k++ )
    {
    sqrsum[k] = 0LL;
    }

  for( j = FirstVector(P,Pindex); ! EndOfPartition(j); j = NextVector(P,j) )
    {
    for( k = 0; k < VectorSize(TS); k++ )
      {
      sqrsum[k] += sqr((llong)VectorScalar(TS, j, k)) * VectorFreq(TS, j);
      }
    }

  for( k = 0; k < VectorSize(TS); k++ )
    {
    mean    = (double)CCScalar(P, Pindex, k) / (double)CCFreq(P, Pindex);
    sqrmean = (double)sqrsum[k] / (double)CCFreq(P, Pindex);
    SD[k]   = (llong)round(sqrt( (double)(sqrmean - sqr(mean)) ));
    }

  deallocate(sqrsum);
}


/* ----------------------------------------------------------------- */


static llong PartitionThirdMoment(TRAININGSET*  TS,
                                  PARTITIONING* P,
                                  int           Pindex,
                                  BOOKNODE*     vec,
                                  int           Vsize,
                                  llong         MinDist)
{
  int    j, k;
  llong  errorsum = 0L;
  llong* counter;

  counter = AllocateCounterVector(VectorSize(TS));

  for( j = FirstVector(P,Pindex); ! EndOfPartition(j); j = NextVector(P,j) )
    {
    for( k = 0; k < VectorSize(TS); k++ )
      {
      counter[k] += ThirdMomentDifference(vec->vector[k], VectorScalar(TS,j,k));
      }
    }
  for( k = 0; k < VectorSize(TS); k++ )
    {
    errorsum += ( counter[k] >= 0 ? counter[k] : -counter[k] );
    }
  deallocate(counter);

  return( errorsum );
}


/* ----------------------------------------------------------------- */


static void LloydScalarQuantization(PROJECTION projection[],
                                    int        size,
                                    double*    centroid)
{
  int i;
  double sum = 0.0;
  double lowsum, highsum;
  int    nlow, nhigh;
  double lowmean, highmean, oldcentroid = 0.0;

  for( i = 0; i < size; i++ )
    {
    sum += projection[i].distance;
    }
  *centroid = sum / (double)size;

  while( *centroid != oldcentroid )
    {
    oldcentroid = *centroid;
    lowsum = highsum = 0.0;
    nlow = nhigh = 0;
    for( i = 0; i < size; i++ )
      {
      if( projection[i].distance < *centroid )
        {
        lowsum += projection[i].distance;
        nlow++;
        }
      else
        {
        highsum += projection[i].distance;
        nhigh++;
        }
      }
    lowmean  = lowsum  / (double)nlow;
    highmean = highsum / (double)nhigh;
    *centroid = (lowmean + highmean) / 2.0;

    if( Value(ShowProgress) >= 4 )
      {
      printf("LLoyd lm=%9.4f hm=%9.4f cent=%9.4f\n",
              lowmean, highmean, *centroid);
      }
    }
}


/* ----------------------------------------------------------------- */


int cmpdouble(const void* a, const void* b, const void* info)
{
  return( ((PROJECTION*)a)->distance < ((PROJECTION*)b)->distance
          ? 1 : 0 );
}


/* ----------------------------------------------------------------- */


static void SortDoubles(int size, PROJECTION projection[])
{
  if( size <= 0 )
    {
    printfe("ERROR: Sorting %i doubles\n", size);
    exit( -1 );
    }

  QuickSort(projection, size, sizeof(PROJECTION), NULL, cmpdouble);
}


/* ----------------------------------------------------------------- */


static void SearchOptimalPivotEquitz(TRAININGSET*  TS,
                                     PARTITIONING* P,
                                     int           Pindex,
                                     PROJECTION    projection[],
                                     int*          pivotindex)
{
  int      i, j, k;

  llong*   C1  = AllocateCounterVector(VectorSize(TS));
  llong*   C2  = AllocateCounterVector(VectorSize(TS));
  int      fr1 = 0;
  int      fr2 = CCFreq(P, Pindex);
  double   D1, D2, Dtot, Dmin = DBL_MAX;
  BOOKNODE cent = CreateEmptyNode(VectorSize(TS));
  double   error1, error2;
  int      X, Xfr;

  *pivotindex = -1;
  for( k = 0; k < VectorSize(TS); k++ )
    {
    C1[k] = 0LL;
    C2[k] = CCScalar(P, Pindex, k);
    }

  PartitionCentroid(P, Pindex, &cent);
  D1 = 0LL;
  D2 = PartitionError(TS, P, Pindex, &cent, MAXLLONG);

  /* Check all n-1 split points. */
  for( i = 0; i < UniqueVectors(P, Pindex) - 1; i++ )
    {
    j = projection[i].vindex;
    Xfr = VectorFreq(TS, j);

    /* 'To partition' case */
    error1 = 0;
    if( fr1 > 0 )
      {
      for( k = 0; k < VectorSize(TS); k++ )
        {
        X = VectorScalar(TS, j, k);
        error1 += sqr( ((double)C1[k] / (double)fr1) - (double)X);
        }
      }
    D1 += ((double)(fr1 * Xfr) / (double)(fr1 + Xfr)) * error1;

    /* New centroids */
    for( k = 0; k < VectorSize(TS); k++ )
       {
       X = VectorScalar(TS, j, k);
       C1[k] += X * Xfr;
       C2[k] -= X * Xfr;
       }
    fr1 += VectorFreq(TS, j);
    fr2 -= VectorFreq(TS, j);

    /* 'From partition' case */
    error2 = 0;
    for( k = 0; k < VectorSize(TS); k++ )
      {
      X = VectorScalar(TS, j, k);
      error2 += sqr( ((double)C2[k] / (double)fr2) - (double)X);
      }
    D2 -= ((double)(fr2 * Xfr) / (double)(fr2 + Xfr)) * error2;

    /* Combine errors and compare to the current minimum. */
    Dtot = D1 + D2;
    if( Dtot < Dmin )
      {
      Dmin = Dtot;
      *pivotindex = i;
      }
    }

  FreeNode(cent);
  deallocate(C1);
  deallocate(C2);
}


/* ----------------------------------------------------------------- */


static void SearchOptimalPivotWu(TRAININGSET*  TS,
                                 PARTITIONING* P,
                                 int           Pindex,
                                 PROJECTION    projection[],
                                 int*          pivotindex)
{
  int      i, j, k;

  llong*   Sum1  = AllocateCounterVector(VectorSize(TS));
  llong*   Sum2  = AllocateCounterVector(VectorSize(TS));
  int      fr1 = 0;
  int      fr2 = CCFreq(P, Pindex);
  double   result, maxresult = 0;
  double    S1, S2;
  int      Xchange, Xfr;

  *pivotindex = -1;
  for( k = 0; k < VectorSize(TS); k++ )
    {
    Sum1[k] = 0LL;
    Sum2[k] = CCScalar(P, Pindex, k);
    }

  /* Check all n-1 split points. */
  for( i = 0; i < UniqueVectors(P, Pindex) - 1; i++ )
    {
    j = projection[i].vindex;
    Xfr = VectorFreq(TS, j);

    /* Move from cluster2 to cluster1. */
    for( k = 0; k < VectorSize(TS); k++ )
      {
      Xchange = VectorScalar(TS, j, k) * Xfr;
      Sum1[k] += Xchange;
      Sum2[k] -= Xchange;
      }
    fr1 += Xfr;
    fr2 -= Xfr;

    S1 = 0.0;
    S2 = 0.0;
    for( k = 0; k < VectorSize(TS); k++ )
      {
      S1 += sqr((double)Sum1[k]);
      S2 += sqr((double)Sum2[k]);
      }
    result = S1 / fr1 + S2 / fr2;

    /* Compare result to the current maximum. */
    if( result > maxresult )
      {
      maxresult = result;
      *pivotindex = i;
      }

    if( Value(ShowProgress) >= 4 )
      {
      printf("OptWu res=%11.4f max=%11.4f pivotind=%4i\n",
              result, maxresult, *pivotindex);
      }
    }

  deallocate(Sum1);
  deallocate(Sum2);
}


/* ----------------------------------------------------------------- */


static void SearchOptimalPivot(TRAININGSET*  TS,
                               PARTITIONING* P,
                               int           Pindex,
                               PROJECTION    projection[],
                               int*          pivotindex)
{
  int      i, j, k;

  llong*   C1  = AllocateCounterVector(VectorSize(TS));
  llong*   C2  = AllocateCounterVector(VectorSize(TS));
  int      fr1 = 0;
  int      fr2 = CCFreq(P, Pindex);
  double   D1, D2, Dtot, Dmin = DBL_MAX; 
  BOOKNODE cent = CreateEmptyNode(VectorSize(TS));
  int      X, Xfr;
  int      l, m;

  *pivotindex = -1;
  for( k = 0; k < VectorSize(TS); k++ )
    {
    C1[k] = 0LL;
    C2[k] = CCScalar(P, Pindex, k);
    }

  /* Check all n-1 split points. */
  for( i = 0; i < UniqueVectors(P, Pindex) - 1; i++ )
    {
    j = projection[i].vindex;
    Xfr = VectorFreq(TS, j);

    /* New centroids */
    for( k = 0; k < VectorSize(TS); k++ )
       {
       X = VectorScalar(TS, j, k);
       C1[k] += X * Xfr;
       C2[k] -= X * Xfr;
       }
    fr1 += Xfr;
    fr2 -= Xfr;

    D1 = 0;
    for( k = 0; k < VectorSize(TS); k++ )
      {
      cent.vector[k] = (int)round( C1[k] / (double)(fr1));
      }
    for( l = 0; l <= i; l++ )
      {
      m = projection[l].vindex;
      D1 += MSEDifference(&cent, &Node(TS, m), VectorSize(TS), MAXLLONG)
            * VectorFreq(TS, m);
      }

    D2 = 0;
    for( k = 0; k < VectorSize(TS); k++ )
       {
       cent.vector[k] = (int)round( C2[k] / (double)(fr2));
       }
    for( l = i+1; l < UniqueVectors(P, Pindex); l++ )
      {
      m = projection[l].vindex;
      D2 += MSEDifference(&cent, &Node(TS, m), VectorSize(TS), MAXLLONG)
            * VectorFreq(TS, m);
      }

    /* Combine errors and compare to the current minimum. */
    Dtot = D1 + D2;
    if( Dtot < Dmin )
      {
      Dmin = Dtot;
      *pivotindex = i;
      }
    }

  FreeNode(cent);
  deallocate(C1);
  deallocate(C2);
}


/* ----------------------------------------------------------------- */


static int DimensionWithHighestVariance(TRAININGSET*  TS,
                                        PARTITIONING* P,
                                        int           Pindex,
                                        double        avg[])
{
  int j, k;
  double tmp;
  int    hdim;
  double* varsum = (double*)allocate(VectorSize(TS) * sizeof(double));

  for( k = 0; k < VectorSize(TS); k++ )
    {
    varsum[k] = 0.0;
    }

  /* For each dimension calculate the sum of square differences
     from the average. (~variance) */
  for( j = FirstVector(P,Pindex); ! EndOfPartition(j); j = NextVector(P,j) )
    {
    for( k = 0; k < VectorSize(TS); k++ )
      {
      tmp = VectorScalar(TS, j, k) - avg[k];
      varsum[k] += sqr(tmp) * VectorFreq(TS, j);
      }
    }

  /* Select dimension with highest variance. */
  hdim = 0;
  for( k = 1; k < VectorSize(TS); k++ )
    {
    if( varsum[k] > varsum[hdim] )
      {
      hdim = k;
      }
    }

  deallocate(varsum);

  return( hdim );
}


/* ----------------------------------------------------------------- */


static void RadiusWeightedMean(TRAININGSET*  TS,
                               PARTITIONING* P,
                               int           Pindex,
                               double        avg[],
                               double        rwm[])
{
  int j, k;
  double r, rsum = 0.0;

  for( k = 0; k < VectorSize(TS); k++ )
    {
    rwm[k] = 0.0;
    }
  for( j = FirstVector(P,Pindex); ! EndOfPartition(j); j = NextVector(P,j) )
    {
    r = 0.0;
    for( k = 0; k < VectorSize(TS); k++ )
      {
      r += sqr(VectorScalar(TS, j, k) - avg[k]);
      }
    r = sqrt(r);
    rsum += r * VectorFreq(TS, j);
    for( k = 0; k < VectorSize(TS); k++ )
      {
      rwm[k] += r * VectorFreq(TS, j) * VectorScalar(TS, j, k);
      }
    }
  for( k = 0; k < VectorSize(TS); k++ )
    {
    rwm[k] /= rsum;

    if( Value(ShowProgress) >= 4 )
      {
      printf("k=%2i avg=%9.4f rwm=%9.4f rsum=%9.4f UniVectors=%5i\n",
              k, avg[k], rwm[k], rsum, UniqueVectors(P, Pindex));
      }
    }
}


/* ================================================================= */


static double* MxV(double* m[], double v[], double v2[], int dim)
{
  int r, c;

  for( r = 0; r < dim; r++ )
    {
    v2[r] = 0.0;
    for( c = 0; c < dim; c++ )
      {
      v2[r] += m[r][c] * v[c];
      }
    }
  return( v2 );
}


/* ----------------------------------------------------------------- */


static double Vnorm(double v[], int dim, int p)
/* Calculates the Lp norm of vector v.
   When p == 0, infinite norm is calculated. */
{
  int    k;

  switch( p )
    {
    case 1:
      {
      double sum = 0.0;
      for( k = 0; k < dim; k++ )
        {
        sum += absd(v[k]);
        }
      return( sum );
      }
    case 2:
      {
      double sum = 0.0;
      for( k = 0; k < dim; k++ )
        {
        sum += v[k] * v[k];
        }
      return( sqrt(sum) );
      }
    case 0:
      {
      double maxelem = absd(v[0]);
      for( k = 1; k < dim; k++ )
        {
        if( absd(v[k]) > maxelem )
          {
          maxelem = absd(v[k]);
          }
        }
      return( maxelem );
      }
    default:
      {
      printf("ERROR: In Vnorm p==%i is not supported.\n", p);
      exit( -1 );
      }
    }
}


/* ----------------------------------------------------------------- */

static double* VxS(double v[], int dim, double s)
{
  int k;

  for( k = 0; k < dim; k++ )
    {
    v[k] *= s;
    }
  return( v );
}


/* ----------------------------------------------------------------- */


static double* VsubV(double v1[], double v2[], double diff[], int dim)
{
  int k;

  for( k = 0; k < dim; k++ )
    {
    diff[k] = v1[k] - v2[k];
    }
  return( diff );
}


/* ----------------------------------------------------------------- */


static double* Vcopy(double source[], double dest[], int dim)
{
  int k;

  for( k = 0; k < dim; k++ )
    {
    dest[k] = source[k];
    }
  return( dest );
}


/* ----------------------------------------------------------------- */


static void CalculateCovarianceMatrix(TRAININGSET*  TS,
                                      PARTITIONING* P,
                                      int           Pindex,
                                      double*       covariance[],
                                      double        avg[])
{
  int     j, k;
  int     r, c;
  double* diff = (double*)allocate(VectorSize(TS) * sizeof(double));

  /* Initialize */
  for( k = 0; k < VectorSize(TS); k++ )
    {
    avg[k] = (double)CCScalar(P, Pindex, k) / (double)CCFreq(P, Pindex);

    for( c = 0; c < VectorSize(TS); c++ )
      {
      covariance[k][c] = 0.0;
      }
    }

  /* Calculate covariance matrix. */
  for( j = FirstVector(P,Pindex); ! EndOfPartition(j); j = NextVector(P,j) )
    {
    for( k = 0; k < VectorSize(TS); k++ )
      {
      diff[k] = (double)VectorScalar(TS, j, k) - avg[k];
      }

    for( r = 0; r < VectorSize(TS); r++ )
      {
      for( c = r; c < VectorSize(TS); c++ )
        {
        covariance[r][c] += diff[r] * diff[c] * VectorFreq(TS, j);
        }
      }
    }

  /* Copy upper triangular matrix to lower triangular matrix. */
  /* (Loop through the target area.) */
  for( r = 1; r < VectorSize(TS); r++ )
    {
    for( c = 0; c < r; c++ )
      {
      covariance[r][c] = covariance[c][r];
      }
    }

  /* Normalize */
  for( r = 0; r < VectorSize(TS); r++ )
    {
    for( c = 0; c < VectorSize(TS); c++ )
      {
      covariance[r][c] /= CCFreq(P, Pindex);
      }
    }

  if( Value(ShowProgress) >= 4 )
    {
    printf("Avg %i:\n", Pindex);
    for( k = 0; k < VectorSize(TS); k++ )
      {
      printf("%9.4f ", avg[k]);
      }
    printf("\nCovariance %i:\n", Pindex);
    for( r = 0; r < VectorSize(TS); r++ )
      {
      for( c = 0; c < VectorSize(TS); c++ )
        {
        printf("%9.4f ", covariance[r][c]);
        }
      printf("\n");
      }
    }
  deallocate(diff);
}


/* ----------------------------------------------------------------- */


static void InitEigenVector(double v[], int dim, int norm)
{
  int k;
  double vnorm;

  for( k = 0; k < dim; k++ )
    {
    v[k] = drand();
    }
  vnorm = Vnorm(v, dim, norm);
  VxS(v, dim, 1/vnorm);
}


/* ----------------------------------------------------------------- */


static void PowerMethod(double*  covariance[],
                        int      dim,
                        double   eigenvector[])
{
  double  errorlimit = 1e-10;
  double  error      = DBL_MAX;
  int     k;
  double* diff = (double*)allocate(dim * sizeof(double));
  double* y = (double*)allocate(dim * sizeof(double));
  double  ynorm;
  int     i = 0;
  double  initeigenelem = 1.0 / sqrt((double)dim);

  /* Initialize. */
  for( k = 0; k < dim; k++ )
    {
    eigenvector[k] = initeigenelem;
    }

  while( error > errorlimit )
    {
    MxV(covariance, eigenvector, y, dim);
/*     eigenvalue = VtxV(eigenvector, y, dim); */
    ynorm = Vnorm(y, dim, 2);
    if( ynorm == 0.0 )
      {
      InitEigenVector(eigenvector, dim, 2);
      error = DBL_MAX;
      }
    else
      {
      VxS(y, dim, 1/ynorm);
      VsubV(eigenvector, y, diff, dim);
      error = Vnorm(diff, dim, 2);
      Vcopy(y, eigenvector, dim);
      }

    if( Value(ShowProgress) >= 4 )
      {
      printf("PM %2i: error=%11.8f ynorm=%13.8f, ", i, error, ynorm);
      for( k = 0; k < dim; k++ )
        {
        printf("%11.8f ", eigenvector[k]);
        }
      printf("\n");
      }
    i++;
    }

  deallocate(diff);
  deallocate(y);
}


/* ----------------------------------------------------------------- */


static void CalculatePCA(TRAININGSET*  TS,
                         PARTITIONING* P,
                         int           Pindex,
                         double        pca[],
                         double        avg[])
{
  int      k;
  double** covariance = (double**)allocate(VectorSize(TS) * sizeof(double*));

/*
  if( UniqueVectors(P, Pindex) == 2 )
    {
    printf("CalcPCA\n");
    PrintPartition(TS, P, Pindex);
    }
*/

  for( k = 0; k < VectorSize(TS); k++ )
    {
    covariance[k] = (double*)allocate(VectorSize(TS) * sizeof(double));
    }
  CalculateCovarianceMatrix(TS, P, Pindex, covariance, avg);
  PowerMethod(covariance, VectorSize(TS), pca);

  if( Value(ShowProgress) >= 3 )
    {
    double infitenorm = Vnorm(pca, VectorSize(TS), 0);

    VxS(pca, VectorSize(TS), 1/infitenorm);
    printf("PCA:\n");
    for( k = 0; k < VectorSize(TS); k++ )
      {
      printf("k=%2i: pca=%11.8f avg=%9.4f\n", k, pca[k], avg[k]);
      }
    }

  for( k = 0; k < VectorSize(TS); k++ )
    {
    deallocate(covariance[k]);
    }
  deallocate(covariance);
}


/* ----------------------------------------------------------------- */


static void CalculateRegression(TRAININGSET*  TS,
                                PARTITIONING* P,
                                int           Pindex,
                                double        regA[],
                                double        regB[],
                                double        avg[])
{
  int     j, k;
  double* product = (double*)allocate(VectorSize(TS) * sizeof(double));
  double  tmp1, tmp2;
  int     basek;


  if( UniqueVectors(P, Pindex) <= 1 )
    {
    printfe("Regression of small partition. Freq=%i unique=%i\n",
            CCFreq(P, Pindex), UniqueVectors(P, Pindex));
    exit( -1 );
    }

  for( k = 0; k < VectorSize(TS); k++ )
    {
    product[k] = 0.0;
    avg[k] = (double)CCScalar(P, Pindex, k) / (double)CCFreq(P, Pindex);
    }

  basek = DimensionWithHighestVariance(TS, P, Pindex, avg);

  for( j = FirstVector(P,Pindex); ! EndOfPartition(j); j = NextVector(P,j) )
    {
    for( k = 0; k < VectorSize(TS); k++ )
      {
      product[k] += (double)VectorScalar(TS, j, basek)
                    * (double)VectorScalar(TS, j, k)
                    * VectorFreq(TS, j);
      }
    }

  regA[basek] = 1.0;
  regB[basek] = 0.0;
  tmp1 = product[basek] - ( sqr((double)CCScalar(P, Pindex, basek)) /
                            (double)CCFreq(P, Pindex) );
  if( tmp1 == 0.0 )
    {
    printfe("Divisor of regression is zero.\n");
    printf("Pindex=%4i tmp1=%9.4f Psum=%f Pfreq=%i basek=%i\n",
           Pindex, tmp1, (double)CCScalar(P, Pindex, basek), CCFreq(P, Pindex), basek);
    for( k = 0; k < VectorSize(TS); k++ )
      {
      printf("d=%2i avg=%9.4f product=%13.4f ", k, avg[k], product[k]);
      printf("A=%9.4f B=%9.4f\n", regA[k], regB[k]);
      }
    PrintPartition(TS, P, Pindex);
    exit( -1 );
    }

  for( k = 0; k < VectorSize(TS); k++ )
    {
    tmp2 = (double)CCScalar(P, Pindex, basek) * (double)CCScalar(P, Pindex, k)
           / (double)CCFreq(P, Pindex);
    regA[k] = (product[k] - tmp2) / tmp1;
    regB[k] = avg[k] - regA[k] * avg[basek];
    }

  if( Value(ShowProgress) >= 3 )
    {
    printf("basek=%2i\n", basek);
    for( k = 0; k < VectorSize(TS); k++ )
      {
      printf("d=%2i avg=%9.4f product=%13.4f ", k, avg[k], product[k]);
      printf("RegA=%9.4f RegB=%9.4f\n", regA[k], regB[k]);
      }
    }

  deallocate(product);
}


/* ----------------------------------------------------------------- */


static void CalculateFurthestLine(TRAININGSET*  TS,
                                  CODEBOOK*     CB,
                                  PARTITIONING* P,
                                  int           Pindex,
                                  double        fline[],
                                  double        avg[])
{
  int    k;
  int    V1, V2;
  double dist;
  double numerator, denominator;

  TwoFurthestVectorsInPartition(TS, CB, P, Pindex, &V1, &V2);

  for( k = 0; k < VectorSize(TS); k++ )
    {
    fline[k] = VectorScalar(TS, V2, k) - VectorScalar(TS, V1, k);
    }

  numerator = 0.0;
  denominator = 0.0;
  for( k = 0; k < VectorSize(TS); k++ )
    {
    numerator += fline[k] * (VectorScalar(CB, Pindex, k) - VectorScalar(TS, V1, k));
    denominator += sqr((double)VectorScalar(TS, V1, k));
    }
  if( denominator == 0.0 )
    {
    printfe("ERROR: denominator=%f\n", denominator);
    exit( -1 );
    }

  dist = numerator / denominator;

  for( k = 0; k < VectorSize(TS); k++ )
    {
    avg[k] = VectorScalar(TS, V1, k) + dist * fline[k];
    }

  if( Value(ShowProgress) >= 3 )
    {
    printf("Furthest (%4i):\n", Pindex);
    for( k = 0; k < VectorSize(TS); k++ )
      {
      printf("d=%2i avg=%9.4f fline=%13.4f\n", k, avg[k], fline[k]);
      }
    }
}


/* ----------------------------------------------------------------- */


static void CalculateProjectionDistances(TRAININGSET*  TS,
                                         PARTITIONING* P,
                                         int           Pold,
                                         double        pca[],
                                         double        avg[],
                                         PROJECTION    projection[],
                                         double*       mindist,
                                         double*       maxdist)
{
  int     j, k;
  int     dindex = 0;
  double  sumpca2;
  double  sumpcapivot;
  double  sumpcapoint;
  double  dist2pivot;
  double* pivot = (double*)allocate(VectorSize(TS) * sizeof(double));

  if( Value(HyperplanePivot) == HPPWeightedCentroid )
    {
    RadiusWeightedMean(TS, P, Pold, avg, pivot);
    }
  else
    {
    for( k = 0; k < VectorSize(TS); k++ )
      {
      pivot[k] = avg[k];
      }
    }

  sumpca2     = 0.0;
  sumpcapivot = 0.0;
  for( k = 0; k < VectorSize(TS); k++ )
    {
    sumpca2     += sqr(pca[k]);
    sumpcapivot += pca[k] * pivot[k];
    }
  if( sumpca2 == 0.0 )
    {
    printfe("ERROR: sumpca2=%f\n", sumpca2);
    exit( -1 );
    }

  *mindist = 0.0;
  *maxdist = 0.0;
  for( j = FirstVector(P, Pold); ! EndOfPartition(j); j = NextVector(P, j) )
    {
    sumpcapoint = 0.0;
    for( k = 0; k < VectorSize(TS); k++ )
      {
      sumpcapoint += pca[k] * VectorScalar(TS, j , k);
      }

    dist2pivot = (sumpcapoint - sumpcapivot) / sumpca2;

    projection[dindex].distance = dist2pivot;
    projection[dindex].vindex = j;
    dindex++;

    if( dist2pivot < *mindist )
      {
      *mindist = dist2pivot;
      }
    else
      {
      if( dist2pivot > *maxdist )
        {
        *maxdist = dist2pivot;
        }
      }

    if( Value(ShowProgress) >= 4 )
      {
      printf("j=%4i pcapoint=%9.4f pcapivot=%9.4f pca2=%9.4f dist[%4i]=%9.4f mind=%9.4f maxd=%9.4f\n",
             j, sumpcapoint, sumpcapivot, sumpca2, dindex-1, projection[dindex-1].distance, *mindist, *maxdist);
      }
    }
  deallocate(pivot);
}


/* ----------------------------------------------------------------- */


static void SplitByRegressionPivot(TRAININGSET*    TS,
                                   CODEBOOK*       CB,
                                   PARTITIONING*   P,
                                   CLUSTERINFO     CI[],
                                   PREPARTITIONING PP[],
                                   int             Pindex,
                                   PROJECTION      projection[],
                                   double          newzero)
/* Determinates the values of CI[Pold].Centroid1 & 2 and PP. */
{
  int    j, k;
  llong* C1 = AllocateCounterVector(VectorSize(TS));
  llong* C2 = AllocateCounterVector(VectorSize(TS));
  int    freq1 = 0, freq2 = 0;
  int    p;
  int    origPsize = UniqueVectors(P, Pindex);

  for( k = 0; k < VectorSize(TS); k++ )
    {
    C1[k] = 0LL;
    C2[k] = 0LL;
    }

  for( p = 0; p < origPsize; p++ )
    {
    j = projection[p].vindex;
    /* Vectors projected to newzero belong to the 'new' cluster (C2). */
    if( projection[p].distance > newzero )
      {
      PP[j] = 1;
      for( k = 0; k < VectorSize(TS); k++ )
        {
        C1[k] += VectorScalar(TS, j, k) * VectorFreq(TS, j);
        }
      freq1 += VectorFreq(TS, j);
      }
    else
      {
      PP[j] = 2;
      for( k = 0; k < VectorSize(TS); k++ )
        {
        C2[k] += VectorScalar(TS, j, k) * VectorFreq(TS, j);
        }
      freq2 += VectorFreq(TS, j);
      }
    if( Value(ShowProgress) >= 4 )
      {
      printf("p=%4i distance[]=%9.4f vindex[]=%4i cluster=%i\n",
             p, projection[p].distance, projection[p].vindex, PP[j]);
      }
    }

  if( freq1 == 0LL || freq2 == 0LL )
    {
    printf("Error: Pindex=%4i freq1=%5i freq2=%5i newzero=%f\n",
           Pindex, freq1, freq2, newzero);
    for( k = 0; k < VectorSize(TS); k++ )
      {
      printf("k=%2i C1=%10.1f C2=%10.1f\n", k, (double)C1[k], (double)C2[k]);
      }
    for( p = 0; p < origPsize; p++ )
      {
      printf("p=%4i distance[]=%10.7f vindex[]=%4i cluster=%i\n",
             p, projection[p].distance, projection[p].vindex, PP[projection[p].vindex]);
      }
    PrintPartition(TS, P, Pindex);
    exit( -1 );
    }

  for( k = 0; k < VectorSize(TS); k++ )
    {
    CI[Pindex].Centroid1[k] = (VECTORELEMENT)round((double)C1[k]/(double)freq1);
    CI[Pindex].Centroid2[k] = (VECTORELEMENT)round((double)C2[k]/(double)freq2);
    }

  deallocate(C1);
  deallocate(C2);
}


/* ----------------------------------------------------------------- */


static void SplitByRegressionIndex(TRAININGSET*    TS,
                                   CODEBOOK*       CB,
                                   PARTITIONING*   P,
                                   CLUSTERINFO     CI[],
                                   PREPARTITIONING PP[],
                                   int             Pindex,
                                   PROJECTION      projection[],
                                   int             pivotindex)
/* Determinates the values of CI[Pold].Centroid1 & 2 and PP. */
{
  int    j, k;
  llong* C1 = AllocateCounterVector(VectorSize(TS));
  llong* C2 = AllocateCounterVector(VectorSize(TS));
  int    freq1 = 0, freq2 = 0;
  int    p;
  int    origPsize = UniqueVectors(P, Pindex);

  assert( pivotindex < origPsize - 1 );

  for( k = 0; k < VectorSize(TS); k++ )
    {
    C1[k] = 0LL;
    C2[k] = 0LL;
    }

  for( p = 0; p <= pivotindex; p++ )
    {
    j = projection[p].vindex;

    PP[j] = 1;
    for( k = 0; k < VectorSize(TS); k++ )
      {
      C1[k] += VectorScalar(TS, j, k) * VectorFreq(TS, j);
      }
    freq1 += VectorFreq(TS, j);

    if( Value(ShowProgress) >= 4 )
      {
      printf("p=%4i distance[]=%9.4f vindex[]=%4i cluster=%i\n",
             p, projection[p].distance, projection[p].vindex, PP[j]);
      }
    }

  for( p = pivotindex+1; p < origPsize; p++ )
    {
    j = projection[p].vindex;

    PP[j] = 2;
    for( k = 0; k < VectorSize(TS); k++ )
      {
      C2[k] += VectorScalar(TS, j, k) * VectorFreq(TS, j);
      }
    freq2 += VectorFreq(TS, j);

    if( Value(ShowProgress) >= 4 )
      {
      printf("p=%4i distance[]=%9.4f vindex[]=%4i cluster=%i\n",
             p, projection[p].distance, projection[p].vindex, PP[j]);
      }
    }

  if( freq1 == 0LL || freq2 == 0LL )
    {
    printf("Error: Pindex=%4i freq1=%5i freq2=%5i pivotindex=%4i\n",
           Pindex, freq1, freq2, pivotindex);
    for( k = 0; k < VectorSize(TS); k++ )
      {
      printf("k=%2i C1=%10.1f C2=%10.1f\n", k, (double)C1[k], (double)C2[k]);
      }
    for( p = 0; p < origPsize; p++ )
      {
      printf("p=%4i distance[]=%10.7f vindex[]=%4i cluster=%i\n",
             p, projection[p].distance, projection[p].vindex, PP[projection[p].vindex]);
      }
    PrintPartition(TS, P, Pindex);
    exit( -1 );
    }

  for( k = 0; k < VectorSize(TS); k++ )
    {
    CI[Pindex].Centroid1[k] = (VECTORELEMENT)round((double)C1[k]/(double)freq1);
    CI[Pindex].Centroid2[k] = (VECTORELEMENT)round((double)C2[k]/(double)freq2);
    }

  deallocate(C1);
  deallocate(C2);
}


/* ----------------------------------------------------------------- */


static void SplitByHyperPlane(TRAININGSET*    TS,
                              CODEBOOK*       CB,
                              PARTITIONING*   P,
                              CLUSTERINFO     CI[],
                              PREPARTITIONING PP[],
                              int             Pold)
{
  double*  pca  = (double*)allocate(VectorSize(TS) * sizeof(double));
  double*  regB = (double*)allocate(VectorSize(TS) * sizeof(double));
  double*  avg  = (double*)allocate(VectorSize(TS) * sizeof(double));
  PROJECTION* projection = (PROJECTION*)allocate(UniqueVectors(P, Pold) * sizeof(PROJECTION));
  double   mindist, maxdist;
  double   pivot;
  int      pivotindex;

  switch( Value(HyperplaneLine) )
    {
    case HPLPCA:
      CalculatePCA(TS, P, Pold, pca, avg);
      break;
    case HPLRegression:
      CalculateRegression(TS, P, Pold, pca, regB, avg);
      break;
    case HPLFurthest:
      CalculateFurthestLine(TS, CB, P, Pold, pca, avg);
      break;
    }

  CalculateProjectionDistances(TS, P, Pold, pca, avg,
                               projection, &mindist, &maxdist);

  switch( Value(HyperplanePivot) )
    {
    case HPPCentroid:
      pivot = 0.0;
      SplitByRegressionPivot(TS, CB, P, CI, PP, Pold, projection, pivot);
      break;
    case HPPOptimal:
      SortDoubles(UniqueVectors(P, Pold), projection);
      SearchOptimalPivot(TS, P, Pold, projection, &pivotindex);
      SplitByRegressionIndex(TS, CB, P, CI, PP, Pold, projection, pivotindex);
      break;
    case HPPOptimalEquitz:
      SortDoubles(UniqueVectors(P, Pold), projection);
      SearchOptimalPivotEquitz(TS, P, Pold, projection, &pivotindex);
      SplitByRegressionIndex(TS, CB, P, CI, PP, Pold, projection, pivotindex);
      break;
    case HPPWeightedCentroid:
      /* Calculated in 'CalculateProjectionDistances' */
      pivot = 0.0;
      SplitByRegressionPivot(TS, CB, P, CI, PP, Pold, projection, pivot);
      break;
    case HPPGoeddelBass:
      pivot = (mindist + maxdist) / 2.0;
      SplitByRegressionPivot(TS, CB, P, CI, PP, Pold, projection, pivot);
      break;
    case HPPMomentPreserving:
      printf("HPPMomentPreserving is not implemented.\n");
      exit( -1 );
    case HPPLloydScalarQuantization:
      LloydScalarQuantization(projection, UniqueVectors(P, Pold), &pivot);
      SplitByRegressionPivot(TS, CB, P, CI, PP, Pold, projection, pivot);
      break;
    case HPPOptimalWu:
      SortDoubles(UniqueVectors(P, Pold), projection);
      SearchOptimalPivotWu(TS, P, Pold, projection, &pivotindex);
      SplitByRegressionIndex(TS, CB, P, CI, PP, Pold, projection, pivotindex);
      break;
    }

  deallocate(pca);
  deallocate(regB);
  deallocate(avg);
  deallocate(projection);
}


/* ================================================================= */


static int SelectTwoCodevectors(TRAININGSET*  TS,
                                CODEBOOK*     CB,
                                PARTITIONING* P,
                                CLUSTERINFO   CI[],
                                int           Pindex)
/* This just selects the values for CI[Pindex].Centroid1 and
   CI[Pindex].Centroid2.
   Returns non-zero if success. */
{
  llong*     SD;
  int        tmpV1 = -1, tmpV2 = -1;
  VECTORTYPE tmp1 = CI[Pindex].Centroid1;
  VECTORTYPE tmp2 = CI[Pindex].Centroid2;

  /* Ensure that there is something to split. */
  if( UniqueVectors(P, Pindex) < 2 )
  {
    return 0;
  } 

  switch( Value(NewVectors) )
    {
    case CurrentAndSigma:
      /* Current... */
      CopyVector(Vector(CB, Pindex), tmp1, VectorSize(CB));

      /* ...and centroid + standard deviation of partition. */
      SD = AllocateCounterVector(VectorSize(TS));
      PartitionSD(TS, P, Pindex, SD);
      CopyVector(Vector(CB, Pindex), tmp2, VectorSize(CB));
      AddCounterToVector(SD, tmp2, VectorSize(TS), TS->MaxValue);
      deallocate(SD);
      break;
    case CurrentPlusMinusSigma:
      /* Initialize new vectors with centroid. */
      CopyVector(Vector(CB, Pindex), tmp1, VectorSize(CB));
      CopyVector(Vector(CB, Pindex), tmp2, VectorSize(CB));

      /* Calculate the SD vector of the partition. */
      SD = AllocateCounterVector(VectorSize(TS));
      PartitionSD(TS, P, Pindex, SD);

      /* Modify codevectors with SD vector. */
      AddCounterToVector(SD, tmp1, VectorSize(TS), TS->MaxValue);
      SubtractCounterFromVector(SD, tmp2, VectorSize(TS), TS->MinValue);
      deallocate(SD);
      break;
    case CurrentAndFurthest:
      /* Current... */
      CopyVector(Vector(CB, Pindex), tmp1, VectorSize(CB));

      /* ...and furthest. */
      FurthestVectorInPartition(TS, P, Pindex, &Node(CB, Pindex), &tmpV2);
      CopyVector(Vector(TS, tmpV2), tmp2, VectorSize(CB));
      break;
    case nv4:
      printf("Not used.\n");
      exit( -1 );
    case CurrentAndRandom:
      /* Find random vector from partition. */
      tmpV1 = RandomVectorFromPartition(TS, P, Pindex, CCFreq(P, Pindex), ENDPARTITION);
      tmpV2 = RandomVectorFromPartition(TS, P, Pindex, CCFreq(P, Pindex), tmpV1);

      CopyVector(Vector(CB, Pindex), tmp1, VectorSize(CB));
      if( EqualVectors(Vector(TS, tmpV1), Vector(CB, Pindex), VectorSize(TS)) )
        {
        CopyVector(Vector(TS, tmpV2), tmp2, VectorSize(CB));
        }
      else
        {
        CopyVector(Vector(TS, tmpV1), tmp2, VectorSize(CB));
        }
      break;
    case TwoFurthest:
      TwoFurthestVectorsInPartition(TS, CB, P, Pindex, &tmpV1, &tmpV2);
      CopyVector(Vector(TS, tmpV1), tmp1, VectorSize(CB));
      CopyVector(Vector(TS, tmpV2), tmp2, VectorSize(CB));
      break;
    case TwoRandom:
      /* We need two vectors. But if there are only two vectors
         selection is simple. */
      if( UniqueVectors(P, Pindex) == 2 )
        {
        tmpV1 = FirstVector(P, Pindex);
        tmpV2 = NextVector(P, tmpV1);
        }
      else
        {
        /* Two random vectors... */
        tmpV1 = RandomVectorFromPartition(TS, P, Pindex, CCFreq(P, Pindex), ENDPARTITION);
        tmpV2 = RandomVectorFromPartition(TS, P, Pindex, CCFreq(P, Pindex), tmpV1);
        }
      CopyVector(Vector(TS, tmpV1), tmp1, VectorSize(CB));
      CopyVector(Vector(TS, tmpV2), tmp2, VectorSize(CB));
      break;
    case MeansOfCurrentAndFurthest:
    case FasterVariant:
      TwoFurthestVectorsInPartition(TS, CB, P, Pindex, &tmpV1, &tmpV2);
      MeanOfTwoVectors(Vector(CB, Pindex), Vector(TS, tmpV1), tmp1, VectorSize(CB));
      MeanOfTwoVectors(Vector(CB, Pindex), Vector(TS, tmpV2), tmp2, VectorSize(CB));

      /* Are means equal? */
      if( EqualVectors(tmp1, tmp2, VectorSize(TS)) )
        {
        /* Yes, use instead the furthest vectors in partition. */
        CopyVector(Vector(TS, tmpV1), tmp1, VectorSize(CB));
        CopyVector(Vector(TS, tmpV2), tmp2, VectorSize(CB));
        }
      break;
    default:
      break;
    }
  if( Value(ShowProgress) >= 4 )
    {
    printf("Centroid: ");PrintVector(Vector(CB, Pindex), VectorSize(TS), 1);
    printf("Selected1:");PrintVector(tmp1, VectorSize(TS), 1);
    printf("Selected2:");PrintVector(tmp2, VectorSize(TS), 1);
    }

  if( EqualVectors(tmp1, tmp2, VectorSize(TS)) )
    {
    printfe("ERROR: SelectTwoCodevectors Pindex=%i tmpV1=%i tmpV2=%i pfreq=%i Size=%i\n",
           Pindex, tmpV1, tmpV2, CCFreq(P, Pindex), BookSize(CB));
    exit( -1 );
    }

  return( 1 );
}

/* ----------------------------------------------------------------- */


static llong TentativeError(TRAININGSET*    TS,
                            CODEBOOK*       CB,
                            PARTITIONING*   P,
                            CLUSTERINFO     CI[],
                            PREPARTITIONING PP[],
                            int             Pindex)
{
  int   j;
  llong currerror, error;
  llong totalerror = 0LL;

  if( Value(ShowProgress) >= 4 )
    {
    printf("Tent.C: ");PrintVector(Vector(CB, Pindex), VectorSize(TS), 1);
    printf("Tent.C1:");PrintVector(CI[Pindex].Centroid1, VectorSize(TS), 1);
    printf("Tent.C2:");PrintVector(CI[Pindex].Centroid2, VectorSize(TS), 1);
    }

  for( j = FirstVector(P,Pindex); ! EndOfPartition(j); j = NextVector(P,j) )
    {
    currerror = VectorDist(Vector(TS, j), Vector(CB, Pindex), VectorSize(TS));
    error = VectorDist(Vector(TS, j),
                       PP[j] == 1 ? CI[Pindex].Centroid1 : CI[Pindex].Centroid2,
                       VectorSize(TS));
    totalerror += (currerror - error) * VectorFreq(TS, j);

    if( Value(ShowProgress) >= 4 )
      {
      printf("curre=%10.1f err=%10.1f toterr=%10.1f PP=%i ",
             (double)currerror, (double)error, (double)totalerror, PP[j]);
      PrintVector(Vector(TS, j), VectorSize(TS), 1);
      }
    }

  return( totalerror );
}


/* ----------------------------------------------------------------- */


static void TentativePartitioning(TRAININGSET*    TS,
                                  PARTITIONING*   P,
                                  CLUSTERINFO     CI[],
                                  PREPARTITIONING PP[],
                                  int             Pindex)
{
  llong E1, E2;
  int   j;

  for( j = FirstVector(P,Pindex); ! EndOfPartition(j); j = NextVector(P,j) )
    {
    E1 = VectorDist(Vector(TS, j), CI[Pindex].Centroid1, VectorSize(TS));
    E2 = VectorDist(Vector(TS, j), CI[Pindex].Centroid2, VectorSize(TS));
    if( E1 < E2 )
      {
      PP[j] = 1;
      }
    else
      {
      PP[j] = 2;
      }
    }
}


/* ----------------------------------------------------------------- */


static void TentativeCentroids(TRAININGSET*    TS,
                               PARTITIONING*   P,
                               CLUSTERINFO     CI[],
                               PREPARTITIONING PP[],
                               int             Pindex)
{
  int      j, k;
  llong*   C1 = AllocateCounterVector(VectorSize(TS));
  llong*   C2 = AllocateCounterVector(VectorSize(TS));
  int      freq1 = 0, freq2 = 0;

  /* Initialize local centroid counters. */
  for( k = 0; k < VectorSize(TS); k++ )
    {
    C1[k] = 0LL;
    C2[k] = 0LL;
    }
  for( j = FirstVector(P,Pindex); ! EndOfPartition(j); j = NextVector(P,j) )
    {
    if( PP[j] == 1 )
      {
      for( k = 0; k < VectorSize(TS); k++ )
        {
        C1[k] += VectorScalar(TS, j, k) * VectorFreq(TS, j);
        }
      freq1 += VectorFreq(TS, j);
      }
    else /* PP[j] == 2 */
      {
      for( k = 0; k < VectorSize(TS); k++ )
        {
        C2[k] += VectorScalar(TS, j, k) * VectorFreq(TS, j);
        }
      freq2 += VectorFreq(TS, j);
      }
    }

  /* ...and centroids. */
  for( k = 0; k < VectorSize(TS); k++ )
    {
    if( freq1 > 0 ) CI[Pindex].Centroid1[k] = (VECTORELEMENT)round((double)C1[k]/(double)freq1);
    if( freq2 > 0 ) CI[Pindex].Centroid2[k] = (VECTORELEMENT)round((double)C2[k]/(double)freq2);
    }
  deallocate(C1);
  deallocate(C2);
}


/* ----------------------------------------------------------------- */


static void PerformLocalGLA(TRAININGSET*    TS,
                            CODEBOOK*       CB,
                            PARTITIONING*   P,
                            CLUSTERINFO     CI[],
                            PREPARTITIONING PP[],
                            int             Pindex)
{
  int      j, k;
  int      nchanged = 0;
  llong    E1, E2;
  int      iter;
  int      toPartition; /* 1 or 2 */

  llong*   C1 = AllocateCounterVector(VectorSize(TS));
  llong*   C2 = AllocateCounterVector(VectorSize(TS));
  int      freq1 = 0, freq2 = 0;

  if( Value(ShowProgress) >= 3 )
    {
    printf("Local GLA iter:");
    }

  /* Initialize local centroid counters. */
  for( k = 0; k < VectorSize(TS); k++ )
    {
    C1[k] = 0LL;
    C2[k] = 0LL;
    }
  for( j = FirstVector(P,Pindex); ! EndOfPartition(j); j = NextVector(P,j) )
    {
    if( PP[j] == 1 )
      {
      for( k = 0; k < VectorSize(TS); k++ )
        {
        C1[k] += VectorScalar(TS, j, k) * VectorFreq(TS, j);
        }
      freq1 += VectorFreq(TS, j);
      }
    else /* PP[j] == 2 */
      {
      for( k = 0; k < VectorSize(TS); k++ )
        {
        C2[k] += VectorScalar(TS, j, k) * VectorFreq(TS, j);
        }
      freq2 += VectorFreq(TS, j);
      }
    }

  /* ...and centroids. */
  for( k = 0; k < VectorSize(TS); k++ )
    {
    if( freq1 > 0 ) CI[Pindex].Centroid1[k] = (VECTORELEMENT)round((double)C1[k]/(double)freq1);
    if( freq2 > 0 ) CI[Pindex].Centroid2[k] = (VECTORELEMENT)round((double)C2[k]/(double)freq2);
    }

  iter = 0;
  nchanged = 1;
  while( nchanged > 0 && iter < Value(LocalGLAIterations) )
    {
    /* Partitioning. */
    nchanged = 0;
    for( j = FirstVector(P,Pindex); ! EndOfPartition(j); j = NextVector(P,j) )
      {
      E1 = VectorDist(Vector(TS, j), CI[Pindex].Centroid1, VectorSize(TS));
      E2 = VectorDist(Vector(TS, j), CI[Pindex].Centroid2, VectorSize(TS));
      toPartition = ( E1 < E2 ? 1 : 2 );
      if( toPartition != PP[j] )
        {
        PP[j] = toPartition;
        nchanged++;
        if( toPartition == 1 )
          {
          for( k = 0; k < VectorSize(TS); k++ )
            {
            C1[k] += VectorScalar(TS, j, k) * VectorFreq(TS, j);
            C2[k] -= VectorScalar(TS, j, k) * VectorFreq(TS, j);
            }
          freq1 += VectorFreq(TS, j);
          freq2 -= VectorFreq(TS, j);
          }
        else /* toPartition == 2 */
          {
          for( k = 0; k < VectorSize(TS); k++ )
            {
            C1[k] -= VectorScalar(TS, j, k) * VectorFreq(TS, j);
            C2[k] += VectorScalar(TS, j, k) * VectorFreq(TS, j);
            }
          freq1 -= VectorFreq(TS, j);
          freq2 += VectorFreq(TS, j);
          }
        }
      }

    /* Calculate new codevectors if necessary. */
    if( nchanged > 0 )
      {
      for( k = 0; k < VectorSize(TS); k++ )
        {
        if( freq1 > 0 ) CI[Pindex].Centroid1[k] = (VECTORELEMENT)round((double)C1[k]/(double)freq1);
        if( freq2 > 0 ) CI[Pindex].Centroid2[k] = (VECTORELEMENT)round((double)C2[k]/(double)freq2);
        }
      }

    iter++;
    }

  if( Value(ShowProgress) >= 3 )
    {
    printf("%7i\n", iter);
    }
  deallocate(C1);
  deallocate(C2);
}


/* ----------------------------------------------------------------- */


static void SplitClusterTentatively(TRAININGSET*    TS,
                                    CODEBOOK*       CB,
                                    PARTITIONING*   P,
                                    CLUSTERINFO     CI[],
                                    PREPARTITIONING PP[],
                                    int             Pindex)
{
  if( Value(Hyperplane) == YES )
    {
    SplitByHyperPlane(TS, CB, P, CI, PP, Pindex);
    }
  else
    {
    /* Try to find new codevectors for splited partition */
    SelectTwoCodevectors(TS, CB, P, CI, Pindex);
    TentativePartitioning(TS, P, CI, PP, Pindex);
    }

  if( Value(LocalGLA) == YES )
    {
    PerformLocalGLA(TS, CB, P, CI, PP, Pindex);
    }

  /* Now we have partitioning in PP and
     codevectors in CI[Pindex].Centroid1 and CI[Pindex].Centroid2. */
}


/* ----------------------------------------------------------------- */

#define sign(a) ( (a) < 0 ? -1 : ((a) == 0 ? 0 : 1) )

int CompareClusters(void* c1, void* c2, void* info)
{
  int result = sign( ((CLUSTERINFO*)c1)->DistortionValue -
                     ((CLUSTERINFO*)c2)->DistortionValue);
  if( result == 0  )
    {
    result = sign( ((CLUSTERINFO*)c1)->whoami - ((CLUSTERINFO*)c2)->whoami );
    }
/*
  printf("c1=%.1f c2=%.1f %i %i res=%i\n",
         (double)((CLUSTERINFO*)c1)->DistortionValue,
         (double)((CLUSTERINFO*)c2)->DistortionValue,
         ((CLUSTERINFO*)c1)->whoami,
         ((CLUSTERINFO*)c2)->whoami,
         result);
*/
  return( result );
}


/* ----------------------------------------------------------------- */


static void EvaluateCluster(TRAININGSET*    TS,
                            CODEBOOK*       CB,
                            PARTITIONING*   P,
                            CLUSTERINFO     CI[],
                            PREPARTITIONING PP[],
                            int             Pindex)
/* For deciding if cluster should be split. */
{
  int   maxj1, maxj2;

  if( UniqueVectors(P, Pindex) < 2 ) /* Is there something to split? */
    { /* No. */
    CI[Pindex].DistortionValue =
    CI[Pindex].SquareError = 0LL;
    }
  else
    { /* Yes. */
    switch( Value(PartitionErrorType) )
      {
      case SplitAllPartitions:
        {
        /* This is just for producing output information.
           Binary splitting does not need this. */
        CI[Pindex].DistortionValue =
        CI[Pindex].SquareError =
          PartitionError(TS, P, Pindex, &Node(CB, Pindex), MAXLLONG);
        break;
        }
      case SquareError:
        {
        CI[Pindex].DistortionValue =
        CI[Pindex].SquareError =
          PartitionError(TS, P, Pindex, &Node(CB, Pindex), MAXLLONG);
        break;
        }
      case SizeOfPartition:
        {
        /* If there is only one kind of vectors in partition
           we can not split it. So, let the distortion be zero. */
        CI[Pindex].DistortionValue =
          ( UniqueVectors(P, Pindex) > 1 ? CCFreq(P, Pindex): 0L );
        break;
        }
      case FurthestTwoVectors:
        {
        /* TODO: Store maxj1 and maxj2 for later use. */

        CI[Pindex].DistortionValue =
        CI[Pindex].FurthestDistance =
          TwoFurthestVectorsInPartition(TS, CB, P, Pindex, &maxj1, &maxj2);
        CI[Pindex].TrVector1 = maxj1;
        CI[Pindex].TrVector2 = maxj2;
        break;
        }
      case ThirdMoment:
        {
        CI[Pindex].DistortionValue =
        CI[Pindex].ThirdMoment =
          PartitionThirdMoment(TS, P, Pindex, &Node(CB, Pindex), VectorSize(TS), MAXLLONG);
        break;
        }
      case SplitRandom:
        {
        /* This is not an even distribution!!! */
        CI[Pindex].SquareError =
          PartitionError(TS, P, Pindex, &Node(CB, Pindex), MAXLLONG);
        CI[Pindex].DistortionValue =
          (llong)( drand() * (double)CI[Pindex].SquareError );
        break;
        }
      case OptimalSelection:
        {
        SplitClusterTentatively(TS, CB, P, CI, PP, Pindex);
        CI[Pindex].DistortionValue =
        CI[Pindex].DistortionChange =
          TentativeError(TS, CB, P, CI, PP, Pindex);
        break;
        }
      default:
        {
        break;
        }
      }
    }

  if( Value(ShowProgress) >= 4 )
    {
    printf("EvaluateCluster:%4i CBfreq=%4i freq=%4i uniq=%4i dist=%.1f\n",
           Pindex, VectorFreq(CB, Pindex), CCFreq(P, Pindex),
           UniqueVectors(P, Pindex), (double)CI[Pindex].DistortionValue);
    }
}


/*-------------------------------------------------------------------*/


static void Remapping(TRAININGSET*    TS,
                      CODEBOOK*       CB,
                      PARTITIONING*   P,
                      CLUSTERINFO     CI[],
                      PREPARTITIONING PP[],
                      BINTREE*        Tree,
                      int             codevectorindices[],
                      int             codevectorcount,
                      DISTANCETYPE    disttype)
{
  int      i, j, nextj, prev, part;
  int      best;
  llong    olddist, newdist;
  YESNO* reevaluate = (YESNO*)allocate(BookSize(CB) * sizeof(YESNO));

  assert( codevectorcount == 2 );

  for( i = 0; i < BookSize(CB); i++ )
    {
    reevaluate[i] = NO;
    }

  for( part = 0; part < BookSize(CB); part++ )
    {
    j = prev = FirstVector(P, part);
    while (! EndOfPartition(j))
      {
      if( UniqueVectors(P,part) >= 2 )
        {
        best = -1;
        olddist = VectorDistance(Vector(CB, part),
                               Vector(TS, j),
                               VectorSize(TS),
                               MAXLLONG,
                               disttype);

        for( i = 0; i < codevectorcount; i++)
          {
          newdist = VectorDistance(Vector(CB, codevectorindices[i]),
                                   Vector(TS, j),
                                   VectorSize(TS),
                                   olddist,
                                   disttype);
          if( newdist < olddist )
            {
            best = codevectorindices[i];
            olddist = newdist;
            }
          } /* for( i = 0; i < codevectorcount; i++) */

        if( best < 0 )
          {
          prev = j;
          j = NextVector(P, j);
          }
        else
          {
          if( reevaluate[part] == NO )
            {
            if( Value(PartitionErrorType) != SplitAllPartitions &&
                part != codevectorindices[0] &&
                part != codevectorindices[1] )
              {
              if( DeleteNodeFromBintree(Tree, &(CI[part]), NULL) == NULL )
                {
                printf("Del: j=%i map=%i best=%i\n", j, part, best);
                printf("dist=%f\n", (double)CI[part].DistortionValue);
    /*             PrintBintree(Tree); */
                }
              }
            reevaluate[part] = YES;
            }

          nextj = NextVector(P, j);
          ChangePartitionFast(TS, P, best, j, prev);
          j = nextj;
          }
        }
      }
    }

  for( i = 0; i < BookSize(CB); i++ )
    {
    if( reevaluate[i] == YES )
      {
      PartitionCentroid(P, i, &Node(CB, i));
      if( i != codevectorindices[0] && i != codevectorindices[1] )
        {
        EvaluateCluster(TS, CB, P, CI, PP, i);
        if( Value(PartitionErrorType) != SplitAllPartitions )
          {
          InsertToBintree(Tree, &(CI[i]), NULL);
          }
        }
      }
    }

  deallocate(reevaluate);
}


/* ----------------------------------------------------------------- */


static void SplitPartition(TRAININGSET*    TS,
                           CODEBOOK*       CB,
                           PARTITIONING*   P,
                           CLUSTERINFO     CI[],
                           PREPARTITIONING PP[],
                           BINTREE*        Tree,
                           int             Vold,
                           int             Vnew)
{
  int j;
  int prev;
  int nextj;

  if( Value(ShowProgress) >= 3 )
    {
    printf("PSize:%7i ", VectorFreq(CB, Vold));
    }

  /* Make partitioning for cluster based on PP. */
  j = prev = FirstVector(P, Vold);
  while( ! EndOfPartition(j) )
    {
    if( PP[j] == 1 )
      {
      prev = j;
      j = NextVector(P, j);
      }
    else
      {
      nextj = NextVector(P, j);
      /* Experimenting faster variant (PF) */
      if ( Value(NewVectors) == FasterVariant )
        {
        ChangePartitionFast(TS, P, Vnew, j, prev); 
        }
      else
        {
        ChangePartition(TS, P, Vnew, j);
        }
      j = nextj;
      }
    }

  /* Use precalculated centroids. */
  CopyVector(CI[Vold].Centroid1, Vector(CB, Vold), VectorSize(CB));
  CopyVector(CI[Vold].Centroid2, Vector(CB, Vnew), VectorSize(CB));
  VectorFreq(CB, Vold) = CCFreq(P, Vold);
  VectorFreq(CB, Vnew) = CCFreq(P, Vnew);

/*
printf("new freqs:%7i %7i XXX\n", VectorFreq(CB, Vold), VectorFreq(CB, Vnew));
PrintVector(Vector(CB, Vold), VectorSize(CB), 1);
PrintVector(Vector(CB, Vnew), VectorSize(CB), 1);
*/

  if( Value(PartitionRemapping) == YES )
    {
    int clusters[2] = { Vold, Vnew };
    Remapping(TS, CB, P, CI, PP, Tree, clusters, 2, EUCLIDEANSQ);
    }

  if( Value(ShowProgress) >= 3 )
    {
    printf("new freqs:%7i %7i\n", VectorFreq(CB, Vold), VectorFreq(CB, Vnew));
    }
}


/* ----------------------------------------------------------------- */


static YESNO SplitCluster(TRAININGSET*    TS,
                          CODEBOOK*       CB,
                          PARTITIONING*   P,
                          CLUSTERINFO     CI[],
                          PREPARTITIONING PP[],
                          BINTREE*        Tree,
                          int             Pindex)
/* Split partition Pindex to Pindex and BookSize(CB).
   Returns YES if success, NO otherwise. */
{
  if( Value(ShowProgress) >= 4 )
    {
    printf("SplitCluster: start Pindex=%i\n", Pindex);
    }

  /* Is there something to split? */
  if( UniqueVectors(P, Pindex) <= 1 )
    {
    /* No. */
    return( NO );
    }

  if( Value(PartitionErrorType) != OptimalSelection )
    {
    SplitClusterTentatively(TS, CB, P, CI, PP, Pindex);
    }

  /* Now we have partitioning in PP and
     codevectors in CI[Pindex].Centroid1 and CI[Pindex].Centroid2. */

  IncreaseCodebookSize(CB, BookSize(CB)+1);
  IncreaseNumberOfPartitions(P, PartitionCount(P)+1);

  /* Split partition based on prepartitioning (PP). */
  SplitPartition(TS, CB, P, CI, PP, Tree, Pindex, BookSize(CB)-1);

  /* Recalculate centroids. Added by JS for v.0.37 */
  PartitionCentroid(P, Pindex, &Node(CB, Pindex));
  PartitionCentroid(P, BookSize(CB)-1, &Node(CB, BookSize(CB)-1));

  /* Calculate distortions of 'Pindex' and 'BookSize(CB)-1'. */
  EvaluateCluster(TS, CB, P, CI, PP, Pindex);
  EvaluateCluster(TS, CB, P, CI, PP, BookSize(CB)-1);

  if( Value(PartitionErrorType) != SplitAllPartitions )
    {
    InsertToBintree(Tree, &(CI[Pindex]), NULL);
    InsertToBintree(Tree, &(CI[BookSize(CB)-1]), NULL);
    }

  if( Value(ShowProgress) >= 4 )
    {
    printf("SplitCluster: OK P=%4i PF=%5i P2=%4i P2F=%5i\n",
           Pindex, VectorFreq(CB, Pindex),
           BookSize(CB)-1, VectorFreq(CB, BookSize(CB)-1));
    }
  return( YES );
}


/* ----------------------------------------------------------------- */


static int SelectCluster(BINTREE* Tree)
{
  CLUSTERINFO* cluster = DeleteMaximumFromBintree(Tree);
  assert( cluster != NULL );

  if( Value(ShowProgress) >= 4 )
    {
    printf("SelectedCluster:whoami=%i dv=%f\n",
           cluster->whoami,
           (double)cluster->DistortionValue);
    PrintBintree(Tree);
    }
  return( cluster->whoami );
}


/* ----------------------------------------------------------------- */


static void IterateGLA(TRAININGSET*    TS,
                       CODEBOOK*       CB,
                       SASchedule*     SAS,
                       PARTITIONING*   P,
                       CLUSTERINFO     CI[],
                       PREPARTITIONING PP[],
                       BINTREE*        Tree)
{
  int   i;

  if( Value(GLAIterations) > 0 )
    {
    for( i = 0; i < Value(GLAIterations); i++ )
      {
      GenerateOptimalPartitioning(TS, CB, P);
      GenerateOptimalCodebook(TS, CB, P);
      FillEmptyPartitions(TS, CB, P);
      }

    /* Reevaluate clusters. */
    FreeBintree(Tree);
    for( i = 0; i < BookSize(CB); i++ )
      {
      EvaluateCluster(TS, CB, P, CI, PP, i);
      if( Value(PartitionErrorType) != SplitAllPartitions )
        {
        InsertToBintree(Tree, &(CI[i]), NULL);
        }
      }
    }
}


/* ================================================================= */

#define min(a,b) ((a) < (b) ? (a) : (b))

static void BinarySplit(TRAININGSET*    TS,
                        CODEBOOK*       CB,
                        SASchedule*     SAS,
                        PARTITIONING*   P,
                        CLUSTERINFO     CI[],
                        PREPARTITIONING PP[],
                        BINTREE*        Tree,
                        int             ResultCBSize)
{
  int i;
  int loopsize = min(BookSize(CB), ResultCBSize - BookSize(CB));
  int nsplit = 0;

  for( i = 0; i < BookSize(CB) && nsplit < loopsize; i++ )
    {
    if( SplitCluster(TS, CB, P, CI, PP, Tree, i) )
      {
      nsplit++;
      }
    }

  IterateGLA(TS, CB, SAS, P, CI, PP, Tree);
}


/* ================================================================= */


static void InitializeSPLIT(TRAININGSET*     TS,
                            CODEBOOK*        CB,
                            PARTITIONING*    P,
                            CLUSTERINFO**    CI,
                            PREPARTITIONING* PP[],
                            BINTREE*         Tree,
                            int              ResultCBSize)
{
  int i, j;

  if( BookSize(CB) < 1 )
    {
    /* We have not an initial codebook and initial partition. */
    printfe("ERROR: Split needs an initial codebook (BookSize(CB)==%i)\n",
           BookSize(CB));
    exit( -1 );
    }

  if( Value(ShowProgress) >= 4 )
    {
    printf("InitializeSPLIT starts. CBsize=%i ResCB=%i\n", BookSize(CB), ResultCBSize);
    }

  SetAllocatedCodebookSize(CB, ResultCBSize);
  SetNumberOfAllocatedPartitions(P, ResultCBSize);

  *CI = (CLUSTERINFO*)allocate(ResultCBSize * sizeof(CLUSTERINFO));
  *PP = (PREPARTITIONING*)allocate(BookSize(TS) * sizeof(PREPARTITIONING));
  for( j = 0; j < BookSize(TS); j++ ) /* Necessary initializing? */
    {
    (*PP)[j] = 0;
    }

  /* For binary tree, set identification labels and allocate memory. */
  for( i = 0; i < ResultCBSize; i++ )
    {
    (*CI)[i].whoami = i;
    (*CI)[i].Centroid1 = CreateEmptyVector(VectorSize(CB));
    (*CI)[i].Centroid2 = CreateEmptyVector(VectorSize(CB));
    }

  InitBintree(Tree, CompareClusters);
  for( i = 0; i < BookSize(CB); i++ )
    {
    if( Value(ShowProgress) >= 4 )
      {
      printf("Initialize %i\n", i);
      }

    /* Ensure codevector frequencies. */
    VectorFreq(CB, i) = CCFreq(P, i);

    EvaluateCluster(TS, CB, P, *CI, *PP, i);
    if( Value(PartitionErrorType) != SplitAllPartitions )
      {
      InsertToBintree(Tree, &((*CI)[i]), NULL);
      }
    }

  if( Value(ShowProgress) >= 4 )
    {
    printf("InitializeSPLIT done.\n");
    }
}


/* ----------------------------------------------------------------- */


static void ShutDownSPLIT(CLUSTERINFO     CI[],
                          PREPARTITIONING PP[],
                          BINTREE*        Tree,
                          int             ResultCBSize)
{
  int i;

  for( i = 0; i < ResultCBSize; i++ )
    {
    FreeVector(CI[i].Centroid1);
    FreeVector(CI[i].Centroid2);
    }
  deallocate(CI);
  deallocate(PP);
  FreeBintree(Tree);
}


/*====================================================================*/


void Split(TRAININGSET*  TS,
           CODEBOOK*     CB,
           SASchedule*   SAS,
           PARTITIONING* P,
           int           ResultCBSize)
/* Assume: BookSize(TS) >= 1.
           BookSize(CB) >= 1.
           Partition P is current.
*/
{
  CLUSTERINFO*     CI;
  BINTREE          Tree;
  PREPARTITIONING* PP;
  int              cluster;

  InitializeSPLIT(TS, CB, P, &CI, &PP, &Tree, ResultCBSize);

  while( BookSize(CB) < ResultCBSize )
    {
    if( Value(PartitionErrorType) == SplitAllPartitions )
      {
      BinarySplit(TS, CB, SAS, P, CI, PP, &Tree, ResultCBSize);
      }
    else
      {
      cluster = SelectCluster(&Tree);

/* printf("SF= %5i \n", VectorFreq(CB, cluster)); */

      SplitCluster(TS, CB, P, CI, PP, &Tree, cluster);
      IterateGLA(TS, CB, SAS, P, CI, PP, &Tree);
      }

    if( Value(ShowProgress) >= 2 )
      {
      printf("Split CB:%5i %9.4f\n",
             BookSize(CB), PrintableError(SE2MSE(TS, CB, CI), CB));
      if( Value(ShowProgress) >= 4 )
        {
        int i;
        for( i = 0; i < BookSize(CB); i++ )
          {
          printf("%3i = %9.4f\n", i, PrintableError(PartitionError(TS, P, i, &Node(CB, i), MAXLLONG), CB));
          }
        printf("MSE = %9.4f\n", PrintableError(AverageErrorCBFast(TS, CB, P, MSE), CB));
        }
      fflush(stdout);
      }
    }
  ShutDownSPLIT(CI, PP, &Tree, ResultCBSize);
}


/*====================================================================*/


static SEARCHTREE* CreateSTNode(CODEBOOK* CB)
{
  SEARCHTREE* t;

  t = allocate( sizeof(SEARCHTREE) );
  t->Leaf = YES;
  t->Left = NULL;
  t->Right = NULL;
  t->ClusterIndex = -1;
  t->Centroid = CreateEmptyVector(VectorSize(CB));
  return t;
}


/* ----------------------------------------------------------------- */


static void SplitSTNode(CODEBOOK* CB, SEARCHTREE* SplitNode,
                        int changed, int new)
{
  SplitNode->Left = CreateSTNode(CB);
  SplitNode->Right = CreateSTNode(CB);
  CopyVector(Vector(CB, changed), SplitNode->Left->Centroid, VectorSize(CB));
  CopyVector(Vector(CB, new), SplitNode->Right->Centroid, VectorSize(CB));
  SplitNode->Leaf = NO;
}


/* ----------------------------------------------------------------- */


static void FindClusterIndex(CODEBOOK* CB, SEARCHTREE* LeafNode)
{
  int i = 0;

  while (i < BookSize(CB))
      {
      if EqualVectors(Vector(CB, i), LeafNode->Centroid, VectorSize(CB))
          {
          LeafNode->ClusterIndex = i;
          return;
          }
      i++;
      }
  printf("WARNING: FindClusterIndex failed\n");
  LeafNode->ClusterIndex = 0;
}


/* ----------------------------------------------------------------- */


SEARCHTREE* GenerateSearchTreeForCodebook(CODEBOOK* CB)
{
  int cluster;
  CODEBOOK splitCB;
  PARTITIONING P;
  SEARCHTREE** LeafPointers;        /* pointers to leaf nodes */
  SEARCHTREE* ST;
  BINTREE Tree;
  CLUSTERINFO* CI;
  PREPARTITIONING* PP;

  LeafPointers = allocate(BookSize(CB) * sizeof(SEARCHTREE*));
  CreateNewCodebook(&splitCB, 1, CB);
  CreateNewPartitioning(&P, CB, 1);
  GenerateOptimalCodebook(CB, &splitCB, &P);
  InitializeSPLIT(CB, &splitCB, &P, &CI, &PP, &Tree, BookSize(CB));
  ST = CreateSTNode(CB);
  CopyVector(Vector(CB, 0), ST->Centroid, VectorSize(CB));
  LeafPointers[0] = ST;

  cluster = 0;
  while (BookSize(&splitCB) < BookSize(CB))
      {
      /* cluster = SelectCluster(&Tree); */
      if ( CCFreq(&P, cluster) == 1)
          cluster++;
      else
          {
          SplitCluster(CB, &splitCB, &P, CI, PP, &Tree, cluster);
          SplitSTNode(&splitCB, LeafPointers[cluster], cluster,
                      BookSize(&splitCB)-1);
          LeafPointers[BookSize(&splitCB)-1] = LeafPointers[cluster]->Right;
          LeafPointers[cluster] = LeafPointers[cluster]->Left;
          if (CCFreq(&P, cluster) == 1)
              FindClusterIndex(CB, LeafPointers[cluster]);
          if (CCFreq(&P, BookSize(&splitCB)-1) == 1)
              FindClusterIndex(CB, LeafPointers[BookSize(&splitCB)-1]);
          }
      if (cluster >= BookSize(&splitCB)) cluster = 0;
      }

  deallocate(LeafPointers);
  FreeCodebook(&splitCB);
  FreePartitioning(&P);
  ShutDownSPLIT(CI, PP, &Tree, BookSize(CB));

  return ST;
}
