/*-------------------------------------------------------------------*/
/* SEARCH.C        Pasi Fr„nti                                       */
/*                 Juha Kivij„rvi  (PCA-tree)                        */
/*                                                                   */
/* Search structures for codebook (CB):                              */
/*   - Search using vector-mean                                      */
/*   - PCA-based binary tree (PCA-tree)                              */
/*                                                                   */
/*-------------------------------------------------------------------*/

#define ProgName           "SS"
#define VersionNumber      "Version 0.03a"
#define LastUpdated        "28.7.99"
#define SearchStructID     "VQ SEARCH STRUCTURE"
#define FileFormatVersion  "1.0"

/* ----------------------------------------------------------------- */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <values.h>

#include "bintree.h"
#include "memctrl.h"
#include "file.h"
#include "cb.h"
#include "owntypes.h"
#include "search.h"
#include "sortcb.h"
#include "split.h"


#define EPSILON 0.0000001

/*=====================  Vector mean search  =========================*/


static int NearestVectorMeanByBinarySearch(CODEBOOK* CB, int value)
{
  int begin=0;
  int end=BookSize(CB)-1;
  int middle;

  while(begin!=end)
     {
     middle = (begin+end) >> 1;
     if( value<=VectorMean(CB,middle) ) end = middle;
     else begin = middle + 1;
     }
  return(begin);
}


/* ----------------------------------------------------------------- */


int FindNearestUsingVectorMean(BOOKNODE* v, CODEBOOK* CB, int* tested)
{
  int	Begin,i;
  int   tests = 1;
  int   MinIndex;
  float MinMean,MaxMean;
  float kroot=sqrt(VectorSize(CB));
  llong MinError;
  llong e;

  Begin = MinIndex = NearestVectorMeanByBinarySearch(CB,v->vmean);
  MinError = VectorDist(Vector(CB,MinIndex), v->vector, VectorSize(CB));

/*printf("\nTraining vector: "); PrintVector(v->vector,VectorSize(CB),1);
  printf("Initial=%i  Mean=%i  MinError=%li  Search area = %.1f..%.1f\n",
          MinIndex,v->vmean,(long)MinError,MinMean,MaxMean); */

  i = Begin-1;
  MinMean  = v->vmean - sqrt(MinError)/kroot - 1;
  while( VectorMean(CB,i)>MinMean && i>=0 )
     {
     e = VectorDist(Vector(CB,i), v->vector, VectorSize(CB));
     if( e<MinError )
        {
        MinError=e; MinIndex=i;
        MinMean  = v->vmean - sqrt(MinError)/kroot - 1;
        }
  /* printf("Low: Mean=%i dist=%li Min=%i\n", VectorMean(CB,i), (long)e, MinIndex); */
     i--; tests++;
     }

  i = Begin+1;
  MaxMean  = v->vmean + sqrt(MinError)/kroot - 1;
  while( VectorMean(CB,i)<MaxMean && i<BookSize(CB) )
	 {
     e = VectorDist(Vector(CB,i), v->vector, VectorSize(CB));
     if( e<MinError )
        {
        MinError=e; MinIndex=i;
        MaxMean  = v->vmean + sqrt(MinError)/kroot - 1;
        }
  /* printf("High: Mean=%i dist=%li Min=%i\n", VectorMean(CB,i), (long)e, MinIndex); */
     i++; tests++;
     }

  *tested = tests;
  return( MinIndex );
}


/* ----------------------------------------------------------------- */


void GenerateOptimalPartitionUsingVmeans(TRAININGSET*  TS,
                                         CODEBOOK*     CB,
                                         PARTITIONING* P)
{
  int nearest,tested,i;
  long total=0;
/*int nearest0;
  llong error; */

  CalculateVectorMeans(CB);
  SortCodebook(CB, VECTOR_MEAN);
/* PrintCodebook(CB); */
  for( i=0; i<BookSize(TS); i++ )
     {
/*   nearest0 = FindNearestVector( &Node(TS,i), CB, &error, Map(P,i), EUCLIDEANSQ); */
     nearest = FindNearestUsingVectorMean( &Node(TS,i), CB, &tested );
     total += tested;
/*   if( nearest!=nearest0 &&
        VectorDist( Vector(TS,i), Vector(CB,nearest0), VectorSize(CB) ) !=
        VectorDist( Vector(TS,i), Vector(CB,nearest),  VectorSize(CB) ) )
      { printf("ERROR: nearest(%i)<>nearest0(%i)\n",nearest,nearest0);
        printf("       d1=%li, d0=%li \n",
         (long) VectorDist( Vector(TS,i), Vector(CB,nearest),  VectorSize(CB) ),
         (long) VectorDist( Vector(TS,i), Vector(CB,nearest0), VectorSize(CB) ) );
        ExitProcessing(-1); } */
     if(nearest!=Map(P,i))  ChangePartition( TS, P, nearest, i);
     }
  printf("Tested %7.1f vectors on average\n", total/BookSize(TS) + 0.5 );
}


/*========================  Search tree  =============================*/


SEARCHTREE* GenerateSearchTree(CODEBOOK* CB)
{
  int i;

  for (i=0; i<BookSize(CB); i++)
      VectorFreq(CB, i) = 1;
  return GenerateSearchTreeForCodebook(CB);
}


/* ----------------------------------------------------------------- */


void FreeSearchTree(SEARCHTREE* ST)
{
  if (!(ST->Leaf))
      {
      FreeSearchTree(ST->Left);
      FreeSearchTree(ST->Right);
      }
  FreeVector(ST->Centroid);
  deallocate(ST);
}


/* ----------------------------------------------------------------- */


static void PrintSearchTreeIndented(SEARCHTREE* ST,
                                    int VectorSize,
                                    int spaces)
{
  int i;

  for (i=0; i<spaces; i++) printf(" ");
  PrintVector(ST->Centroid, VectorSize, 1);
  if (ST->Leaf == NO)
      {
      PrintSearchTreeIndented(ST->Left, VectorSize, spaces + 2);
      PrintSearchTreeIndented(ST->Right, VectorSize, spaces + 2 );
      }
}


/* ----------------------------------------------------------------- */

void PrintSearchTree(SEARCHTREE* ST, int VectorSize)
{
   PrintSearchTreeIndented(ST, VectorSize, 0);
}


/* ----------------------------------------------------------------- */


int SearchFromTree(SEARCHTREE* ST, BOOKNODE* node, int VectorSize)
{
SEARCHTREE* tnode = ST;
llong dist1, dist2;

  while( tnode->Leaf == NO )
     {
     dist1 = VectorDist(node->vector, tnode->Left->Centroid, VectorSize);
     dist2 = VectorDistance(node->vector, tnode->Right->Centroid,
                            VectorSize, dist1, EUCLIDEANSQ);
     if (dist1 <= dist2)
         tnode = tnode->Left;
     else
         tnode = tnode->Right;
     }

return tnode->ClusterIndex;
}


/* ----------------------------------------------------------------- */


double RunGLATreeSearch(TRAININGSET*  TS,
                        CODEBOOK*     CB,
                        PARTITIONING* P,
                        int           iterations)
/* starts always from the codebook, supports only MSE */
{
  SEARCHTREE* ST;
  int i;
  double newerror = MAXDOUBLE, olderror;
  int icount = 0;
/*
  PARTITIONING tp; int wrong; int item;
  CreateNewPartitioning(&tp, TS, BookSize(CB));
*/
  if (iterations == 0) return AverageErrorForSolution(TS, CB, P, MSE);

  do
     {
     olderror = newerror;
     /* 1. Generate search tree: O(M*logM) */
     ST = GenerateSearchTree(CB);
     /* 2. Suboptimal partitioning: O(N*logM) */
     for(i=0; i<BookSize(TS); i++)
         {
         ChangePartition(TS, P, SearchFromTree(ST, &Node(TS, i),
                                               VectorSize(CB)), i);
         }

/*
     PrintSearchTree(ST, VectorSize(TS));

  GenerateOptimalPartitioning(TS, CB, &tp); wrong = 0;
  for (i=0; i<BookSize(TS); i++)
      if (VectorScalar(CB, Map(P, i), 0) != VectorScalar(CB, Map(&tp, i), 0)
          || VectorScalar(CB, Map(P, i), 1) != VectorScalar(CB, Map(&tp, i), 1))
          printf("vector (%d,%d) goes to (%d,%d) instead of (%d,%d)\n",
             VectorScalar(TS, i, 0), VectorScalar(TS, i, 1),
             VectorScalar(CB, Map(P, i), 0), VectorScalar(CB, Map(P, i), 1),
             VectorScalar(CB, Map(&tp, i), 0), VectorScalar(CB, Map(&tp, i), 1));

  printf("The following mappings exist:\n");
  for (i=0; i<BookSize(CB); i++)
  {   printf("*** ");
      PrintVector(Vector(CB,i), VectorSize(CB), 1);
  item = FirstVector(P, i);
  while (!(EndOfPartition(item)))
      PrintVector(Vector(TS,item), VectorSize(CB), 1), item = NextVector(P, item);
  } */
     /* 3. Optimal codebook: O(N) */
     GenerateOptimalCodebook(TS, CB, P);
     FillEmptyPartitions(TS, CB, P);

     newerror = AverageErrorForSolution(TS, CB, P, MSE);
     icount++;
     /* 4. Optional noise adding */
     /* AddNoiseToCodebook(CB, SAS); */
     }
  while ( olderror - newerror > EPSILON && icount < iterations);
  
  FreeSearchTree(ST);

  return newerror;
}

