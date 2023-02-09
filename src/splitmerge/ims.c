/*-------------------------------------------------------------------*/
/* IMS.C          Timo Kaukoranta                                    */
/*                                                                   */
/*                                                                   */
/* - Module of Iterative Merge and Split                             */
/*                                                                   */
/*-------------------------------------------------------------------*/

#define ProgName       "IMS"
#define VersionNumber  "V.0.09"
#define LastUpdated    "16.6.09"  /* PF */

/* ----------------------------------------------------------------- */
/*      Changelog:                                                   */
/*                                                                   */
/* PF   0.09  BugFix: CB size was wrongly initialized to double.     */
/*      0.08  Original version by TK (14.10.97)                      */
/*                                                                   */
/* ----------------------------------------------------------------- */

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <values.h>

#include "cb.h"
#include "memctrl.h"
#include "pnn.h"
#include "random.h"
#include "sa.h"
#include "split.h"

/* -------------------- Parameters for the IMS module --------------- */

static struct { int ShowProgress;
                int CodebookSize;
                int AlternatingOrder;
                int AlternatingAmount;
                int MinCodebookSize;
                int MaxCodebookSize;
                int StepSizeChange;
                int Iterations;
                int IMSGLAIterations;
                int RequireDistortionDecrease;
              } ModuleParameters;

#define Value(x)      (ModuleParameters.x)

/* These enum definitions are from the *.fac file. */
enum { SplitMerge, MergeSplit, FiftyFifty, BetterOne };

/*-----------------------------  M i s c  ----------------------------*/

typedef struct { llong LocalDist;
                 int   V1;
                 int   V2;
               } NEARESTTYPE;

#define round(a)      ((a) + 0.5)
#define fpositive(a)  ((a) > FLT_EPSILON)
#define printfe       printf("\nFile %s line %i ", __FILE__, __LINE__);printf

#if defined(min) || defined(max)
#undef min
#undef max
#endif
#define min(a,b)      ((a) < (b) ? (a) : (b))
#define max(a,b)      ((a) > (b) ? (a) : (b))

#define NewPNNSize(size,step)   ((size) - (step))
#define NewSplitSize(size,step) ((size) + (step))


/*
#define NewPNNSize(size)   ( Value(AlternatingAmount) == 0       \
                             ? size / 2                          \
                             : size - Value(AlternatingAmount) )
#define NewSplitSize(size) ( Value(AlternatingAmount) == 0       \
                             ? size * 2                          \
                             : size + Value(AlternatingAmount) )
*/

/* ================== Parameters for the IMS module ================ */


void SetIMSParameters(int ShowProgress,
                      int CodebookSize,
                      int AlternatingOrder,
                      int AlternatingAmount,
                      int MinCodebookSize,
                      int MaxCodebookSize,
                      int StepSizeChange,
                      int Iterations,
                      int IMSGLAIterations,
                      int RequireDistortionDecrease)
{
  Value(ShowProgress)              = ShowProgress;
  Value(CodebookSize)              = CodebookSize;
  Value(AlternatingOrder)          = AlternatingOrder;
  Value(AlternatingAmount)         = AlternatingAmount;
  Value(MinCodebookSize)           = MinCodebookSize;
  Value(MaxCodebookSize)           = MaxCodebookSize;
  Value(StepSizeChange)            = StepSizeChange;
  Value(Iterations)                = Iterations;
  Value(IMSGLAIterations)          = IMSGLAIterations;
  Value(RequireDistortionDecrease) = RequireDistortionDecrease;
}


/* ================================================================= */


static void InitializeIMS(TRAININGSET*  TS,
                          CODEBOOK*     CB,
                          PARTITIONING* P,
                          int*          MinCBSize,
                          int*          MaxCBSize,
                          int*          StepSize,
                          CODEBOOK*     CB2)
{
  /* What are the minimum and the maximum sizes of the codebook? */
  switch( Value(AlternatingOrder) )
    {
    case SplitMerge:
      {
      /* We first increase the size by split. */
      if( Value(AlternatingAmount) == 0 )
        {
        *MinCBSize = BookSize(CB);
        *MaxCBSize = 2 * BookSize(CB);
        *StepSize  = BookSize(CB);
        }
      else
        {
        *MinCBSize = BookSize(CB) - Value(AlternatingAmount);
        *MaxCBSize = BookSize(CB) + Value(AlternatingAmount);
        *StepSize  = Value(AlternatingAmount);
        }
      break;
      }
    case MergeSplit:
      {
      /* We first decrease the size by merge. */
      if( Value(AlternatingAmount) == 0 )
        {
        *MinCBSize = BookSize(CB) / 2;
        *MaxCBSize = BookSize(CB);
        *StepSize  = BookSize(CB) / 2;
        }
      else
        {
        *MinCBSize = BookSize(CB) - Value(AlternatingAmount);
        *MaxCBSize = BookSize(CB) + Value(AlternatingAmount);
        *StepSize  = Value(AlternatingAmount);
        }
      break;
      }
    case FiftyFifty:
    case BetterOne:
      {
      *MinCBSize = Value(MinCodebookSize);
      *MaxCBSize = Value(MaxCodebookSize);
      *StepSize  = Value(AlternatingAmount);
      break;
      }
    default:
      {
      printfe("ERROR: unknown AlternatingOrder=%i\n", Value(AlternatingOrder));
      exit( -1 );
      }
    }

  /* Check validity */
  if( BookSize(CB) < *MinCBSize ||
      *MaxCBSize < BookSize(CB) )
    {
    printf("ERROR: BookSize=%i Min=%i Max=%i\n",
           BookSize(CB), *MinCBSize, *MaxCBSize);
    exit( -1 );
    }

/* New variant: 16.6.09 (PF) */
  CreateNewCodebook(CB2, *MaxCBSize, TS);
  CopyCodebook(CB, CB2);
  IncreaseNumberOfPartitions(P, *MaxCBSize);
  PartitionCount(P) = BookSize(CB);

/* Old variant: wrong CB size (double) */
/*
  printf("CB=%i minCB=%i maxCB=%i\n", BookSize(CB), *MinCBSize, *MaxCBSize);
  IncreaseCodebookSize(CB2, *MaxCBSize);
  IncreaseNumberOfPartitions(P, *MaxCBSize);
  PartitionCount(P) = BookSize(CB2);
  CreateNewCodebook(CB2, *MaxCBSize, TS);
  CopyCodebook(CB, CB2);
  printf("CB2=%i minCB=%i maxCB=%i\n", BookSize(CB), *MinCBSize, *MaxCBSize); 
*/

}


/* ----------------------------------------------------------------- */


static void ShutDownIMS(CODEBOOK* CB2)
{
  FreeCodebook(CB2);
}



/*=========================  CODEBOOK HANDLING  ============================*/


static void PrintPartition(TRAININGSET* TS, PARTITIONING* P, int Pindex)
{
  int j, k;

  printf("PP%5i: Freq=%5i Unique=%5i\n",
         Pindex, CCFreq(P,Pindex), UniqueVectors(P, Pindex));
  for( j = FirstVector(P,Pindex); ! EndOfPartition(j); j = NextVector(P,j) )
    {
    printf("PP%5i: j=%5i freq=%5i: ", Pindex, j, VectorFreq(TS, j));
    for( k = 0; k < VectorSize(TS); k++ )
      {
      printf("%2i ", VectorScalar(TS, j, k));
      }
    printf("\n");
    }
  printf("\n");
}


/* ----------------------------------------------------------------- */


static llong PartitioningError(TRAININGSET*  TS,
                               CODEBOOK*     CB,
                               PARTITIONING* P)
{
  int   j;
  llong e, totalerror = 0LL;

  for( j = 0; j < BookSize(TS); j++ )
    {
    e = VectorDist(Vector(TS, j), Vector(CB, Map(P, j)), VectorSize(TS));
    totalerror += e * VectorFreq(TS, j);
    }
  return( totalerror );
}


/* ----------------------------------------------------------------- */


static double llong2doubleError(TRAININGSET* TS, llong error)
{
  assert( TotalFreq(TS) > 0 );

  return( PrintableError( (double)error /
                          ((double)TotalFreq(TS) *  (double)VectorSize(TS)),
                          TS ) );
}


/* ================================================================= */


static void IterateSplitMerge(TRAININGSET*  TS,
                              CODEBOOK*     CB,
                              SASchedule*   SAS,
                              PARTITIONING* P,
                              int           StepSize)
{
  int SplitSize, PNNSize;
  double d0 = 0, d1 = 0, d2 = 0;

  SplitSize = min(NewSplitSize(BookSize(CB), StepSize), BookSize(TS));
  PNNSize   = max(BookSize(CB),                         1);

  if( Value(ShowProgress) >= 3 )
    {
    printf("Split starting:\n");
    d0 = AverageErrorCB(TS, CB, MSE);
    }

  Split(TS, CB, SAS, P, SplitSize);

  if( Value(ShowProgress) >= 3 )
    {
    printf("PNN starting:\n");
    d1 = AverageErrorCB(TS, CB, MSE);
    }

  PairwiseNearestNeighbour(TS, CB, SAS, P, PNNSize);

  if( Value(ShowProgress) >= 3 )
    {
    d2 = AverageErrorCB(TS, CB, MSE);
    printf("S= %9.4f ", PrintableError(d1-d0, CB));
    printf("P= %9.4f ", PrintableError(d2-d1, CB));
    }
}


/* ----------------------------------------------------------------- */


static void IterateMergeSplit(TRAININGSET*  TS,
                              CODEBOOK*     CB,
                              SASchedule*   SAS,
                              PARTITIONING* P,
                              int           StepSize)
{
  int PNNSize, SplitSize;
  double d0 = 0, d1 = 0, d2 = 0;

  SplitSize = min(BookSize(CB),                      BookSize(TS));
  PNNSize   = max(NewPNNSize(BookSize(CB),StepSize), 1);

  if( Value(ShowProgress) >= 3 )
    {
    d0 = AverageErrorCB(TS, CB, MSE);
    }

  PairwiseNearestNeighbour(TS, CB, SAS, P, PNNSize);

  if( Value(ShowProgress) >= 3 )
    {
    d1 = AverageErrorCB(TS, CB, MSE);
    }

  Split(TS, CB, SAS, P, SplitSize);

  if( Value(ShowProgress) >= 3 )
    {
    d2 = AverageErrorCB(TS, CB, MSE);
    printf("P= %9.4f ", PrintableError(d1-d0, CB));
    printf("S= %9.4f ", PrintableError(d2-d1, CB));
    }
}


/* ----------------------------------------------------------------- */


static void IterateFiftyFifty(TRAININGSET*  TS,
                              CODEBOOK*     CB,
                              SASchedule*   SAS,
                              PARTITIONING* P,
                              int           MinCBSize,
                              int           MaxCBSize,
                              int           StepSize)
{
  int SplitSize, PNNSize;
  int r;

  assert( MinCBSize <= BookSize(CB) );
  assert( BookSize(CB) <= MaxCBSize );

  /* Try to keep BookSize as an average of min and max sizes. */
  r = irand(MinCBSize + 1, MaxCBSize);
  if( r <= BookSize(CB) )
    {
    PNNSize = max(NewPNNSize(BookSize(CB), StepSize), MinCBSize);
    PairwiseNearestNeighbour(TS, CB, SAS, P, PNNSize);
    }
  else
    {
    SplitSize = min(NewSplitSize(BookSize(CB), StepSize), MaxCBSize);
    Split(TS, CB, SAS, P, SplitSize);
    }
}


/* ----------------------------------------------------------------- */


static void IterateBetterOne(TRAININGSET*  TS,
                             CODEBOOK*     CB,
                             SASchedule*   SAS,
                             PARTITIONING* P,
                             int           MinCBSize,
                             int           MaxCBSize,
                              int           StepSize)
{
  int          SplitSize, PNNSize;
  double       Enow, Esplit, Epnn;
  double       Dsplit, Dpnn;
  CODEBOOK     CBsplit, CBpnn;
  PARTITIONING Psplit, Ppnn;

  assert( MinCBSize <= BookSize(CB) );
  assert( BookSize(CB) <= MaxCBSize );

  SplitSize = min(NewSplitSize(BookSize(CB), StepSize), MaxCBSize);
  PNNSize   = max(NewPNNSize(BookSize(CB), StepSize),   MinCBSize);

  CreateNewCodebook(&CBsplit, SplitSize,    TS);
  CreateNewCodebook(&CBpnn,   BookSize(CB), TS);

  CopyCodebook(CB, &CBsplit);
  CopyCodebook(CB, &CBpnn);

  CreateNewPartitioning(&Psplit, TS, SplitSize);
  CreateNewPartitioning(&Ppnn,   TS, BookSize(CB));

  CopyPartitioning(P, &Psplit);
  CopyPartitioning(P, &Ppnn);

  Split(TS, &CBsplit, SAS, &Psplit, SplitSize);
  PairwiseNearestNeighbour(TS, &CBpnn, SAS, &Ppnn, PNNSize);

  Enow   = AverageErrorCBFast(TS, CB,       P,       MSE);
  Esplit = AverageErrorCBFast(TS, &CBsplit, &Psplit, MSE);
  Epnn   = AverageErrorCBFast(TS, &CBpnn,   &Ppnn,   MSE);

  Dsplit = (Enow - Esplit);
  Dpnn   = (Epnn - Enow);

  if( Value(ShowProgress) >= 4 )
    {
    printf("DS/DP %9.4f%c %9.4f ",
           Dsplit,
           (Dsplit < Dpnn ? '<' : (Dsplit == Dpnn ? '=' : '>') ),
           Dpnn);
    fflush(stdout);
    }

  if( Dsplit > Dpnn || ( Dsplit == Dpnn && irand(0,1) ) )
    {
    CopyCodebook(&CBsplit, CB);
    CopyPartitioning(&Psplit, P);
    }
  else
    {
    CopyCodebook(&CBpnn, CB);
    CopyPartitioning(&Ppnn, P);
    }

  FreeCodebook(&CBsplit);
  FreeCodebook(&CBpnn);
  FreePartitioning(&Psplit);
  FreePartitioning(&Ppnn);
}


/* ================================================================= */


static void ChangeStep(int* StepSize)
{
  double change = (double)(*StepSize) * (double)Value(StepSizeChange) / 100.0;

  *StepSize = max(0, *StepSize - (int)round(change));
}


/* ----------------------------------------------------------------- */


static void RestoreCodebookSize(TRAININGSET*  TS,
                                CODEBOOK*     CB,
                                SASchedule*   SAS,
                                PARTITIONING* P,
                                int           InitialCBSize)
{
  while( BookSize(CB) != InitialCBSize )
    {
    if( Value(ShowProgress) >= 5 )
      {
      printf("Restore: %4i->%4i %9.4f\n",
             BookSize(CB), InitialCBSize,
             llong2doubleError(TS, PartitioningError(TS, CB, P)));
      fflush(stdout);
      }

    while( BookSize(CB) < InitialCBSize )
      {
      Split(TS, CB, SAS, P, InitialCBSize);
      }
    while( BookSize(CB) > InitialCBSize )
      {
      PairwiseNearestNeighbour(TS, CB, SAS, P, InitialCBSize);
      }
    }
}


/* ----------------------------------------------------------------- */


static void Iterate(TRAININGSET*  TS,
                    CODEBOOK*     CB,
                    SASchedule*   SAS,
                    PARTITIONING* P,
                    int           MinCBSize,
                    int           MaxCBSize,
                    int           StepSize)
{
  switch( Value(AlternatingOrder) )
    {
    case SplitMerge: IterateSplitMerge(TS, CB, SAS, P, StepSize);
                     break;
    case MergeSplit: IterateMergeSplit(TS, CB, SAS, P, StepSize);
                     break;
    case FiftyFifty: IterateFiftyFifty(TS, CB, SAS, P, MinCBSize, MaxCBSize, StepSize);
                     break;
    case BetterOne:  IterateBetterOne(TS, CB, SAS, P, MinCBSize, MaxCBSize, StepSize);
                     break;
    default:
      {
      printfe("ERROR: unknown AlternatingOrder=%i\n", Value(AlternatingOrder));
      exit( -1 );
      }
    }
}


/* ----------------------------------------------------------------- */


static void IterateGLA(TRAININGSET*  TS,
                       CODEBOOK*     CB,
                       PARTITIONING* P,
                       llong*        error)
{
  int    i = 0;
  double d0=0, d1=0;
  llong  e = MAXLLONG;

  if( Value(ShowProgress) >= 3 )
    {
    d0 = AverageErrorCB(TS, CB, MSE);
    }

  *error = PartitioningError(TS, CB, P);
  GenerateOptimalPartitioning(TS, CB, P);

  while( i < Value(IMSGLAIterations) && *error < e )
    {
    e = *error;
    GenerateOptimalCodebook(TS, CB, P);
    FillEmptyPartitions(TS, CB, P);
    GenerateOptimalPartitioning(TS, CB, P);
    *error = PartitioningError(TS, CB, P);
    i++;
    }

  if( Value(ShowProgress) >= 3 )
    {
    d1 = AverageErrorCB(TS, CB, MSE);
    printf("GLA= %9.4f %3i ", PrintableError(d1 - d0, CB), i);
    }
}


/* ----------------------------------------------------------------- */


static void SAStoCodebook(CODEBOOK*     CB,
                          SASchedule*   SAS)
{
  int i;

  if( SASUseToCB(SAS) )
    {
    for( i = 0; i < BookSize(CB); i++ )
      {
      RandomizeVectorBySA(SAS, &Node(CB, i), &Node(CB, i),
                          VectorSize(CB), CB->MaxValue);
      }
    DecreaseTemperature(SAS);
    }
}


/* ================================================================= */


void IMS(TRAININGSET*  TS,
         CODEBOOK*     CB,
         SASchedule*   SAS,
         PARTITIONING* P,
         int*          iteration)
/* Assume: BookSize(TS) >= 1.
           BookSize(CB) >= 1.
           Partition P is current.

  P is not necessary current after execution.
*/
{
  CODEBOOK CB2;
  int      InitialCBSize = BookSize(CB);
  llong    d1, d0;
  int      MinCBSize, MaxCBSize;
  int      StepSize;

  InitializeIMS(TS, CB, P, &MinCBSize, &MaxCBSize, &StepSize, &CB2);

  d0 = MAXLLONG;
  d1 = PartitioningError(TS, &CB2, P);
  if( Value(ShowProgress) >= 2 )
    {
    printf("Iter:%4i %9.4f ", 0, llong2doubleError(TS, d1));
    if( Value(ShowProgress) >= 3 )
      {
      printf(" %9.4f", AverageErrorCB(TS, &CB2, MSE));
      }
    printf("\n");
    fflush(stdout);
    }

  *iteration = 0;
  while( *iteration < Value(Iterations) &&
         (Value(RequireDistortionDecrease) == NO || d1 < d0) )
    {
    if( BookSize(&CB2) == BookSize(CB) )
      {
      if( d1 < d0 ) /* Error is decreased? */
        {
        CopyCodebook(&CB2, CB);
        d0 = d1;
        }
      }

    Iterate(TS, &CB2, SAS, P, MinCBSize, MaxCBSize, StepSize);

    SAStoCodebook(&CB2, SAS);

    if( Value(IMSGLAIterations) > 0 )
      {
      IterateGLA(TS, &CB2, P, &d1);
      }
    else
      {
      d1 = PartitioningError(TS, &CB2, P);
      }

    (*iteration)++;

    ChangeStep(&StepSize);

    if( Value(ShowProgress) >= 2 )
      {
      printf("Iter:%4i %9.4f %4i %4i",
             *iteration, llong2doubleError(TS, d1), BookSize(&CB2), StepSize);
      if( Value(ShowProgress) >= 3 )
        {
        printf(" %9.4f", AverageErrorCB(TS, &CB2, MSE));
        }
      printf("\n");
      fflush(stdout);
      }
    }

  RestoreCodebookSize(TS, &CB2, SAS, P, InitialCBSize);
  d1 = PartitioningError(TS, &CB2, P);
  if( d1 < d0 ) /* Error is decreased? */
    {
    CopyCodebook(&CB2, CB);
    }

  ShutDownIMS(&CB2);

  if( Value(ShowProgress) >= 2 )
    {
    printf("\n");
    }
}


