/*-------------------------------------------------------------------*/
/* GLA.C          Timo Kaukoranta                                    */
/*                                                                   */
/* - Generalized Lloyd Algorithm.                                    */
/*                                                                   */
/*-------------------------------------------------------------------*/

#define ProgName       "GLA"
#define VersionNumber  "V.0.19a"
#define LastUpdated    "28.7.99" /* P.F. */

/* ----------------------------------------------------------------- */

#include <assert.h>
#include <stdio.h>
#include <values.h>

#include "cb.h"
#include "gla.h"
#include "interfc.h"
#include "memctrl.h"
/* #include "search.h" */
#include "sort.h"
#include "sortcb.h"
#include "sa.h"


/* ----------------------------------------------------------------- */
/* To verify the correctness of the search routines,
   change '0' on the following line to '1'. If not, vice versa. */

#if 0

static void VerifyNearest(TRAININGSET*  TS,
                          CODEBOOK*     CB,
                          PARTITIONING* P,
                          int           j,
                          int           nearest)
{
  llong error;
  int   n0 = FindNearestVector(&Node(TS,j),
                               CB,
                               &error,
                               Map(P,j),
                               EUCLIDEANSQ);
  if( nearest != n0 &&
      VectorDist( Vector(TS,j), Vector(CB,n0), VectorSize(CB) ) !=
      VectorDist( Vector(TS,j), Vector(CB,nearest),  VectorSize(CB) ) )
    { printf("ERROR VerifyNearest: j=%4i nearest=%i n0=%i\n", j, nearest, n0);
      printf("d(nearest)=%.0f, d(n0)=%.0f\n",
           (double) VectorDist(Vector(TS,j), Vector(CB,nearest), VectorSize(CB)),
           (double) VectorDist(Vector(TS,j), Vector(CB,n0), VectorSize(CB)) );
    ExitProcessing(-1);
    }
}

#else

#define VerifyNearest(TS, CB, P, j, nearest)

#endif

/* ------------------ Parameters for the GLA module ----------------- */

static struct { int ShowProgress;
                int GLAMethod;
                int UsePDS;
                int PDSInitialGuess;
                int CodevectorCalculation;
              } ModuleParameters;

#define Value(x)      (ModuleParameters.x)

/* These enum definitions are from the *.fac file. */
enum { GLAFullSearch, GLALocal, GLAWu, GLATIE, GLALocalTIE,
       GLAMPS, GLALocalMPS };
enum { PDSoff, PDSon };
enum { FirstVector, CurrentVector };
enum { Centroid, NearestToCentroid };


/* ----------------------------------------------------------------- */


typedef struct
  {
  YESNO* ischanged;      /* Is item 'i' changed */
  int*   set;            /* Indices of changed items */
  int    size;           /* Number of changed items */
  int    maxsize;        /* For future */
  } CHANGESET;


/* ----------------------------------------------------------------- */


typedef struct
  {
/*  CHANGESET* cluster;*/    /* Information about changed clusters */
/*  CHANGESET* tmp;*/        /* Temporary for 'cluster' */
  CHANGESET* codevector; /* Changed code vectors */
  CODEBOOK   prevCB;     /* Codebook of the previous round */
  int*       changedclustersize; /* Sizes of the changed clusters */
  llong*     nearestd;   /* Distance from trainingv to current codev */
  YESNO*     nowcloser;  /* Code vector has come closer to trainingv */
  llong**    DM;         /* Distance matrix for code vectors */
  int**      sDM;        /* Each row gives sorted order of 'DM' */
  int**      scDM;       /* Same as 'sDM' but only for changed clusters */
  int*       codesum;    /* Element sums of the code vectors */
  int*       codeorder;  /* Order for code vectors */
  int*       changedcodeorder;  /* Order for changed code vectors */
/*  llong*  furthestd;*/ /* Distance from furthest trainingv to codev */
  } DATA;


/*-----------------------------  M i s c  ----------------------------*/

#define round(a)      ((a) + 0.5)

llong (*distancef)(VECTORTYPE   v1,
                   VECTORTYPE   v2,
                   int          Vsize,
                   llong        maxdist,
                   DISTANCETYPE disttype);
ERRORFTYPE   errorfunction;
DISTANCETYPE distancefunction;

llong ksum = 0LL;
llong ksumtotal = 0LL;
llong ndistcalcs = 0LL;
llong ndistcalcstotal = 0LL;
int   Nnearer = 0;
int   ndistcalcssavedby0 = 0;

llong* dist2code = NULL;

#define MethodTracksChangedCodeVectors()  ( Value(GLAMethod) == GLALocal || \
                                            Value(GLAMethod) == GLALocalTIE || \
                                            Value(GLAMethod) == GLALocalMPS )


/* =================== Parameters for the GLA module =============== */


void SetGLAParameters(int ShowProgress,
                      int GLAMethod,
                      int UsePDS,
                      int PDSInitialGuess,
                      int CodevectorCalculation,
                      ERRORFTYPE errorf)
{
  Value(ShowProgress)          = ShowProgress;
  Value(GLAMethod)             = GLAMethod;
  Value(UsePDS)                = UsePDS;
  Value(PDSInitialGuess)       = PDSInitialGuess;
  Value(CodevectorCalculation) = CodevectorCalculation;
  errorfunction = errorf;
  distancefunction = DistType(errorf);
}


/* ====================  SIMULATED ANNEALING ======================= */


static void SAStoTS(SASchedule*   SAS,
                    BOOKNODE*     origvector,
                    BOOKNODE*     noisevector,
                    BOOKNODE**    resultvector,
                    int           dim,
                    int           maxvalue)
{
  if( SASUseToTS(SAS) )
    {
    RandomizeVectorBySA(SAS, origvector, noisevector, dim, maxvalue);
    *resultvector = noisevector;
    }
  else
    {
    *resultvector = origvector;
    }
}


/* ======================  CHANGESET  ============================== */

#define CS_size(CS)          ((CS)->size)
#define CS_ischanged(CS,who) ((CS)->ischanged[who])

/* ----------------------------------------------------------------- */

static CHANGESET* CS_make(int size)
{
  CHANGESET* CS = allocate(sizeof(CHANGESET));
  int i;

  CS->ischanged = allocate(size * sizeof(YESNO));
  CS->set       = allocate(size * sizeof(int));
  CS->size      = 0;
  CS->maxsize   = size;

  for( i = 0; i < size; i++ )
    {
    CS->ischanged[i] = NO;
    }
  return( CS );
}


/* ----------------------------------------------------------------- */


static void CS_free(CHANGESET* CS)
{
  deallocate(CS->ischanged);
  deallocate(CS->set);
  deallocate(CS);
}


/* ----------------------------------------------------------------- */


static void CS_clear(CHANGESET* CS)
{
  int i;

  for( i = 0; i < CS_size(CS); i++ )
    {
    CS->ischanged[CS->set[i]] = NO;
    }
  CS->size = 0;
}


/* ----------------------------------------------------------------- */


static void CS_insert(CHANGESET* CS, int who)
{
  if( ! CS_ischanged(CS, who) )
    {
    CS->ischanged[who] = YES;
    CS->set[CS->size] = who;
    CS->size++;
    }
}


/* ----------------------------------------------------------------- */


static int CS_get(CHANGESET* CS, int index)
{
  if( index < CS_size(CS) )
    {
    return( CS->set[index] );
    }
  else
    {
    ErrorMessage("Index (%i) out of CHANGESET size (%i)\n", index, CS_size(CS));
    ExitProcessing( -1 );
    return( -1 );
    }
}


/* ----------------------------------------------------------------- */

#if 0
static void SwapChanged(DATA* D)
{
  CHANGESET* tmp;

  CS_clear(D->cluster);
  tmp = D->cluster;
  D->cluster = D->tmp;
  D->tmp = tmp;
}
#endif


/* ================  VECTOR DISTANCE FUNCTIONS  ==================== */


static llong VectorDistanceTrivial(VECTORTYPE   v1,
                                   VECTORTYPE   v2,
                                   int          Vsize,
                                   llong        maxdist,  /* Foo parameter */
                                   DISTANCETYPE disttype) /* Foo parameter */
/* Precondition: Vsize >= 1 */
{
  int   i = 0;
  llong diff;
  llong dist = 0LL;

  do
    {
    diff = (llong)v1[i] - (llong)v2[i];
    dist += diff * diff;

    ksum++; /* Statistical information */

    } while( ++i < Vsize );

  ndistcalcs++; /* Statistical information */

  return( dist );
}


/* ----------------------------------------------------------------- */


static llong VectorDistancePDS(VECTORTYPE   v1,
                               VECTORTYPE   v2,
                               int          Vsize,
                               llong        maxdist,
                               DISTANCETYPE disttype) /* Foo parameter */
/* Precondition: Vsize >= 1 */
/* if VectorDistance > maxdist, returns at least maxdist */
{
  int   i = 0;
  llong diff;
  llong dist = 0LL;

  do
    {
    diff = (llong)v1[i] - (llong)v2[i];
    dist += diff * diff;

    ksum++; /* Statistical information */

    } while( dist < maxdist && ++i < Vsize );

  ndistcalcs++; /* Statistical information */

  return( dist );
}


/* =====================  DISTANCE MATRIX  ========================= */


static void DistanceMatrixMake(CODEBOOK* CB, llong*** D)
{
  int i;
  
  *D = (llong**)allocate(BookSize(CB) * sizeof(llong*));
  for( i = 0; i < BookSize(CB); i++ )
    {
    /* Sentinel is at the end of the row. */
    (*D)[i] = (llong*)allocate((BookSize(CB)+1) * sizeof(llong));
    }
}


/* ----------------------------------------------------------------- */


static void DistanceMatrixFree(CODEBOOK* CB, llong*** D)
{
  int i;

  for( i = 0; i < BookSize(CB); i++ )
    {
    deallocate((*D)[i]);
    }
  deallocate(*D);
}


/* ----------------------------------------------------------------- */


static void DistanceMatrixCalculate(CODEBOOK*    CB,
                                    llong*       D[],
                                    DISTANCETYPE disttype)
{
  int   i, j;
  llong dist;

  for( i = 0; i < BookSize(CB); i++ )
    {
    D[i][i] = 0LL;
    D[i][BookSize(CB)] = MAXLLONG;  /* Sentinel */
    for( j = i+1; j < BookSize(CB); j++ )
      {
      dist = distancef(Vector(CB, i), Vector(CB, j), VectorSize(CB),
                       MAXLLONG, disttype);
      D[i][j] = D[j][i] = dist;
      }
    }
}


/* ----------------------------------------------------------------- */


static void DistanceMatrixCalculateChanged(CODEBOOK*    CB,
                                           llong*       D[],
                                           CHANGESET*   CS,
                                           DISTANCETYPE disttype)
{
  int   i, j, l, h;
  llong dist;

  /* For each changed code vector i. */
  for( l = 0; l < CS_size(CS); l++ )
    {
    i = CS_get(CS, l);
    /* For each unchanged code vector j before first changed code vector. */
    for( j = 0;  j < CS_get(CS, 0); j++ )
      {
      dist = distancef(Vector(CB, i), Vector(CB, j), VectorSize(CB),
                       MAXLLONG, disttype);
      D[i][j] = D[j][i] = dist;
      }

    /* For each unchanged code vector j,
       i.e. code vectors between changed code vectors. */
    for( h = 0; h < CS_size(CS)-1; h++ )
      {
      for( j = CS_get(CS, h)+1;  j < CS_get(CS, h+1); j++ )
        {
        dist = distancef(Vector(CB, i), Vector(CB, j), VectorSize(CB),
                         MAXLLONG, disttype);
        D[i][j] = D[j][i] = dist;
        }
      }

    /* For each unchanged code vector j after last changed code vector. */
    for( j = CS_get(CS, CS_size(CS)-1)+1;  j < BookSize(CB); j++ )
      {
      dist = distancef(Vector(CB, i), Vector(CB, j), VectorSize(CB),
                       MAXLLONG, disttype);
      D[i][j] = D[j][i] = dist;
      }

    /* For each changed code vector j. */
    for( h = l+1; h < CS_size(CS); h++ )
      {
      j = CS_get(CS, h);
      dist = distancef(Vector(CB, i), Vector(CB, j), VectorSize(CB),
                       MAXLLONG, disttype);
      D[i][j] = D[j][i] = dist;
      }
    }

#if 0
  for( l = 0; l < CS_size(CS); l++ )
    {
    i = CS_get(CS, l);
    for( j = 0; j < BookSize(CB); j++ )
      {
      dist = distancef(Vector(CB, i), Vector(CB, j), VectorSize(CB),
                       MAXLLONG, disttype);
      D[i][j] = D[j][i] = dist;
      }
    }
#endif
}


/* ======================  ORDER MATRIX  =========================== */


static void OrderMatrixMake(CODEBOOK* CB, int*** L)
{
  int i;
  
  *L = (int**) allocate(BookSize(CB) * sizeof(int*));
  for( i = 0; i < BookSize(CB); i++ )
    {
    /* Sentinel is at the end of the row. */
    (*L)[i] = (int*) allocate((BookSize(CB)+1) * sizeof(int));
    }
}


/* ----------------------------------------------------------------- */


static void OrderMatrixFree(CODEBOOK* CB, int*** L)
{
  int i;

  for( i = 0; i < BookSize(CB); i++ )
    {
    deallocate((*L)[i]);
    }
  deallocate(*L);
}


/* ----------------------------------------------------------------- */


static int cmpDistanceMatrix(const void* e1, const void* e2, const void* info)
{
  return( ((llong*)info)[*((int*)e1)] < ((llong*)info)[*((int*)e2)] ? 1 : 0 );
}


/* ----------------------------------------------------------------- */


static void OrderMatrixInitialize(CODEBOOK* CB, int* L[])
{
  int   i, j;

  for( i = 0; i < BookSize(CB); i++ )
    {
    /* Sentinel is initialized also. */
    for( j = 0; j < BookSize(CB)+1; j++ )
      {
      L[i][j] = j;
      }
    }
}


/* ----------------------------------------------------------------- */


static void OrderMatrixQuickSort(CODEBOOK*    CB,
                                 llong*       D[],
                                 int*         L[],
                                 DISTANCETYPE disttype)
{
  int i;

  for( i = 0; i < BookSize(CB); i++ )
    {
    /* Sentinel is not sorted. */
    QuickSort(L[i], BookSize(CB), sizeof(int), D[i], cmpDistanceMatrix);
    }
}


/* ----------------------------------------------------------------- */


static void OrderMatrixInsertSort(CODEBOOK*    CB,
                                  llong*       D[],
                                  int*         L[],
                                  DISTANCETYPE disttype)
{
  int i;

  for( i = 0; i < BookSize(CB); i++ )
    {
    /* Sentinel is not sorted. */
    InsertSort(L[i], BookSize(CB), sizeof(int), D[i], cmpDistanceMatrix);
    }
}


/* ----------------------------------------------------------------- */


static void OrderMatrixPickChanged(CODEBOOK*  CB,
                                   CHANGESET* CS,
                                   int*       sDM[],
                                   int*       scDM[])
{
  int  i, j, l;
  int* sDMi;
  int* scDMi;

  for( i = 0; i < BookSize(CB); i++ )
    {
    sDMi = sDM[i];
    scDMi = scDM[i];
    for( j = 0, l = 0; l < CS_size(CS); j++ )
      {
      if( CS_ischanged(CS, sDMi[j]) )
        {
        scDMi[l++] = sDMi[j];
        }
      }
    /* Reference to the sentinel. */
    scDMi[l] = BookSize(CB);
    }
}


/* ================================================================= */

#if 0
static void CollectChangedClusters(CODEBOOK*     CB,
                                   PARTITIONING* P,
                                   DATA*         D)
{
  int i;
  
  D->nchangedset = 0;
  for( i = 0; i < BookSize(CB); i++ )
    {
    if( D->ischanged[i] == YES )
      {
      D->changedset[D->nchangedset] = i;
      D->changedsizeset[D->nchangedset] = UniqueVectors(P,i);
      D->nchangedset++;
      }
    }
}
#endif

/* ----------------------------------------------------------------- */


static llong CurrentDistortion(TRAININGSET* TS, DATA* D)
{
  int    j;
  llong  dist = 0LL;
  llong* distance = D->nearestd;

  for( j = 0; j < BookSize(TS); j++ )
    {
    dist += (llong)distance[j] * (llong)VectorFreq(TS, j);
    }
  return( dist );
}


/* ----------------------------------------------------------------- */


static double Total2AverageDistortion(TRAININGSET* TS, llong totaldist)
{
  return( (double) totaldist / (double) (TotalFreq(TS) * VectorSize(TS)) );
}


/* ----------------------------------------------------------------- */


static void CheckChangedCodeVectors(CODEBOOK*     CB,
                                    PARTITIONING* P,
                                    DATA*         D)
{
  int i;
  
  CS_clear(D->codevector);

  for( i = 0; i < BookSize(CB); i++ )
    {
    if( !(EqualVectors(Vector(CB,i), Vector(&D->prevCB,i), VectorSize(CB))) )
      {
      CS_insert(D->codevector, i);
      D->changedclustersize[i] = CCFreq(P, i);
      }
    }
  CopyCodebook(CB, &D->prevCB);
}


/* ==========================  PRINT  ============================== */


static void PrintIterationInfo(TRAININGSET*  TS,
                               CODEBOOK*     CB,
                               SASchedule*   SAS,
                               PARTITIONING* P,
                               int           Iterations,
                               llong         d1,
                               llong         d2,
                               int           changedmappings,
                               ERRORFTYPE    errorfunction,
                               DATA*         D)
{
  int    i;
  int    Nc = 0;
  double currD = Total2AverageDistortion(TS, d2);
  double deltaD = Total2AverageDistortion(TS, d1-d2);

  printf("Iter= %3i ERR= %9.4f dERR= %9.4f ",
         Iterations,
         PrintableError(currD, CB),
         d1 == MAXLLONG ? 9999.9999 : PrintableError(deltaD, CB));

  if( SASInUse(SAS) )
    {
    printf(" temper.= %8.4f ", SAS->CurrentTemperature);
    }

  if( Value(ShowProgress) >= 3 )
    {
    for( i = 0; i < CS_size(D->codevector); i++ )
      {
      Nc += D->changedclustersize[CS_get(D->codevector, i)];
      }
    printf("Mc/Nc/Nn/Nr= %3i %4i %4i %4i ",
           CS_size(D->codevector), Nc, Nnearer, changedmappings);

    /* Statistical information */
    ksumtotal += ksum;
    ndistcalcstotal += ndistcalcs;

    if( ndistcalcs > 0LL && ndistcalcstotal > 0LL )
      {
      printf("#d= %7.0f %7.0f %5.2f %7.0f %7.0f %5.2f ",
             (double)ndistcalcs,
             (double)ksum,
             (double)ksum/(double)ndistcalcs,
             (double)ndistcalcstotal,
             (double)ksumtotal,
             (double)ksumtotal/(double)ndistcalcstotal);
      }
    printf("%7i ", ndistcalcssavedby0);
    }

  ksum = 0LL;
  ndistcalcs = 0LL;
  ndistcalcssavedby0 = 0;

  printf("\n");
  fflush(stdout);
}


/*=========================  NEAREST SEARCHING  ============================*/


static void CheckNearest(BOOKNODE*     v,
                         CODEBOOK*     CB,
                         int           candidate,
                         int*          minindex,
                         llong*        minerror,
                         DISTANCETYPE  disttype)
{
  llong e = distancef(v->vector,
                      Vector(CB, candidate),
                      VectorSize(CB),
                      *minerror,
                      disttype);

  if( e < *minerror )
    {
    *minerror = e;
    *minindex = candidate;
    }
}


/* ----------------------------------------------------------------- */


static int FindNearestVectorInCB(BOOKNODE*     v,
                                 int           tvector,
                                 CODEBOOK*     CB,
                                 int           current,
                                 DATA*         D,
                                 DISTANCETYPE  disttype)
{
  int    i;
  int    minindex;
  llong* minerror = &D->nearestd[tvector];

  if( Value(PDSInitialGuess) == CurrentVector )
    {
    minindex  = current;
    for( i = 0; i < current && *minerror != 0; i++ )
      {
      CheckNearest(v, CB, i, &minindex, minerror, disttype);
      }
    for( i++; i < BookSize(CB) && *minerror != 0; i++ )
      {
      CheckNearest(v, CB, i, &minindex, minerror, disttype);
      }
    }
  else
    {
    minindex = 0;
    *minerror = MAXLLONG;
    for( i = 0; i < BookSize(CB) && *minerror != 0; i++ )
      {
      CheckNearest(v, CB, i, &minindex, minerror, disttype);
      }
    }

  ndistcalcssavedby0 += (BookSize(CB) - i);

  return( minindex );
}

/* ----------------------------------------------------------------- */


static int FindNearestVectorInSet(BOOKNODE*    v,
                                  int          tvector,
                                  CODEBOOK*    CB,
                                  int          current,
                                  DATA*        D,
                                  DISTANCETYPE disttype)
{
  int    i;
  int    minindex = current;
  llong* minerror = &D->nearestd[tvector];

  for( i = 0; i < CS_size(D->codevector) && *minerror != 0; i++ )
    {
    CheckNearest(v, CB, CS_get(D->codevector, i), &minindex, minerror, disttype);
    }

  ndistcalcssavedby0 += (CS_size(D->codevector) - i);

  return( minindex );
}


/* ----------------------------------------------------------------- */


static void CalculateCurrentDistances(TRAININGSET*  TS,
                                      CODEBOOK*     CB,
                                      PARTITIONING* P,
                                      DATA*         D,
                                      DISTANCETYPE  disttype)
{
  int   i, j;
  llong olddist;

  Nnearer = 0;
  if( Value(PDSInitialGuess) == CurrentVector )
    {
    if( MethodTracksChangedCodeVectors() )
      {
      for( i = 0; i < BookSize(CB); i++ )
        {
        if( CS_ischanged(D->codevector, i) )
          {
          for( j = FirstVector(P, i); ! EndOfPartition(j); j = NextVector(P, j) )
            {
            olddist = D->nearestd[j];
            D->nearestd[j] = distancef(Vector(TS, j),
                                       Vector(CB, i),
                                       VectorSize(CB),
                                       MAXLLONG,
                                       disttype);
            if( D->nearestd[j] < olddist )
              {
              D->nowcloser[j] = YES;
              Nnearer += VectorFreq(TS, j);
              }
            else
              {
              D->nowcloser[j] = NO;
              }
            }
          }
        else
          {
          for( j = FirstVector(P, i); ! EndOfPartition(j); j = NextVector(P, j) )
            {
            D->nowcloser[j] = NO;
            }
          }
        }
      }
    else
      {
      for( j = 0; j < BookSize(TS); j++ )
        {
        D->nearestd[j] = distancef(Vector(TS, j),
                                   Vector(CB, Map(P,j)),
                                   VectorSize(CB),
                                   MAXLLONG,
                                   disttype);
        }
      }
    }
}


/*=========================  CODEBOOK HANDLING  ============================*/


static void FindNearestVectorsGeneral(TRAININGSET*  TS,
                                      CODEBOOK*     CB,
                                      SASchedule*   SAS, /* NULL allowed */
                                      PARTITIONING* P,
                                      DATA*         D,
                                      int*          changedmappings)
{
  int       j;
  int       nearest;
  BOOKNODE  noisevector;
  BOOKNODE* tsvector;

  *changedmappings = 0;
  noisevector = CreateEmptyNode(VectorSize(CB));

  CalculateCurrentDistances(TS, CB, P, D, distancefunction);

  /* Find mapping from training vector to code vector */
  for(j = 0; j < BookSize(TS); j++)
    {
    SAStoTS(SAS, &Node(TS, j), &noisevector, &tsvector,
            VectorSize(TS), TS->MaxValue);

    nearest = FindNearestVectorInCB(tsvector, j, CB, Map(P, j), D,
                                    distancefunction);
    VerifyNearest(TS, CB, P, j, nearest);

    if( nearest != Map(P, j) )
      {
      ChangePartition(TS, P, nearest, j);
      (*changedmappings)++;
      }
    }

  FreeNode(noisevector);
}


/*==================  GENERALIZED LLOYD ALGORITHM  =====================*/


static void FindNearestVectorsLocal(TRAININGSET*  TS,
                                    CODEBOOK*     CB,
                                    SASchedule*   SAS, /* NULL allowed */
                                    PARTITIONING* P,
                                    DATA*         D,
                                    int*          changedmappings)
{
  int       j;
  int       nearest;
  BOOKNODE  noisevector;
  BOOKNODE* tsvector;

  *changedmappings = 0;
  noisevector = CreateEmptyNode(VectorSize(CB));

  CalculateCurrentDistances(TS, CB, P, D, distancefunction);

  /* Find mapping from training vector to code vector */
  for( j = 0; j < BookSize(TS); j++ )
    {
    SAStoTS(SAS, &Node(TS, j), &noisevector, &tsvector,
            VectorSize(TS), TS->MaxValue);

    if( CS_ischanged(D->codevector, Map(P, j)) && ! D->nowcloser[j] )
      {
      nearest = FindNearestVectorInCB(tsvector, j, CB, Map(P, j), D,
                                      distancefunction);
      }
    else
      {
      nearest = FindNearestVectorInSet(tsvector, j, CB, Map(P, j), D,
                                       distancefunction);
      }
    VerifyNearest(TS, CB, P, j, nearest);

    if( nearest != Map(P, j) )
      {
      ChangePartition(TS, P, nearest, j);
      (*changedmappings)++;
      }
    }

  FreeNode(noisevector);
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

  for( j = FirstVector(P, Pindex); ! EndOfPartition(j); j = NextVector(P, j) )
    {
    error = distancef(Vector(TS, j), V->vector, VectorSize(TS),
                      MinError, distancefunction);
    if( error < MinError )
      {
      MinError  = error;
      *NearestV = j;
      }
    }
}


/* ----------------------------------------------------------------- */


static void ComputeNewBook(TRAININGSET*  TS,
                           CODEBOOK*     CB,
                           SASchedule*   SAS,  /* NULL allowed */
                           PARTITIONING* P,
                           DATA*         D)
{
  int i;
  int NearCent;

  GenerateOptimalCodebook(TS, CB, P);
  FillEmptyPartitions(TS, CB, P);

  switch( Value(CodevectorCalculation) )
    {
    case Centroid:
      {
      break;
      }
    case NearestToCentroid:
      {
      for( i = 0; i < BookSize(CB); i++ )
        {
        NearestVectorInPartition(TS, P, i, &Node(CB, i), &NearCent);
        CopyNode(&Node(TS, NearCent), &Node(CB, i), VectorSize(CB));
        }
      break;
      }
    default:
      {
      break;
      }
    }
  if( SASUseToCB(SAS) )
    {
    for( i = 0; i < BookSize(CB); i++ )
      {
      RandomizeVectorBySA(SAS, &Node(CB, i), &Node(CB, i),
                          VectorSize(CB), CB->MaxValue);
      }
    }
  if( MethodTracksChangedCodeVectors() ||
      Value(ShowProgress) >= 3 )
    {
    CheckChangedCodeVectors(CB, P, D);
    }
}


/* ====================  FULL SEARCH  ============================== */


/* Lloyd-iteration:
   1. PARTITIONING: Map training vectors (in TS) to codevectors (in CB).
   2. CENTROIDS: Calculate centroid for each partition. */

static void RunGLAFullSearch(TRAININGSET*  TS,
                             CODEBOOK*     CB,
                             SASchedule*   SAS,
                             PARTITIONING* P,
                             DATA*         D,
                             int           IterLimit,
                             double*       InitError,
                             double*       FinalError,
                             int*          Iterations)
{
  llong d0, d1 = MAXLLONG, d2;
  int   changedmappings;

  *Iterations = 0;
  FindNearestVectorsGeneral(TS, CB, SAS, P, D, &changedmappings);
  d0 = d2 = CurrentDistortion(TS, D);

  if( Value(ShowProgress) >= 2 )
    {
    PrintIterationInfo(TS, CB, SAS, P, *Iterations, d1, d2, changedmappings,
                       errorfunction, D);
    }

  while( ( d2 < d1 && *Iterations < IterLimit )
         ||
         ( SASInUse(SAS) && SASEffective(SAS, VectorSize(CB), CB->MaxValue) ) )
    {
    ComputeNewBook(TS, CB, SAS, P, D);
    FindNearestVectorsGeneral(TS, CB, SAS, P, D, &changedmappings);
    d1 = d2;
    d2 = CurrentDistortion(TS, D);

    DecreaseTemperature(SAS);
    (*Iterations)++;

    if( Value(ShowProgress) >= 2 )
      {
      PrintIterationInfo(TS, CB, SAS, P, *Iterations, d1, d2, changedmappings,
                         errorfunction, D);
      }
    }

  *InitError  = Total2AverageDistortion(TS, d0);
  *FinalError = Total2AverageDistortion(TS, d2);
}


/* =======================  LOCAL  ================================= */


/* Lloyd-iteration:
   1. PARTITIONING: Map training vectors (in TS) to codevectors (in CB).
   2. CENTROIDS: Calculate centroid for each partition. */

static void RunGLALocal(TRAININGSET*  TS,
                        CODEBOOK*     CB,
                        SASchedule*   SAS,
                        PARTITIONING* P,
                        DATA*         D,
                        int           IterLimit,
                        double*       InitError,
                        double*       FinalError,
                        int*          Iterations)
{
  llong d0, d1 = MAXLLONG, d2;
  int   changedmappings;

  *Iterations = 0;
  FindNearestVectorsGeneral(TS, CB, SAS, P, D, &changedmappings);
  d0 = d2 = CurrentDistortion(TS, D);

  if( Value(ShowProgress) >= 2 )
    {
    PrintIterationInfo(TS, CB, SAS, P, *Iterations, d1, d2, changedmappings,
                       errorfunction, D);
    }

  while( ( d2 < d1 && *Iterations < IterLimit )
         ||
         ( SASInUse(SAS) && SASEffective(SAS, VectorSize(CB), CB->MaxValue) ) )
    {
    ComputeNewBook(TS, CB, SAS, P, D);
    FindNearestVectorsLocal(TS, CB, SAS, P, D, &changedmappings);
    d1 = d2;
    d2 = CurrentDistortion(TS, D);

    DecreaseTemperature(SAS);
    (*Iterations)++;

    if( Value(ShowProgress) >= 2 )
      {
      PrintIterationInfo(TS, CB, SAS, P, *Iterations, d1, d2, changedmappings,
                         errorfunction, D);
      }
    }

  *InitError  = Total2AverageDistortion(TS, d0);
  *FinalError = Total2AverageDistortion(TS, d2);
}


/* ======================  GLA WU  ================================= */


static int B(llong Dk[], int Lk[], int start, int end, llong dist)
{
  int half = (start + end) / 2;

  if( dist == 0LL )
    {
    return( 0 );
    }

  dist *= 4;

  assert( Dk[Lk[start]] < dist );
  assert( dist < Dk[Lk[end]] );

  while( start < end - 1 )
    {
    if( Dk[Lk[half]] < dist )
      {
      start = half;
      }
    else /* dist <= Dk[Lk[half]] */
      {
      end = half;
      }
    half = (start + end) / 2;
    }

  assert( Dk[Lk[start]] < dist );
  assert( dist <= Dk[Lk[start+1]] );

  return( start );
}


/* ----------------------------------------------------------------- */


static int FindNearestVectorGLAWu(BOOKNODE*    v,
                                  int          tvector,
                                  CODEBOOK*    CB,
                                  int          t,
                                  DATA*        D,
                                  DISTANCETYPE disttype,
                                  int          tmp[],
                                  int*         label)
{
  int   i;
  int   Bt, Bk;
  int*  S = (int*)allocate(BookSize(CB) * sizeof(int));
  int   Ssize = 0, Soldsize;
  int   k;
  llong distk;
  int*  Lk;
  int   startlabel = *label;
  llong* distt = &D->nearestd[tvector];

  *distt = distancef(v->vector, Vector(CB, t), VectorSize(CB),  /* ????? */
                       MAXLLONG, disttype);
  Bt = B(D->DM[t], D->sDM[t], 0, BookSize(CB), *distt);

  Ssize = Bt;
  for( i = 1; i <= Bt; i++ )
    {
    S[i-1] = D->sDM[t][i];
    }

  if( Value(ShowProgress) >= 4 )
    {
    printf("Bt/t= %3i %3i ", Bt, t);
    }

  while( Ssize > 0 )
    {
#if 0
    int quess = Ssize/2;
    k = S[quess];
#endif
    /* Pick from 'S' */
    k = S[0];
    Lk = D->sDM[k];
    distk = distancef(v->vector, Vector(CB, k), VectorSize(CB),
                      MAXLLONG, disttype);
    Bk = B(D->DM[k], Lk, 0, BookSize(CB), distk);

    if( Value(ShowProgress) >= 5 )
      {
      printf("Ssize/t/Bk/k= %3i %3i %3i %3i\n", Ssize, t , Bk, k);
      }

    (*label)++;
    if( distk < *distt )
      {
      t = k;
      *distt = distk;

      /* Intersection S = S \ Bk(v) */
       for( i = 1; i < Ssize; i++ )
        {
        tmp[S[i]] = *label;
        }
#if 0
      for( i = 0; i < quess; i++ ) tmp[S[i]] = *label;
      for( i++; i < Ssize; i++ ) tmp[S[i]] = *label;
#endif
      Ssize = 0;
      for( i = 1; i <= Bk; i++ )
        {
        if( tmp[Lk[i]] == *label )
          {
          S[Ssize++] = Lk[i];
          }
        }
      }
    else
      {
      /* Intersection S = S \ Bk(v) */
      for( i = 1; i <= Bk; i++ )
        {
        tmp[Lk[i]] = *label;
        }

      Soldsize = Ssize;
      Ssize = 0;
      for( i = 1; i < Soldsize; i++ )
        {
        if( tmp[S[i]] == *label )
          {
          S[Ssize++] = S[i];
          }
        }
#if 0
      for( i = 0; i < quess; i++ ) if( tmp[S[i]] == *label ) S[Ssize++] = S[i];
      for( i++; i < Soldsize; i++ ) if( tmp[S[i]] == *label ) S[Ssize++] = S[i];
#endif
      }
    }

  deallocate(S);

  if( Value(ShowProgress) >= 4 )
    {
    printf("Intersecs= %3i\n", *label-startlabel);
    }

  return( t );
}


/* ----------------------------------------------------------------- */


static void FNVGLAWu(TRAININGSET*  TS,
                                    CODEBOOK*     CB,
                                    SASchedule*   SAS, /* NULL allowed */
                                    PARTITIONING* P,
                                    DATA*         D,
                                    int*          changedmappings)
{
  int       i, j;
  int       nearest;
  BOOKNODE  noisevector;
  BOOKNODE* tsvector;
  int*      tmp = allocate(BookSize(CB) * sizeof(int));
  int       label = 0;

  noisevector = CreateEmptyNode(VectorSize(CB));
  for( i = 0; i < BookSize(CB); i++ )
    {
    tmp[i] = label;
    }

  *changedmappings= 0;
  /* Find mapping from training vector to code vector */
  for(j = 0; j < BookSize(TS); j++)
    {
    SAStoTS(SAS, &Node(TS, j), &noisevector, &tsvector,
            VectorSize(TS), TS->MaxValue);

    nearest = FindNearestVectorGLAWu(tsvector, j, CB, Map(P, j), D,
                                     distancefunction, tmp, &label);
    VerifyNearest(TS, CB, P, j, nearest);

    if( nearest != Map(P, j) )
      {
      ChangePartition(TS, P, nearest, j);
      (*changedmappings)++;
      }
    }

  deallocate(tmp);
  FreeNode(noisevector);
}


/* ----------------------------------------------------------------- */


static void RunGLAWu(TRAININGSET*  TS,
                     CODEBOOK*     CB,
                     SASchedule*   SAS,
                     PARTITIONING* P,
                     DATA*         D,
                     int           IterLimit,
                     double*       InitError,
                     double*       FinalError,
                     int*          Iterations)
{
  llong d0, d1 = MAXLLONG, d2;
  int   changedmappings;

  *Iterations = 0;

  DistanceMatrixCalculate(CB, D->DM, distancefunction);
  OrderMatrixInitialize(CB, D->sDM);
  OrderMatrixQuickSort(CB, D->DM, D->sDM, distancefunction);
  FNVGLAWu(TS, CB, SAS, P, D, &changedmappings);
  d0 = d2 = CurrentDistortion(TS, D);

  if( Value(ShowProgress) >= 2 )
    {
    PrintIterationInfo(TS, CB, SAS, P, *Iterations, d1, d2, changedmappings,
                       errorfunction, D);
    }

  while( ( d2 < d1 && *Iterations < IterLimit )
         ||
         ( SASInUse(SAS) && SASEffective(SAS, VectorSize(CB), CB->MaxValue) ) )
    {
    ComputeNewBook(TS, CB, SAS, P, D);
    DistanceMatrixCalculate(CB, D->DM, distancefunction);
    OrderMatrixQuickSort(CB, D->DM, D->sDM, distancefunction);
    FNVGLAWu(TS, CB, SAS, P, D, &changedmappings);
    d1 = d2;
    d2 = CurrentDistortion(TS, D);

    DecreaseTemperature(SAS);
    (*Iterations)++;

    if( Value(ShowProgress) >= 2 )
      {
      PrintIterationInfo(TS, CB, SAS, P, *Iterations, d1, d2, changedmappings,
                         errorfunction, D);
      }
    }

  *InitError  = Total2AverageDistortion(TS, d0);
  *FinalError = Total2AverageDistortion(TS, d2);
}


/* ========================  TIE  ================================== */


static int FindNearestVectorGLATIE(BOOKNODE*    v,
                                   int          tvector,
                                   CODEBOOK*    CB,
                                   int          current,
                                   DATA*        D,
                                   DISTANCETYPE disttype)
{
  int    i;
  int    minindex;
  llong  radius;
  llong* minerror = &D->nearestd[tvector];
  llong* DMc  = D->DM[current];
  int*   sDMc = D->sDM[current];

  minindex = current;
  radius = 4 * (*minerror);

  /* Utilize the sentinel and note that Lc[0]==current */
  for( i = 1; DMc[sDMc[i]] < radius && *minerror != 0; i++ )
    {
    CheckNearest(v, CB, sDMc[i], &minindex, minerror, disttype);
    }

  return( minindex );
}


/* ----------------------------------------------------------------- */


static void FindNearestVectorsGLATIE(TRAININGSET*  TS,
                                     CODEBOOK*     CB,
                                     SASchedule*   SAS, /* NULL allowed */
                                     PARTITIONING* P,
                                     DATA*         D,
                                     int*          changedmappings)
{
  int       j;
  int       nearest;
  BOOKNODE  noisevector;
  BOOKNODE* tsvector;

  noisevector = CreateEmptyNode(VectorSize(CB));
  CalculateCurrentDistances(TS, CB, P, D, distancefunction);

  *changedmappings = 0;
  /* Find mapping from training vector to code vector */
  for( j = 0; j < BookSize(TS); j++ )
    {
    SAStoTS(SAS, &Node(TS, j), &noisevector, &tsvector,
            VectorSize(TS), TS->MaxValue);

    nearest = FindNearestVectorGLATIE(tsvector, j, CB, Map(P, j), D,
                                      distancefunction);
    VerifyNearest(TS, CB, P, j, nearest);

    if( nearest != Map(P, j) )
      {
      ChangePartition(TS, P, nearest, j);
      (*changedmappings)++;
      }
    }

  FreeNode(noisevector);
}


/* ----------------------------------------------------------------- */


static void InitializeGLATIE(CODEBOOK* CB, llong*** D, int*** L)
{
  DistanceMatrixMake(CB, D);
  OrderMatrixMake(CB, L);
}


/* ----------------------------------------------------------------- */


static void ShutDownGLATIE(CODEBOOK* CB, llong*** D, int*** L)
{
  DistanceMatrixFree(CB, D);
  OrderMatrixFree(CB, L);
}


/* ----------------------------------------------------------------- */


static void RunGLATIE(TRAININGSET*  TS,
                      CODEBOOK*     CB,
                      SASchedule*   SAS,
                      PARTITIONING* P,
                      DATA*         D,
                      int           IterLimit,
                      double*       InitError,
                      double*       FinalError,
                      int*          Iterations)
{
  llong d0, d1 = MAXLLONG, d2;
  int   changedmappings;

  *Iterations = 0;

  DistanceMatrixCalculate(CB, D->DM, distancefunction);
  OrderMatrixInitialize(CB, D->sDM);
  OrderMatrixQuickSort(CB, D->DM, D->sDM, distancefunction);
  FindNearestVectorsGLATIE(TS, CB, SAS, P, D, &changedmappings);
  d0 = d2 = CurrentDistortion(TS, D);

  if( Value(ShowProgress) >= 2 )
    {
    PrintIterationInfo(TS, CB, SAS, P, *Iterations, d1, d2, changedmappings,
                       errorfunction, D);
    }

  while( ( d2 < d1 && *Iterations < IterLimit )
         ||
         ( SASInUse(SAS) && SASEffective(SAS, VectorSize(CB), CB->MaxValue) ) )
    {
    ComputeNewBook(TS, CB, SAS, P, D);
    DistanceMatrixCalculate(CB, D->DM, distancefunction);
    OrderMatrixQuickSort(CB, D->DM, D->sDM, distancefunction);
    FindNearestVectorsGLATIE(TS, CB, SAS, P, D, &changedmappings);
    d1 = d2;
    d2 = CurrentDistortion(TS, D);

    DecreaseTemperature(SAS);
    (*Iterations)++;

    if( Value(ShowProgress) >= 2 )
      {
      PrintIterationInfo(TS, CB, SAS, P, *Iterations, d1, d2, changedmappings,
                         errorfunction, D);
      }
    }

  *InitError  = Total2AverageDistortion(TS, d0);
  *FinalError = Total2AverageDistortion(TS, d2);
}


/* =====================  LOCAL TIE  =============================== */


static int FindNearestVectorGLALocalTIE(BOOKNODE*    v,
                                        int          tvector,
                                        CODEBOOK*    CB,
                                        int          current,
                                        DATA*        D,
                                        DISTANCETYPE disttype)
{
  int    i;
  int    minindex;
  llong  radius;
  llong* minerror = &D->nearestd[tvector];
  llong* DMc = D->DM[current];
  int*   scDMc = D->scDM[current];

  minindex = current;
  radius = 4 * (*minerror);

  /* Utilize the sentinel and note that Lc[0]==current */
  for( i = 0; DMc[scDMc[i]] < radius && *minerror != 0; i++ )
    {
    CheckNearest(v, CB, scDMc[i], &minindex, minerror, disttype);
    }

  return( minindex );
}


/* ----------------------------------------------------------------- */


static void FindNearestVectorsGLALocalTIE(TRAININGSET*  TS,
                                          CODEBOOK*     CB,
                                          SASchedule*   SAS, /* NULL allowed */
                                          PARTITIONING* P,
                                          DATA*         D,
                                          int*          changedmappings)
{
  int       j;
  int       nearest;
  BOOKNODE  noisevector;
  BOOKNODE* tsvector;

  noisevector = CreateEmptyNode(VectorSize(CB));
  CalculateCurrentDistances(TS, CB, P, D, distancefunction);

  *changedmappings = 0;
  /* Find mapping from training vector to code vector */
  for( j = 0; j < BookSize(TS); j++ )
    {
    SAStoTS(SAS, &Node(TS, j), &noisevector, &tsvector,
            VectorSize(TS), TS->MaxValue);

    if( CS_ischanged(D->codevector, Map(P, j)) && ! D->nowcloser[j] )
      {
      nearest = FindNearestVectorGLATIE(tsvector, j, CB, Map(P, j), D,
                                        distancefunction);
      }
    else
      {
      nearest = FindNearestVectorGLALocalTIE(tsvector, j, CB, Map(P, j), D,
                                             distancefunction);
      }
    VerifyNearest(TS, CB, P, j, nearest);

    if( nearest != Map(P, j) )
      {
      ChangePartition(TS, P, nearest, j);
      (*changedmappings)++;
      }
    }

  FreeNode(noisevector);
}


/* ----------------------------------------------------------------- */


static void RunGLALocalTIE(TRAININGSET*  TS,
                           CODEBOOK*     CB,
                           SASchedule*   SAS,
                           PARTITIONING* P,
                           DATA*         D,
                           int           IterLimit,
                           double*       InitError,
                           double*       FinalError,
                           int*          Iterations)
{
  llong d0, d1 = MAXLLONG, d2;
  int   changedmappings;

  *Iterations = 0;

  DistanceMatrixCalculate(CB, D->DM, distancefunction);
  OrderMatrixInitialize(CB, D->sDM);
  OrderMatrixQuickSort(CB, D->DM, D->sDM, distancefunction);
  FindNearestVectorsGLATIE(TS, CB, SAS, P, D, &changedmappings);
  d0 = d2 = CurrentDistortion(TS, D);

  if( Value(ShowProgress) >= 2 )
    {
    PrintIterationInfo(TS, CB, SAS, P, *Iterations, d1, d2, changedmappings,
                       errorfunction, D);
    }

  while( ( d2 < d1 && *Iterations < IterLimit )
         ||
         ( SASInUse(SAS) && SASEffective(SAS, VectorSize(CB), CB->MaxValue) ) )
    {
    ComputeNewBook(TS, CB, SAS, P, D);
    DistanceMatrixCalculateChanged(CB, D->DM, D->codevector, distancefunction);
    OrderMatrixQuickSort(CB, D->DM, D->sDM, distancefunction);
    OrderMatrixPickChanged(CB, D->codevector, D->sDM, D->scDM);
    FindNearestVectorsGLALocalTIE(TS, CB, SAS, P, D, &changedmappings);
    d1 = d2;
    d2 = CurrentDistortion(TS, D);

    DecreaseTemperature(SAS);
    (*Iterations)++;

    if( Value(ShowProgress) >= 2 )
      {
      PrintIterationInfo(TS, CB, SAS, P, *Iterations, d1, d2, changedmappings,
                         errorfunction, D);
      }
    }

  *InitError  = Total2AverageDistortion(TS, d0);
  *FinalError = Total2AverageDistortion(TS, d2);
}


/* ========================  MPS  ================================== */


static int cmpCodeSum(const void* e1, const void* e2, const void* info)
{
  return( ((int*)info)[*((int*)e1)] < ((int*)info)[*((int*)e2)] ? 1 : 0 );
}


/* ----------------------------------------------------------------- */


static int VectorSum(int* vector, int dim)
{
  int sum = 0;
  int k;

  for( k = 0; k < dim; k++ )
    {
    sum += vector[k];
    }
  return( sum );
}


/* ----------------------------------------------------------------- */


static void CalculateCodevectorSums(CODEBOOK* CB, DATA* D)
{
  int i;

  for( i = 0; i < BookSize(CB); i++ )
    {
    D->codesum[i] = VectorSum(Vector(CB, i), VectorSize(CB));
    }
}


/* ----------------------------------------------------------------- */


static llong SMD(int sum1, int sum2)
{
  return( (llong)(sum1 - sum2) * (llong)(sum1 - sum2) );
}


/* ----------------------------------------------------------------- */


static int NearestBySMD(int vector[],
                        int order[],
                        int start,
                        int end,
                        int value)
{
  int half;

  while( start < end - 1 )
    {
    half = (start + end) / 2;
    if( vector[order[half]] < value )
      {
      start = half;
      }
    else /* value <= vector[order[half]] */
      {
      end = half;
      }
    }

  /* This works also when 'value' is smaller or greater than
     any value in 'vector'. */
  if( value - vector[order[start]] < vector[order[end]] - value )
    {
    return( start );
    }
  else
    {
    return( end );
    }
}

/* ----------------------------------------------------------------- */


static int MPSsearch(BOOKNODE*    v,
                     CODEBOOK*    CB,
                     int*         minindex,
                     llong*       minerror,
                     DISTANCETYPE disttype,
                     int          tsum,
                     int          SMDminindex,
                     int          start,
                     int          end,
                     int          sum[],
                     int          order[])
{
  int   down = SMDminindex - 1;
  int   up   = SMDminindex + 1;

  while( start <= down || up <= end )
    {
    if( start <= down )
      {
      if( SMD(tsum, sum[order[down]]) <= VectorSize(CB) * (*minerror) )
        {
        CheckNearest(v, CB, order[down], minindex, minerror, disttype);
        down--;
        }
      else
        {
        down = start - 1;
        }
      }
    if( up <= end )
      {
      if( SMD(tsum, sum[order[up]]) <= VectorSize(CB) * (*minerror) )
        {
        CheckNearest(v, CB, order[up], minindex, minerror, disttype);
        up++;
        }
      else
        {
        up = end + 1;
        }
      }
    }

  return( *minindex );
}


/* ----------------------------------------------------------------- */


static int FindNearestVectorGLAMPS(BOOKNODE*    v,
                                   int          tvector,
                                   CODEBOOK*    CB,
                                   int          current,
                                   DATA*        D,
                                   DISTANCETYPE disttype)
{
  int    minindex, SMDminindex;
  llong* minerror = &D->nearestd[tvector];
  int    tsum;

  tsum = VectorSum(v->vector, VectorSize(CB));

  SMDminindex = NearestBySMD(D->codesum, D->codeorder, 0, BookSize(CB)-1, tsum);
  minindex = D->codeorder[SMDminindex];
  *minerror = distancef(v->vector, Vector(CB, minindex),
                        VectorSize(CB), MAXLLONG, disttype);

  MPSsearch(v, CB, &minindex, minerror, disttype, tsum,
            SMDminindex, 0, BookSize(CB)-1, D->codesum, D->codeorder);

  return( minindex );
}


/* ----------------------------------------------------------------- */


static void FindNearestVectorsGLAMPS(TRAININGSET*  TS,
                                     CODEBOOK*     CB,
                                     SASchedule*   SAS, /* NULL allowed */
                                     PARTITIONING* P,
                                     DATA*         D,
                                     int*          changedmappings)
{
  int       j;
  int       nearest;
  BOOKNODE  noisevector;
  BOOKNODE* tsvector;

  noisevector = CreateEmptyNode(VectorSize(CB));

  *changedmappings = 0;
  /* Find mapping from training vector to code vector */
  for( j = 0; j < BookSize(TS); j++ )
    {
    SAStoTS(SAS, &Node(TS, j), &noisevector, &tsvector,
            VectorSize(TS), TS->MaxValue);

    nearest = FindNearestVectorGLAMPS(tsvector, j, CB, Map(P, j), D,
                                      distancefunction);
    VerifyNearest(TS, CB, P, j, nearest);

    if( nearest != Map(P, j) )
      {
      ChangePartition(TS, P, nearest, j);
      (*changedmappings)++;
      }
    }

  FreeNode(noisevector);
}


/* ----------------------------------------------------------------- */


static void RunGLAMPS(TRAININGSET*  TS,
                      CODEBOOK*     CB,
                      SASchedule*   SAS,
                      PARTITIONING* P,
                      DATA*         D,
                      int           IterLimit,
                      double*       InitError,
                      double*       FinalError,
                      int*          Iterations)
{
  llong d0, d1 = MAXLLONG, d2;
  int   changedmappings;

  *Iterations = 0;

  CalculateCodevectorSums(CB, D);
  QuickSort(D->codeorder, BookSize(CB), sizeof(int), D->codesum, cmpCodeSum);
  FindNearestVectorsGLAMPS(TS, CB, SAS, P, D, &changedmappings);
  d0 = d2 = CurrentDistortion(TS, D);

  if( Value(ShowProgress) >= 2 )
    {
    PrintIterationInfo(TS, CB, SAS, P, *Iterations, d1, d2, changedmappings,
                       errorfunction, D);
    }

  while( ( d2 < d1 && *Iterations < IterLimit )
         ||
         ( SASInUse(SAS) && SASEffective(SAS, VectorSize(CB), CB->MaxValue) ) )
    {
    ComputeNewBook(TS, CB, SAS, P, D);
    CalculateCodevectorSums(CB, D);
    QuickSort(D->codeorder, BookSize(CB), sizeof(int), D->codesum, cmpCodeSum);

    FindNearestVectorsGLAMPS(TS, CB, SAS, P, D, &changedmappings);
    d1 = d2;
    d2 = CurrentDistortion(TS, D);

    DecreaseTemperature(SAS);
    (*Iterations)++;

    if( Value(ShowProgress) >= 2 )
      {
      PrintIterationInfo(TS, CB, SAS, P, *Iterations, d1, d2, changedmappings,
                         errorfunction, D);
      }
    }

  *InitError  = Total2AverageDistortion(TS, d0);
  *FinalError = Total2AverageDistortion(TS, d2);
}


/* =====================  Local MPS  =============================== */



static void CodevectorSumsPickChanged(CODEBOOK*  CB,
                                      CHANGESET* CS,
                                      int        codeorder[],
                                      int        changedcodeorder[])
{
  int i, l;

  for( i = 0, l = 0; l < CS_size(CS); i++ )
    {
    if( CS_ischanged(CS, codeorder[i]) )
      {
      changedcodeorder[l++] = codeorder[i];
      }
    }
}


/* ----------------------------------------------------------------- */


static int FindNearestVectorGLALocalMPS(BOOKNODE*    v,
                                        int          tvector,
                                        CODEBOOK*    CB,
                                        int          current,
                                        DATA*        D,
                                        DISTANCETYPE disttype)
{
  int    minindex = current;
  int    SMDminindex;
  llong* minerror = &D->nearestd[tvector];
  int    tsum;

  if( CS_size(D->codevector) > 0 )
    {
    tsum = VectorSum(v->vector, VectorSize(CB));

    SMDminindex = NearestBySMD(D->codesum, D->changedcodeorder,
                               0, CS_size(D->codevector)-1, tsum);
    if( D->changedcodeorder[SMDminindex] != current )
      {
      CheckNearest(v, CB, D->changedcodeorder[SMDminindex],
                   &minindex, minerror, disttype);
      }

    MPSsearch(v, CB, &minindex, minerror, disttype, tsum,
              SMDminindex, 0, CS_size(D->codevector)-1,
              D->codesum,  D->changedcodeorder);
    }

  return( minindex );
}


/* ----------------------------------------------------------------- */


static void FindNearestVectorsGLALocalMPS(TRAININGSET*  TS,
                                          CODEBOOK*     CB,
                                          SASchedule*   SAS, /* NULL allowed */
                                          PARTITIONING* P,
                                          DATA*         D,
                                          int*          changedmappings)
{
  int       j;
  int       nearest;
  BOOKNODE  noisevector;
  BOOKNODE* tsvector;

  noisevector = CreateEmptyNode(VectorSize(CB));
  CalculateCurrentDistances(TS, CB, P, D, distancefunction);

  *changedmappings = 0;
  /* Find mapping from training vector to code vector */
  for( j = 0; j < BookSize(TS); j++ )
    {
    SAStoTS(SAS, &Node(TS, j), &noisevector, &tsvector,
            VectorSize(TS), TS->MaxValue);

    if( CS_ischanged(D->codevector, Map(P, j)) && ! D->nowcloser[j] )
      {
      nearest = FindNearestVectorGLAMPS(tsvector, j, CB, Map(P, j), D,
                                        distancefunction);
      }
    else
      {
      nearest = FindNearestVectorGLALocalMPS(tsvector, j, CB, Map(P, j), D,
                                             distancefunction);
      }
    VerifyNearest(TS, CB, P, j, nearest);

    if( nearest != Map(P, j) )
      {
      ChangePartition(TS, P, nearest, j);
      (*changedmappings)++;
      }
    }

  FreeNode(noisevector);
}


/* ----------------------------------------------------------------- */


static void RunGLALocalMPS(TRAININGSET*  TS,
                           CODEBOOK*     CB,
                           SASchedule*   SAS,
                           PARTITIONING* P,
                           DATA*         D,
                           int           IterLimit,
                           double*       InitError,
                           double*       FinalError,
                           int*          Iterations)
{
  llong d0, d1 = MAXLLONG, d2;
  int   changedmappings;

  *Iterations = 0;

  CalculateCodevectorSums(CB, D);
  QuickSort(D->codeorder, BookSize(CB), sizeof(int), D->codesum, cmpCodeSum);

  FindNearestVectorsGLAMPS(TS, CB, SAS, P, D, &changedmappings);
  d0 = d2 = CurrentDistortion(TS, D);

  if( Value(ShowProgress) >= 2 )
    {
    PrintIterationInfo(TS, CB, SAS, P, *Iterations, d1, d2, changedmappings,
                       errorfunction, D);
    }

  while( ( d2 < d1 && *Iterations < IterLimit )
         ||
         ( SASInUse(SAS) && SASEffective(SAS, VectorSize(CB), CB->MaxValue) ) )
    {
    ComputeNewBook(TS, CB, SAS, P, D);
    CalculateCodevectorSums(CB, D);
    QuickSort(D->codeorder, BookSize(CB), sizeof(int), D->codesum, cmpCodeSum);
    CodevectorSumsPickChanged(CB, D->codevector, D->codeorder, D->changedcodeorder);

    FindNearestVectorsGLALocalMPS(TS, CB, SAS, P, D, &changedmappings);
    d1 = d2;
    d2 = CurrentDistortion(TS, D);

    DecreaseTemperature(SAS);
    (*Iterations)++;

    if( Value(ShowProgress) >= 2 )
      {
      PrintIterationInfo(TS, CB, SAS, P, *Iterations, d1, d2, changedmappings,
                         errorfunction, D);
      }
    }

  *InitError  = Total2AverageDistortion(TS, d0);
  *FinalError = Total2AverageDistortion(TS, d2);
}


/* ================================================================= */


static void InitializeGeneral(TRAININGSET*   TS,
                              CODEBOOK*      CB,
                              PARTITIONING*  OrigP,
                              PARTITIONING** P,
                              DATA*          D)
{
  int i;

  if( OrigP != NULL )
    {
    /* Use already initialized partitioning. */
    *P = OrigP;
    }
  else
    {
    /* All training vectors are mapped to first partition (0). */
    *P = allocate(sizeof(PARTITIONING));
    CreateNewPartitioning(*P, TS, BookSize(CB));
    }

  if( errorfunction == MSE )
    {
    switch( Value(UsePDS) )
      {
      case PDSoff: distancef = VectorDistanceTrivial; break;
      case PDSon:  distancef = VectorDistancePDS;     break;
      default:     ErrorMessage("Unknown UsePDS=%i\n", Value(UsePDS));
                   ExitProcessing( -1 );
      }
    }
  else
    {
    distancef = VectorDistance; /* from cb.c */
    }

  if( MethodTracksChangedCodeVectors() ||
      Value(ShowProgress) >= 3 )
    {
    CreateNewCodebook(&D->prevCB, BookSize(CB), TS);
    CopyCodebook(CB, &D->prevCB);

    D->codevector = CS_make(BookSize(CB));
    D->changedclustersize = allocate(BookSize(CB) * sizeof(int));
    for( i = 0; i < BookSize(CB); i++ )
      {
      CS_insert(D->codevector, i);
      D->changedclustersize[i] = CCFreq(*P, i);
      }
    D->nowcloser = allocate(BookSize(TS) * sizeof(YESNO));
    }
  else
    {
    D->codevector         = NULL;
    D->changedclustersize = NULL;
    D->nowcloser          = NULL;
    }

  D->nearestd = allocate(BookSize(TS) * sizeof(llong));
  for( i = 0; i < BookSize(TS); i++ )
    {
    D->nearestd[i] = 0LL;
    }

  if( Value(GLAMethod) == GLAWu ||
      Value(GLAMethod) == GLATIE ||
      Value(GLAMethod) == GLALocalTIE )
    {
    DistanceMatrixMake(CB, &D->DM);
    OrderMatrixMake(CB, &D->sDM);
    OrderMatrixMake(CB, &D->scDM);
    }

  if( Value(GLAMethod) == GLAMPS ||
      Value(GLAMethod) == GLALocalMPS )
    {
    D->codesum   = allocate(BookSize(CB) * sizeof(int));
    D->codeorder = allocate(BookSize(CB) * sizeof(int));
    for( i = 0; i < BookSize(CB); i++ )
      {
      D->codesum[i] = 0;
      D->codeorder[i] = i;
      }
    D->changedcodeorder = allocate(BookSize(CB) * sizeof(int));
    }
  else
    {
    D->codesum   = NULL;
    D->codeorder = NULL;
    D->changedcodeorder = NULL;
    }
}


/* ----------------------------------------------------------------- */


static void ShutDownGeneral(CODEBOOK*     CB,
                            PARTITIONING* OrigP,
                            PARTITIONING* P,
                            DATA*         D)
{
  if( OrigP == NULL )
    {
    FreePartitioning(P);
    deallocate(P);
    }

  if( MethodTracksChangedCodeVectors() ||
      Value(ShowProgress) >= 3 )
    {
    FreeCodebook(&D->prevCB);
    CS_free(D->codevector);
    deallocate(D->changedclustersize);
    deallocate(D->nowcloser);
    }

  deallocate(D->nearestd);

  if( Value(GLAMethod) == GLAWu ||
      Value(GLAMethod) == GLATIE ||
      Value(GLAMethod) == GLALocalTIE )
    {
    DistanceMatrixFree(CB, &D->DM);
    OrderMatrixFree(CB, &D->sDM);
    OrderMatrixFree(CB, &D->scDM);
    }

  if( Value(GLAMethod) == GLAMPS ||
      Value(GLAMethod) == GLALocalMPS )
    {
    deallocate(D->codesum);
    deallocate(D->codeorder);
    deallocate(D->changedcodeorder);
    }
}


/* ----------------------------------------------------------------- */


void GLA(TRAININGSET*  TS,
         CODEBOOK*     CB,
         SASchedule*   SAS,
         PARTITIONING* OrigP,
         int           IterLimit,
         double*       InitError,
         double*       FinalError,
         int*          Iterations)
{
  PARTITIONING* P;
  DATA          D;

  if( IterLimit > 0 )
    {
    InitializeGeneral(TS, CB, OrigP, &P, &D);
    switch( Value(GLAMethod) )
      {
      case GLAFullSearch:
        {
        RunGLAFullSearch(TS, CB, SAS, P, &D, IterLimit,
                         InitError, FinalError, Iterations);
        break;
        }
      case GLALocal:
        {
        RunGLALocal(TS, CB, SAS, P, &D, IterLimit,
                    InitError, FinalError, Iterations);
        break;
        }
      case GLAWu:
        {
        RunGLAWu(TS, CB, SAS, P, &D, IterLimit,
                 InitError, FinalError, Iterations);
        break;
        }
      case GLATIE:
        {
        RunGLATIE(TS, CB, SAS, P, &D, IterLimit,
                  InitError, FinalError, Iterations);
        break;
        }
      case GLALocalTIE:
        {
        RunGLALocalTIE(TS, CB, SAS, P, &D, IterLimit,
                       InitError, FinalError, Iterations);
        break;
        }
      case GLAMPS:
        {
        RunGLAMPS(TS, CB, SAS, P, &D, IterLimit,
                    InitError, FinalError, Iterations);
        break;
        }
      case GLALocalMPS:
        {
        RunGLALocalMPS(TS, CB, SAS, P, &D, IterLimit,
                       InitError, FinalError, Iterations);
        break;
        }
      default:
        {
        ErrorMessage("Unknown GLAMethod=%i\n", Value(GLAMethod));
        ExitProcessing( -1 );
        }
      }
    ShutDownGeneral(CB, OrigP, P, &D);
    }

  if( Value(ShowProgress) >= 3 )
    {
    /* Statistical information */
    printf("Total #d= %7.0f #sq= %7.0f avg= %5.2f\n",
           (double)ndistcalcstotal,
           (double)ksumtotal,
           ndistcalcstotal != 0LL ?
             (double)ksumtotal/(double)ndistcalcstotal :
             0.0);
    }
}

