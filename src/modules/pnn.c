/*-------------------------------------------------------------------*/
/* PNN.C          Timo Kaukoranta                                    */
/*                                                                   */
/*                                                                   */
/* - Module of PNN method                                            */
/*                                                                   */
/*-------------------------------------------------------------------*/

#define ProgName       "PNN"
#define VersionNumber  "V.0.20a"
#define LastUpdated    "28.7.99"  /* P.F. */

/* ----------------------------------------------------------------- */

#include <assert.h>
#include <stdio.h>
#include <values.h>

#include "cb.h"
#include "heap.h"
#include "interfc.h"
#include "memctrl.h"
#include "owntypes.h"
#include "pnn.h"
#include "sa.h"
#include "stack.h"

#define MERGEPERCENT 50

/* -------------------- Parameters for the PNN module -------------- */

static struct { int ShowProgress;
                int MergeErrorType;
                int CentroidCalculation;
                int MergeMethod;
                int PartitionRemapping;
                int GLAIterations;
                int MergePercent;
                int MaxBucket;
                int MinBucket;
              } ModuleParameters;

#define Value(x)      (ModuleParameters.x)

/* These enum definitions are from the *.fac file. */
enum { Equitz, TrainingSetAmount, TrainingSetIncrease };
enum { CodevectorBased, PartitionBased };
enum { Kissing, TrivialON3, Kurita, LazyPNN, FastPNN };

/*-----------------------------  M i s c  ----------------------------*/


typedef struct { llong       VectorDist;
                 llong       LocalDist;
                 int         Vector;
                 int         Recalculate;
                 int         WhoAmI;
                 int         HeapIndex;
               } NEARESTTYPE;


struct KDNODESTRUCT { int                    dim;
                      VECTORELEMENT          threshold;
                      struct KDNODESTRUCT*   left;
                      struct KDNODESTRUCT*   right;
                      struct KDNODESTRUCT*   father;
                      int                    bucket;
                    };

struct KDBUCKETSTRUCT { int   size;
                        int   first;
                        int   cand1;
                        int   cand2;
                        llong dist;
                        struct KDNODESTRUCT* node;
                      };

typedef struct KDNODESTRUCT   KDNODE;
typedef struct KDBUCKETSTRUCT KDBUCKET;

typedef struct { KDNODE*   root;
                 KDBUCKET* bucket;
                 int*      map;
                 int*      next;
                 int*      bucketorder;
                 int       inorder;
                 YESNO*    toberemoved;
                 int       nfreebucket;
                 int*      freebucket;
               } KDTREE;


typedef struct { llong dist;
                 int   from, to;
                 int   heapindex;
               } DISTELEM;

typedef struct { DISTELEM** DM;
                 HEAP       H;
                 int*       matrix2cb;
                 int*       cb2matrix;
               } KURITADATA;

#define cb2m(KD,i)    ((KD)->cb2matrix[i])
#define m2cb(KD,i)    ((KD)->matrix2cb[i])


#define round(a)      ((a) + 0.5)
#define fpositive(a)  ((a) > FLT_EPSILON)
#define swapint(a,b)  { int tmp = (a); (a) = (b); (b) = tmp; }

int   origCBsize = 0;
YESNO MaintainPartitions = NO;


/* =================== Parameters for the PNN module =============== */


void SetPNNParameters(int ShowProgress,
                      int MergeErrorType,
                      int CentroidCalculation,
                      int MergeMethod,
                      int PartitionRemapping,
                      int GLAIterations,
                      int MergePercent,
                      int MaxBucket,
                      int MinBucket)
{
  Value(ShowProgress)        = ShowProgress;
  Value(MergeErrorType)      = MergeErrorType;
  Value(CentroidCalculation) = CentroidCalculation;
  Value(MergeMethod)         = MergeMethod;
  Value(PartitionRemapping)  = PartitionRemapping;
  Value(GLAIterations)       = GLAIterations;
  Value(MergePercent)        = MergePercent;
  Value(MaxBucket)           = MaxBucket;
  Value(MinBucket)           = MinBucket;
}


/* ================================================================= */

#define Scale(KD, from, to) \
  int i = cb2m(KD, from), j =  cb2m(KD, to), tmp; \
  if( i < j ) { tmp = i; i = j;  j = tmp; }


/* ----------------------------------------------------------------- */

static DISTELEM* DM_elem(KURITADATA* KD, int from, int to)
{
  Scale(KD, from, to);
  return( &(KD->DM[i][j]) );
}

/* ----------------------------------------------------------------- */

static int DM_from(KURITADATA* KD, int from, int to)
{
  Scale(KD, from, to);
  return( m2cb(KD, KD->DM[i][j].from) );
}

/* ----------------------------------------------------------------- */

static int DM_to(KURITADATA* KD, int from, int to)
{
  Scale(KD, from, to);
  return( m2cb(KD, KD->DM[i][j].to) );
}

/* ----------------------------------------------------------------- */

static int DM_heapindex(KURITADATA* KD, int from, int to)
{
  Scale(KD, from, to);
  return( KD->DM[i][j].heapindex );
}

/* ----------------------------------------------------------------- */

static int* DM_heapindexaddr(KURITADATA* KD, int from, int to)
{
  Scale(KD, from, to);
  return( &(KD->DM[i][j].heapindex) );
}

/* ----------------------------------------------------------------- */

static int DM_dist(KURITADATA* KD, int from, int to)
{
  Scale(KD, from, to);
  return( KD->DM[i][j].dist );
}

/* ----------------------------------------------------------------- */

static void DM_distset(KURITADATA* KD, int from, int to, llong dist)
{
  Scale(KD, from, to);
  KD->DM[i][j].dist = dist;
/*   PrintMessage("Dset f,t=%i,%i -> i,j=%i,%i\n", from, to, i, j); */
}


/* ================================================================= */


int CmpDistKurita(void* a, void* b, void* info)
{
  if( ((DISTELEM*)a)->dist > ((DISTELEM*)b)->dist )
    {
    return( 1 );
    }
  else
    {
    if( ((DISTELEM*)a)->dist == ((DISTELEM*)b)->dist )
      {
      return( 0 );
      }
    else
      {
      return( -1 );
      }
    }
}


/* ----------------------------------------------------------------- */


int CmpDistLazy(void* a, void* b, void* info)
{
  llong diff = ((NEARESTTYPE*)a)->VectorDist - ((NEARESTTYPE*)b)->VectorDist;

  return( diff > 0 ? 1 : diff == 0 ? 0 : -1 );
}


/* ================================================================= */


static int GetBucketIndex(KDTREE* tree)
{
  assert( 0 < tree->nfreebucket );

  return( tree->freebucket[--(tree->nfreebucket)] );
}


/* ----------------------------------------------------------------- */


static void ReleaseBucketIndex(KDTREE* tree, int index)
{
  tree->freebucket[tree->nfreebucket] = index;
  tree->nfreebucket++;
}


/* ----------------------------------------------------------------- */


static KDNODE* NewTreeNode(void)
{
  KDNODE* node = allocate(sizeof(KDNODE));

  node->dim       = -1;
  node->threshold = -1;
  node->left      = NULL;
  node->right     = NULL;
  node->father    = NULL;
  node->bucket    = -1;

  return( node );
}


/* ----------------------------------------------------------------- */


static void FreeTreeNode(KDNODE* node)
{
  deallocate(node);
}


/* ----------------------------------------------------------------- */


static void AddVectorToBucket(KDTREE* tree, int vindex, int bindex)
{
  KDBUCKET* buc = &tree->bucket[bindex];

  tree->map[vindex] = bindex;
  tree->next[vindex] = buc->first;
  buc->first = vindex;
  buc->size++;
}

/* ----------------------------------------------------------------- */


static void RemoveVectorFromBucket(KDTREE* tree, int vindex)
{
  KDBUCKET* buc = &tree->bucket[tree->map[vindex]];
  int*      prev = &buc->first;
  int       i = buc->first;

  while( i != -1 && i != vindex )
    {
    prev = &tree->next[i];
    i = tree->next[i];
    }
  if( i == vindex )
    {
    *prev = tree->next[i];
    tree->next[i] = -1;
    tree->map[i] = -1;
    buc->size--;
    }
  else
    {
    ErrorMessage("RemoveVectorFromBucket: '%i' is not in the bucket '%i'\n",
             vindex, tree->map[vindex]);
    ExitProcessing( -1 );
    }
}


/* ----------------------------------------------------------------- */


static void JoinBuckets(KDTREE* tree, int to, int from)
{
  KDBUCKET* bucto = &tree->bucket[to];
  KDBUCKET* bucfrom = &tree->bucket[from];
  int       i = bucfrom->first;

  if( i != -1 )
    {
    tree->map[i] = to;
    while( tree->next[i] != -1 )
      {
      i = tree->next[i];
      tree->map[i] = to;
      }
    tree->next[i] = bucto->first;
    bucto->first = bucfrom->first;
    bucfrom->first = -1;
    }
  bucto->size += bucfrom->size;
  bucfrom->size = 0;
  ReleaseBucketIndex(tree, from);
}


/* ----------------------------------------------------------------- */


static void MergeSubtree(KDTREE* tree, KDNODE* subroot)
{
  STACK*    S = S_make();
  KDNODE*   node = subroot;
  int       tobucket = -1;

  S_push(S, subroot);
  while( ! S_empty(S) )
    {
    node = S_pop(S);
    if( node->bucket == -1 )
      {
      assert( node->right != NULL );
      assert( node->left != NULL );
      S_push(S, node->right);
      S_push(S, node->left);
      }
    else /* We have a bucket. */
      {
      if( tobucket == -1 )
        {
        tobucket = node->bucket;
        }
      else
        {
        JoinBuckets(tree, tobucket, node->bucket);
        }
      }
    if( node != subroot )
      {
      FreeTreeNode(node);
      }
    }
  subroot->left = NULL;
  subroot->right = NULL;
  subroot->bucket = tobucket;

  S_free(S);
}


/* ----------------------------------------------------------------- */


static void CalculateVariance(CODEBOOK* CB,
                              KDTREE*   tree,
                              KDNODE*   node,
                              int       dim,
                              double*   var)
{
  int    i;
  llong  sum = 0LL;
  llong  sqsum = 0LL;
  int    freq = 0;
  llong  elem;
  double avg;

  for( i = tree->bucket[node->bucket].first; i != -1; i = tree->next[i] )
    {
    elem = VectorScalar(CB, i, dim);
    sum += elem * VectorFreq(CB, i);
    sqsum += elem * elem * VectorFreq(CB, i);
    freq += VectorFreq(CB, i);
    }

  avg = (double)sum / (double)freq;
  *var = (double)sqsum / (double)freq - avg * avg;
}


/* ----------------------------------------------------------------- */

static void CalculateMedian(CODEBOOK* CB,
                            KDTREE*   tree,
                            KDNODE*   node,
                            int       dim,
                            int*      median)
{
  KDBUCKET* buc = &tree->bucket[node->bucket];
  int*   value = allocate(buc->size * sizeof(int));
  int    i, j;
  int    start = 0;
  int    end = buc->size-1;
  int    medianindex = (start + end) / 2;
  int    left, right;
  int    pivot, tmp;

  for( i = buc->first, j = 0; i != -1; i = tree->next[i], j++ )
    {
    value[j] = VectorScalar(CB, i, dim);
    }

  left = start;
  right = end;
  while( left <= medianindex || medianindex <= right )
    {
    left = start;
    right = end;
    pivot = value[(start + end) / 2];
    while( left <= right )
      {
      while( value[left] < pivot )
        {
        left++;
        }
      while( pivot < value[right] )
        {
        right--;
        }
      if( left <= right )
        {
        tmp = value[left];
        value[left] = value[right];
        value[right] = tmp;
        left++;
        right--;
        }
      }
    if( left <= medianindex )
      {
      start = left;
      }
    if( medianindex <= right )
      {
      end = right;
      }
    }

  *median = value[medianindex];
  deallocate(value);
}


/* ----------------------------------------------------------------- */


static void DivideBucket(CODEBOOK* CB,
                         KDTREE*   tree,
                         KDNODE*   node,
                         KDNODE**  left,
                         KDNODE**  right)
{
  int           i, k;
  int           maxdim = 0;
  double        maxvar, var;
  VECTORELEMENT threshold;
  int*          prev;

  CalculateVariance(CB, tree, node, 0, &maxvar);
  for( k = 1; k < VectorSize(CB); k++ )
    {
    CalculateVariance(CB, tree, node, k, &var);
    if( maxvar < var )
      {
      maxvar = var;
      maxdim = k;
      }
    }
  CalculateMedian(CB, tree, node, maxdim, &threshold);

  node->dim = maxdim;
  node->threshold = threshold;
  *left  = node->left  = NewTreeNode();
  *right = node->right = NewTreeNode();
  (*left)->father  = node;
  (*right)->father = node;
  (*left)->bucket  = node->bucket;
  (*right)->bucket = GetBucketIndex(tree);
  node->bucket = -1;
  tree->bucket[(*right)->bucket].size = 0;
  tree->bucket[(*right)->bucket].first = -1;

  /* Split code vectos into two buckets. */
  prev = &tree->bucket[(*left)->bucket].first;
  i = tree->bucket[(*left)->bucket].first;

  while( i != -1 &&
         tree->bucket[(*right)->bucket].size <
         tree->bucket[(*left)->bucket].size)
    {
    if( VectorScalar(CB, i, maxdim) >= threshold  )
      {
      tree->map[i] = (*right)->bucket;
      *prev = tree->next[i];
      tree->next[i] = tree->bucket[(*right)->bucket].first;
      tree->bucket[(*right)->bucket].first = i;
      tree->bucket[(*left)->bucket].size--;
      tree->bucket[(*right)->bucket].size++;
      i = *prev;
      }
    else
      {
      prev = &tree->next[i];
      i = tree->next[i];
      }
    }
}


/* ----------------------------------------------------------------- */


static void DivideSubtree(CODEBOOK* CB, KDTREE* tree, KDNODE* subtree)
{
  STACK*  S = S_make();
  KDNODE* node;
  KDNODE* left;
  KDNODE* right;

  S_push(S, subtree);
  while( ! S_empty(S) )
    {
    node = S_pop(S);
    if( tree->bucket[node->bucket].size > Value(MaxBucket) )
      {
      DivideBucket(CB, tree, node, &left, &right);
      S_push(S, right);
      S_push(S, left);
      }
    }

  S_free(S);
}


/* ----------------------------------------------------------------- */


static void MakeKDTree(CODEBOOK* CB, KDTREE* tree)
{
  int i;
  int bindex;

  for( i = 0; i < BookSize(CB); i++ )
    {
    tree->freebucket[i] = i;
    }
  tree->nfreebucket = BookSize(CB);

  /* Put all code vectors into the same bucket. */
  bindex = GetBucketIndex(tree);
  tree->root->bucket         = bindex;
  tree->bucket[bindex].size  = BookSize(CB);
  tree->bucket[bindex].first = 0;
  tree->bucket[bindex].cand1 = -1;
  tree->bucket[bindex].cand2 = -1;
  tree->bucket[bindex].node  = NULL;

  for( i = 0; i < BookSize(CB); i++ )
    {
    tree->map[i] = bindex;
    tree->next[i] = i+1;
    tree->toberemoved[i] = NO;
    }
  tree->next[BookSize(CB)-1] = -1;

  DivideSubtree(CB, tree, tree->root);
}


/* ----------------------------------------------------------------- */

static void FreeKDTree(KDTREE* tree)
{
  STACK*  S = S_make();
  KDNODE* node;

  S_push(S, tree->root);
  while( ! S_empty(S) )
    {
    node = S_pop(S);
    if( node->bucket == -1 )
      {
      S_push(S, node->right);
      S_push(S, node->left);
      }
    FreeTreeNode(node);
    }
  S_free(S);

  deallocate(tree->bucket);
  deallocate(tree->map);
  deallocate(tree->next);
  deallocate(tree->bucketorder);
  deallocate(tree->toberemoved);
  deallocate(tree->freebucket);
}


/* ----------------------------------------------------------------- */


static void BalanceKDTree(CODEBOOK* CB, KDTREE* tree)
{
  STACK*  S = S_make();
  KDNODE* node;

  S_push(S, tree->root);
  while( ! S_empty(S) )
    {
    node = S_pop(S);
    if( node->bucket == -1 )
      {
      assert( node->right != NULL );
      assert( node->left != NULL );
      S_push(S, node->right);
      S_push(S, node->left);
      }
    else /* We have a bucket. */
      {
      /* Is the bucket too small and it is not a root? */
      if( tree->bucket[node->bucket].size < Value(MinBucket) &&
          node->father != NULL )
        {
        if( node->father->left == node )
          {
          S_pop(S); /* Remove the right sibling from the stack. */
          }
        node = node->father;
        MergeSubtree(tree, node);
        DivideSubtree(CB, tree, node);
        }
      }
    }

  S_free(S);
}


/* ================================================================= */


static void InitializePNN(TRAININGSET*   TS,
                          CODEBOOK*      CB,
                          PARTITIONING*  OrigP,
                          PARTITIONING** P,
                          NEARESTTYPE**  NN,
                          KURITADATA*    KD,
                          HEAP*          H,
                          KDTREE*        tree)
{
  int i, j;

  origCBsize = BookSize(CB);

  if( OrigP != NULL )
    {
    *P = OrigP;
    }
  else
    {
    *P = allocate(sizeof(PARTITIONING));
    CreateNewPartitioning(*P, TS, BookSize(CB));
    }

  *NN = (NEARESTTYPE*)allocate(BookSize(CB) * sizeof(NEARESTTYPE));
  for( i = 0; i < BookSize(CB); i++ )
    {
    (*NN)[i].Vector      = 0;
    (*NN)[i].VectorDist  = 0LL;
    (*NN)[i].LocalDist   = 0LL;
    (*NN)[i].Recalculate = 0;
    (*NN)[i].WhoAmI      = i;
    (*NN)[i].HeapIndex   = 0;
    }

  switch( Value(MergeMethod) )
    {
    case Kissing:
    case TrivialON3:
      {
      break;
      }
    case Kurita:
      {
      Heap_init(&(KD->H), (BookSize(CB)*(BookSize(CB)+1))/2, CmpDistKurita);

      KD->DM = (DISTELEM**)allocate(BookSize(CB) * sizeof(DISTELEM*));
      KD->DM[0] = NULL;
      for( i = 1; i < BookSize(CB); i++ )
        {
        KD->DM[i] = (DISTELEM*)allocate(i * sizeof(DISTELEM));
        for( j = 0; j < i; j++ )
          {
          KD->DM[i][j].from = i;
          KD->DM[i][j].to = j;
          KD->DM[i][j].heapindex = 0;
          KD->DM[i][j].dist = 0LL;
          }
        }

      KD->matrix2cb = (int*)allocate(BookSize(CB) * sizeof(int));
      KD->cb2matrix = (int*)allocate(BookSize(CB) * sizeof(int));
      for( i = 0; i < BookSize(CB); i++ )
        {
        KD->matrix2cb[i] = i;
        KD->cb2matrix[i] = i;
        }
      break;
      }
    case LazyPNN:
      {
      Heap_init(H, BookSize(CB), CmpDistLazy);
      break;
      }
    case FastPNN:
      {
      tree->root   = NewTreeNode();
      tree->bucket = allocate(BookSize(CB) * sizeof(KDBUCKET));
      tree->map    = allocate(BookSize(CB) * sizeof(int));
      tree->next   = allocate(BookSize(CB) * sizeof(int));
      tree->bucketorder = allocate(BookSize(CB) * sizeof(int));
      tree->toberemoved = allocate(BookSize(CB) * sizeof(YESNO));
      tree->freebucket  = allocate(BookSize(CB) * sizeof(int));
      break;
      }
    default:
      break;
    }

  /* Should we maintain the partitions? */
  if( Value(MergeErrorType) == TrainingSetAmount ||
      Value(MergeErrorType) == TrainingSetIncrease ||
      Value(CentroidCalculation) == PartitionBased ||
      Value(PartitionRemapping) == YES ||
      Value(GLAIterations) == YES )
    {
    MaintainPartitions = YES;
    }
  else
    {
    MaintainPartitions = NO;
    }

}


/* ----------------------------------------------------------------- */


static void ShutDownPNN(TRAININGSET*  TS,
                        CODEBOOK*     CB,
                        PARTITIONING* OrigP,
                        PARTITIONING* P,
                        NEARESTTYPE   NN[],
                        KURITADATA*   KD,
                        HEAP*         H,
                        KDTREE*       tree)
{
  int i;

  if( OrigP == NULL )
    {
    FreePartitioning(P);
    deallocate(P);
    }
  else
    {
    /* We have original partitions but we have not maintained them. */
    if( MaintainPartitions == NO )
      {
      FreePartitioning(OrigP);
      CreateNewPartitioning(OrigP, TS, BookSize(CB));
      GenerateOptimalPartitioning(TS, CB, P);
      }
    }

  deallocate(NN);
  switch( Value(MergeMethod) )
    {
    case Kissing:
    case TrivialON3:
      {
      break;
      }
    case Kurita:
      {
      for( i = 1; i < origCBsize; i++ )
        {
        deallocate(KD->DM[i]);
        }
      deallocate(KD->DM);
      deallocate(KD->matrix2cb);
      deallocate(KD->cb2matrix);
      Heap_free(&(KD->H));
      break;
      }
    case LazyPNN:
      {
      Heap_free(H);
      break;
      }
    case FastPNN:
      {
      FreeKDTree(tree);
      break;
      }
    default:
      break;
    }
}


/*=========================  CODEBOOK HANDLING  ============================*/


static void PrintPartition(TRAININGSET* TS, PARTITIONING* P, int Pindex)
{
  int j;

  PrintMessage("PP%5u: Freq=%5i Unique=%5i\n",
               Pindex, CCFreq(P,Pindex), UniqueVectors(P, Pindex));
  for( j = FirstVector(P,Pindex); ! EndOfPartition(j); j = NextVector(P,j) )
    {
    PrintMessage("j=%5i fr=%5i: ", j, VectorFreq(TS, j));
    PrintVector(Vector(TS, j), VectorSize(TS), 1);
    }
  PrintMessage("\n");
}


/*===================================================================*/


static void GenerateInitialPartition(TRAININGSET*  TS,
                                     CODEBOOK*     CB,
                                     PARTITIONING* P,
                                     NEARESTTYPE   NN[])
{
  int   i, j;
  llong error;

  /* Although MaintainPartitions==NO we guarantee correct frequencies
     of the codevectors (at the end of this procedure)
     by following calculations */

  if( BookSize(TS) == BookSize(CB) )
    {
    /* No initial codebook, start from TS. */
    PutAllInOwnPartition(TS, P);
    }
  else
    {
    /* We have an initial codebook. */

    switch( Value(MergeErrorType) )
      {
      case Equitz:
        {
        break;
        }
      case TrainingSetAmount:
      case TrainingSetIncrease:
        {
        GenerateOptimalPartitioning(TS, CB, P);

        for( j = 0; j < BookSize(TS); j++ )
          {
          error = VectorDistance(Vector(TS, j), Vector(CB, Map(P, j)),
                                 VectorSize(TS), MAXLLONG, EUCLIDEANSQ)
                  * VectorFreq(TS, j);
          NN[Map(P, j)].LocalDist += error;
          }
        break;
        }
      }
    }

  /* Calculations of the distortion increase are based on the
     frequencies of the codevectors. Therefore, we have to ensure
     correctness of them. */
  for( i = 0; i < BookSize(CB); i++ )
    {
    VectorFreq(CB, i) = CCFreq(P, i);
    }
}


/* ================================================================= */


static llong PartitionError(TRAININGSET*  TS,
                            PARTITIONING* P,
                            int           Pindex,
                            BOOKNODE*     vec,
                            llong         MinDist)
{
  int   j;
  llong error = 0LL;

  for( j = FirstVector(P,Pindex);
       ! EndOfPartition(j) && error < MinDist;
       j = NextVector(P,j) )
    {
    error += VectorDistance(vec->vector, Vector(TS, j), VectorSize(TS),
                            MinDist, EUCLIDEANSQ)
             * VectorFreq(TS, j);

    if( Value(ShowProgress) >= 4 )
      {
      PrintMessage("PartitionError: "
                   "Pindex=%i j=%i error=%.0f freq=%i MD=%.0f\n",
                   Pindex, j, (double)error,
                   VectorFreq(TS, j), (double)MinDist);
      }
    }

  assert( error >= 0LL );
  
  return( error );
}


/* ----------------------------------------------------------------- */


static void CalculatePartitionErrors(TRAININGSET*  TS,
                                     CODEBOOK*     CB,
                                     PARTITIONING* P,
                                     NEARESTTYPE   NN[])
{
  int i;

  for( i = 0; i < BookSize(CB); i++ )
    {
    NN[i].LocalDist = PartitionError(TS, P, i, &Node(CB, i), MAXLLONG);
    assert( NN[i].LocalDist >= 0LL );
    }
}


/* ----------------------------------------------------------------- */


static llong ErrorOfPartitions(TRAININGSET*  TS,
                               CODEBOOK*     CB,
                               PARTITIONING* P)
{
  int   j;
  llong error = 0LL;

  for( j = 0; j < BookSize(TS); j++ )
    {
    error += VectorDistance(Vector(TS, j), Vector(CB, Map(P, j)),
                            VectorSize(TS), MAXLLONG, EUCLIDEANSQ)
             * VectorFreq(TS, j);
    }

  return( error );
}


/* ----------------------------------------------------------------- */


static double SE2MSE(TRAININGSET* TS, double error)
{
  assert( TotalFreq(TS) > 0 );

  return( error / ((double)TS->TotalFreq * (double)VectorSize(TS)) );
}


/* ----------------------------------------------------------------- */


static void CentroidOf2Vectors(int       p1,
                               int       p2,
                               CODEBOOK* CB,
                               BOOKNODE* cent)
{
  int k;
  int freqsum;

  assert( 0 <= VectorFreq(CB, p1) );
  assert( 0 <= VectorFreq(CB, p2) );

  if( Value(ShowProgress) >= 5 )
    {
    PrintMessage("\nCentroidOf2Vectors:%i & %i size=%i\n",
                 p1, p2, BookSize(CB));
    }

  freqsum = VectorFreq(CB, p1) + VectorFreq(CB, p2);

  /* Calculate mean values for centroid */
  for( k = 0; k < VectorSize(CB); k++ )
    {
    if( freqsum > 0 )
      {
      cent->vector[k] =
        round( (double)(VectorFreq(CB, p1)*VectorScalar(CB, p1, k) +
                        VectorFreq(CB, p2)*VectorScalar(CB, p2, k))
                       / (double)freqsum );
      }
    else
      {
      cent->vector[k] = 0;
      }
    }
  cent->freq = freqsum;
}


/* ----------------------------------------------------------------- */


static void CentroidOf2Partitions(PARTITIONING* P,
                                  int           p1,
                                  int           p2,
                                  int           Vsize,
                                  BOOKNODE*     cent)
{
  int k;
  int freq;

  assert( 0 <= CCFreq(P, p1) );
  assert( 0 <= CCFreq(P, p2) );

  if( Value(ShowProgress) >= 5 )
    {
    PrintMessage("\nCentroidOf2Partitions:%i & %i\n", p1, p2);
    }

  freq = CCFreq(P, p1) + CCFreq(P, p2);

  /* Calculate mean values for centroid */
  for( k = 0; k < Vsize; k++ )
    {
    if( freq > 0 )
      {
      cent->vector[k] =
        round( (double)(CCScalar(P, p1, k) + CCScalar(P, p2, k))
                       / (double)freq );
      }
    else
      {
      cent->vector[k] = 0;
      }
    }
  cent->freq = freq;
}


/* ----------------------------------------------------------------- */


static void CalculateNewCentroid(int           a,
                                 int           b,
                                 CODEBOOK*     CB,
                                 TRAININGSET*  TS,
                                 PARTITIONING* P,
                                 BOOKNODE*     cent)
{
  switch( Value(CentroidCalculation) )
    {
    case CodevectorBased:
      {
      CentroidOf2Vectors(a, b, CB, cent);
      break;
      }
    case PartitionBased:
      {
      CentroidOf2Partitions(P, a, b, VectorSize(CB), cent);
      break;
      }
    default:
      {
      ErrorMessage("Unknown CentroidCalculation=%i", Value(CentroidCalculation));
      ExitProcessing( -1 );
      }
    }
}


/* ----------------------------------------------------------------- */


static void MergePartitions(int           a,
                            int           b,
                            TRAININGSET*  TS,
                            PARTITIONING* P)
{
  if( MaintainPartitions == YES )
    {
    JoinPartitions(TS, P, a, b);
    }
}


/* ----------------------------------------------------------------- */


static void MoveLastVector(int           a,
                           int           b,
                           CODEBOOK*     CB,
                           TRAININGSET*  TS,
                           PARTITIONING* P,
                           NEARESTTYPE   NN[],
                           KURITADATA*   KD,
                           HEAP*         H,
                           KDTREE*       tree)
{
  int last = BookSize(CB) - 1;
  int i, tmp;

  assert( a < b );

  if( b != last )
    {
    if( Value(ShowProgress) >= 5 )
      {
      PrintMessage("Moving:%i to %i (%i,%i,%.0f)\n",
                   last, b, VectorFreq(CB, b), NN[b].Vector,
                   (double)NN[b].VectorDist);
      }
    CopyNode(&Node(CB, last), &Node(CB, b), VectorSize(CB));
    NN[b] = NN[last];
    NN[b].WhoAmI = b;

    /* Move for partition(s) is done in JoinPartitions. */

    for( i = 0; i < BookSize(CB); i++ )
      {
      if( NN[i].Vector == last )
        {
        NN[i].Vector = b;
        }
      }
    }

  switch( Value(MergeMethod) )
    {
    case Kissing:
    case TrivialON3:
      {
      break;
      }
    case Kurita:
      {
      /* Remove elements from the heap */
      for( i = 0; i < BookSize(CB); i++ )
        {
        if( i != a && i != b )
          {
          Heap_remove(&(KD->H), DM_heapindex(KD, b, i), NULL);
          }
        }

      if( b != last )
        {
        /* Update mappings */
        tmp = cb2m(KD, b);
        cb2m(KD, b) = cb2m(KD, last);
        cb2m(KD, last) = tmp;

        m2cb(KD, cb2m(KD, b)) = b;
        m2cb(KD, cb2m(KD, last)) = last;
        }
#if 0
      for( i = 0; i < BookSize(TS); i++ )
        {
        PrintMessage("i:cb m %i %i %i\n", i, cb2m(KD, i), m2cb(KD, i));
        }
#endif
      break;
      }
    case LazyPNN:
      {
      Heap_moveelem(H, NN[b].HeapIndex, &(NN[b]), &(NN[b].HeapIndex));
      break;
      }
    case FastPNN:
      {
      assert( 0 );
      break;
      }
    }
}

/*--------------------------------------------------------------------*/


static void MergeVectors(TRAININGSET*  TS,
                         CODEBOOK*     CB,
                         PARTITIONING* P,
                         NEARESTTYPE   NN[],
                         KURITADATA*   KD,
                         HEAP*         H,
                         KDTREE*       tree,
                         int           a,
                         int           b,
                         llong         Dist)
/* Merge Node(CB,a) and Node(CB,b). */
{
  llong NewError = -1;
  int   i;
  int   recal = 0;

  if( Value(ShowProgress) >= 3 )
    {
    PrintMessage("\nMergeVectors:%4i & %4i "
                 "F=%5i & %5i LD=%.0f & %.0f Dist=%.0f\n",
                 a, b,
                 VectorFreq(CB, a), VectorFreq(CB, b),
                 (double)NN[a].LocalDist,
                 (double)NN[b].LocalDist,
                 (double)Dist);
    }

  /* We assume ... something. */
  assert( 0 <= a );
  assert( a < b );
  assert( b < BookSize(CB) );

  switch( Value(MergeMethod) )
    {
    case Kissing:
    case LazyPNN:
      {
      /* Mark clusters for recalculation. */
      for( i = 0; i < BookSize(CB); i++ )
        {
        if( NN[i].Vector == a || NN[i].Vector == b )
          {
          NN[i].Recalculate = 1;
          recal++;
          }
        }
      NN[a].Recalculate = 2;

      if( Value(ShowProgress) >= 2 )
        {
        PrintMessage("REC: %4i %4i ", recal, BookSize(CB)-recal); fflush(stdout);
        }
      break;
      }
    case TrivialON3:
    case Kurita:
    case FastPNN:
      {
      break;
      }
    }

  /* New error value for merged partition */
  switch( Value(MergeErrorType) )
    {
    case Equitz:
    case TrainingSetIncrease:
      {
      NewError = NN[a].LocalDist + NN[b].LocalDist + Dist;
      break;
      }
    case TrainingSetAmount:
      {
      NewError = Dist;
      break;
      }
    default:
      {
      ErrorMessage("Unknown MergeErrorType=%i", Value(MergeErrorType));
      ExitProcessing( -1 );
      }
    }

  /* New codevector for merged partition */
  CalculateNewCentroid(a, b, CB, TS, P, &Node(CB, a));

  /* Merge two vector clusters from training set */
  MergePartitions(a, b, TS, P);

  /* Ensure that distances from 'a' are recalculated. */ /* Still needed?? */
  NN[a].Vector = a;
  NN[a].VectorDist = MAXLLONG;
  NN[a].LocalDist = NewError;
  
  if( Value(ShowProgress) >= 4 )
    {
    PrintMessage("MergeVectors (joined): a= %4i aF= %5i\n",
                 a, VectorFreq(CB, a));
    }

  /* Move Node(CB, BookSize(CB) - 1) to Node(CB, b). */
  MoveLastVector(a, b, CB, TS, P, NN, KD, H, tree);

  BookSize(CB)--;

  if( Value(PartitionRemapping) )
    {
/*
    LocalRepartitioningGeneral(TS, CB, P, a, EUCLIDEANSQ);
*/
    int   j, nearest;
    llong error;

    for( j = FirstVector(P,a); ! EndOfPartition(j); j = NextVector(P,j) )
      {
      nearest = FindNearestVector(&Node(TS, j), CB, &error, Map(P, j), EUCLIDEANSQ);
      if( nearest != Map(P, j) )
        {
        NN[nearest].Recalculate = 1;      /* Problem... */
        NN[Map(P, j)].Recalculate = 1;    /* Problem... */
        ChangePartition(TS, P, nearest, j);
        }
      }

    }
}


/* ----------------------------------------------------------------- */


static void MergeVectorsInBuckets(TRAININGSET*  TS,
                                  CODEBOOK*     CB,
                                  PARTITIONING* P,
                                  NEARESTTYPE   NN[],
                                  KURITADATA*   KD,
                                  HEAP*         H,
                                  KDTREE*       tree,
                                  int           ResultCBSize)
{
  int i;
  int mergebuckets = (int)((double)tree->inorder *
                           (double)Value(MergePercent)/100.0 + 0.5);
  int a, b;
  int last, bindex;

  mergebuckets = (mergebuckets > 0 ? mergebuckets : 1);
  if( BookSize(CB) - mergebuckets < ResultCBSize )
    {
    mergebuckets = BookSize(CB) - ResultCBSize;
    }
  for( i = 0; i < mergebuckets; i++ )
    {
    /* Merge two nearest clusters in the bucket. */
    if( tree->bucket[tree->bucketorder[i]].cand1 <
        tree->bucket[tree->bucketorder[i]].cand2 )
      {
      a = tree->bucket[tree->bucketorder[i]].cand1;
      b = tree->bucket[tree->bucketorder[i]].cand2;
      }
    else
      {
      a = tree->bucket[tree->bucketorder[i]].cand2;
      b = tree->bucket[tree->bucketorder[i]].cand1;
      }

    assert( 0 <= a );
    assert( a < b );
    assert( b < BookSize(CB) );

    /* New error value for merged partition */
    NN[a].LocalDist = NN[a].LocalDist +
                      NN[b].LocalDist +
                      tree->bucket[tree->bucketorder[i]].dist;

    /* New codevector for merged partition */
    CalculateNewCentroid(a, b, CB, TS, P, &Node(CB, a));

    /* Merge two vector clusters from training set */
    MergePartitions(a, b, TS, P);

    tree->toberemoved[b] = YES;
    }

  for( last = BookSize(CB)-1; tree->toberemoved[last] == YES; last-- )
    {
    }
  for( i = 0; i < mergebuckets; i++ )
    {
    if( tree->bucket[tree->bucketorder[i]].cand1 <
        tree->bucket[tree->bucketorder[i]].cand2 )
      {
      b = tree->bucket[tree->bucketorder[i]].cand2;
      }
    else
      {
      b = tree->bucket[tree->bucketorder[i]].cand1;
      }

    RemoveVectorFromBucket(tree, b);
    if( b < last )
      {
      bindex = tree->map[last];
      RemoveVectorFromBucket(tree, last);
      AddVectorToBucket(tree, b, bindex);
      CopyNode(&Node(CB, last), &Node(CB, b), VectorSize(CB));
      tree->toberemoved[b] = NO;
      NN[b] = NN[last];
      while( tree->toberemoved[--last] == YES )
        {
        }
      }
    }

  BalanceKDTree(CB, tree);

  BookSize(CB) -= mergebuckets;
}


/* ----------------------------------------------------------------- */


static void Merge(TRAININGSET*  TS,
                  CODEBOOK*     CB,
                  PARTITIONING* P,
                  NEARESTTYPE   NN[],
                  KURITADATA*   KD,
                  HEAP*         H,
                  KDTREE*       tree,
                  int           a,
                  int           b,
                  llong         Dist,
                  int           ResultCBSize)
{
  if( Value(MergeMethod) == FastPNN )
    {
    /* Merge two nearest vectors in the buckets. */
    MergeVectorsInBuckets(TS, CB, P, NN, KD, H, tree, ResultCBSize);
    }
  else
    {
    /* Merge first two nearest clusters. */
    MergeVectors(TS, CB, P, NN, KD, H, tree, a, b, Dist);
    }
}


/* ----------------------------------------------------------------- */


static llong MergeDistortion(TRAININGSET*  TS,
                             CODEBOOK*     CB,
                             PARTITIONING* P,
                             NEARESTTYPE   NN[],
                             int           p1,
                             int           p2)
{
  llong    distortion  = 0LL; /* Return value */
  llong    error       = 0LL;
  double   coefficient = 0.0;
  BOOKNODE cent;

  switch( Value(MergeErrorType) )
    {
    case Equitz:
      {
      if( VectorFreq(CB, p1) == 0 || VectorFreq(CB, p2) == 0 )
        {
        distortion = 0LL;
        }
      else
        {
        error = VectorDistance(Vector(CB, p1), Vector(CB, p2),
                               VectorSize(CB), MAXLLONG, EUCLIDEANSQ);
        coefficient = ((double)VectorFreq(CB,p1) * (double)VectorFreq(CB,p2)) /
                      ((double)VectorFreq(CB,p1) + (double)VectorFreq(CB,p2));
        distortion = round(coefficient * (double)error);
        }
      break;
      }
    case TrainingSetAmount:
      {
      /* PartitionError has to calculate with MAXLLONG,
         because of the structure of CalculateSomeDistances. */

      cent = CreateEmptyNode(VectorSize(CB));

      CentroidOf2Partitions(P, p1, p2, VectorSize(CB), &cent);
      error = 0;
      error += PartitionError(TS, P, p1, &cent, MAXLLONG);
      error += PartitionError(TS, P, p2, &cent, MAXLLONG);

      FreeNode(cent);
      distortion = error;
      break;
      }
    case TrainingSetIncrease:
      {
      /* PartitionError has to calculate with MAXLLONG,
         because of the structure of CalculateSomeDistances. */

      assert( NN != NULL );

      cent = CreateEmptyNode(VectorSize(CB));

      CentroidOf2Partitions(P, p1, p2, VectorSize(CB), &cent);
      error = 0;
      error -= (NN[p1].LocalDist + NN[p2].LocalDist);
      error += PartitionError(TS, P, p1, &cent, MAXLLONG);
      error += PartitionError(TS, P, p2, &cent, MAXLLONG);

      FreeNode(cent);
      distortion = error;
      break;
      }
    default:
      {
      ErrorMessage("Unknown MergeErrorType=%i", Value(MergeErrorType));
      ExitProcessing( -1 );
      }
    }

  if( Value(ShowProgress) >= 5 )
    {
    PrintMessage("MergeDist.:p=%4i,%4i F=%5i,%5i e=%7.0f c=%9.4f di=%7.0f\n",
                 p1, p2, VectorFreq(CB,p1), VectorFreq(CB,p2),
                 (double)error, coefficient, (double)distortion);
    }

  return( distortion );
}


/* ================================================================= */


static void IterateGLA(TRAININGSET*    TS,
                       CODEBOOK*       CB,
                       SASchedule*     SAS,
                       PARTITIONING*   P)
{
  int   i = 0;
  llong e, error = MAXLLONG;

  if( Value(GLAIterations) > 0 )
    {
    e = ErrorOfPartitions(TS, CB, P);
    while( i < Value(GLAIterations) && e < error )
      {
      error = e;
      GenerateOptimalPartitioning(TS, CB, P);
      GenerateOptimalCodebook(TS, CB, P);
      FillEmptyPartitions(TS, CB, P);
      e = ErrorOfPartitions(TS, CB, P);
      i++;
      }

    if( Value(ShowProgress) >= 3 )
      {
      PrintMessage("GLA:%3i ", i);
      }
    }
}


/* ====================== KURITA'S METHOD ========================== */


static void DM_print(KURITADATA* KD, CODEBOOK* CB)
{
  int   i, j;

  PrintMessage("##\n");
  for( i = 0; i < BookSize(CB) - 1; i++ )
    {
    for( j = i + 1; j < BookSize(CB); j++ )
      {
      PrintMessage("%2i,%2i D%5.0f H%2i \n",
                   i, j, (double)DM_dist(KD, i, j), DM_heapindex(KD, i, j));
      }
    PrintMessage("\n");
    }
}

/* ----------------------------------------------------------------- */

static void CalculateKuritaDistances(TRAININGSET*  TS,
                                     CODEBOOK*     CB,
                                     PARTITIONING* P,
                                     KURITADATA*   KD,
                                     int*          a,
                                     int*          b,
                                     llong*        Distab)
{
  int   i, j;
  llong Dist;
  DISTELEM* cell;

  for( i = 0; i < BookSize(CB) - 1; i++ )
    {
    for( j = i + 1; j < BookSize(CB); j++ )
      {
      Dist = MergeDistortion(TS, CB, P, NULL, i, j);
      DM_distset(KD, i, j, Dist);
      Heap_insert(&(KD->H),
                  DM_elem(KD, i, j),
                  NULL,
                  DM_heapindexaddr(KD, i, j));
/* PrintMessage("%2i,%2i D%5.0f H%2i \n", i, j, (double)Dist, DM_heapindex(KD, i, j)); */
      }
/* >PrintMessage("\n"); */
    }

/* DM_print(KD, CB); */

  cell = Heap_removeroot(&(KD->H), NULL);
  *a = m2cb(KD, cell->from);
  *b = m2cb(KD, cell->to);
  *Distab = cell->dist;

/* PrintMessage("f,t=%i,%i -> a,b=%i,%i\n", cell->from, cell->to, *a, *b); */
/* DM_print(KD, CB); */

  if( *a > *b )
    {
    swapint(*a, *b);
    }

  /* We assume ... something. */
  assert( 0 <= *a );
  assert( *a < *b );
  assert( *b < BookSize(CB) );
}


/* ----------------------------------------------------------------- */


static void UpdateKuritaDistances(TRAININGSET*  TS,
                                  CODEBOOK*     CB,
                                  PARTITIONING* P,
                                  KURITADATA*   KD,
                                  int*          a,
                                  int*          b,
                                  llong*        Distab)
{
  int i;
  llong Dist;
  DISTELEM* cell;

  assert( *a < *b );

/* DM_print(KD, CB); */

  /* Update distances from 'a' */
  for( i = 0; i < BookSize(CB); i++ )
    {
    if( i != *a )
      {
      Dist = MergeDistortion(TS, CB, P, NULL, *a, i);
      DM_distset(KD, *a, i, Dist);
      Heap_update(&(KD->H), DM_heapindex(KD, *a, i), NULL);
      }
    }

/* DM_print(KD, CB); */

  cell = Heap_removeroot(&(KD->H), NULL);
  *a = m2cb(KD, cell->from);
  *b = m2cb(KD, cell->to);
  *Distab = cell->dist;

/*   PrintMessage("f,t=%i,%i -> a,b=%i,%i\n", cell->from, cell->to, *a, *b); */
/* DM_print(KD, CB); */

  if( *a > *b )
    {
    swapint(*a, *b);
    }

  /* We assume ... something. */
  assert( 0 <= *a );
  assert( *a < *b );
  assert( *b < BookSize(CB) );
}


/* ================================================================= */


static void UpdateDistance(TRAININGSET*  TS,
                           CODEBOOK*     CB,
                           PARTITIONING* P,
                           NEARESTTYPE   NN[],
                           int           a)
{
  int   i;
  llong Dist;

  NN[a].VectorDist = MAXLLONG;
  for( i = 0; i < BookSize(CB); i++ )
    {
    if( i != a )
      {
      Dist = MergeDistortion(TS, CB, P, NN, a, i);
      if( Dist < NN[a].VectorDist )
        {
        NN[a].Vector     = i;
        NN[a].VectorDist = Dist;
        }
      }
    }
  NN[a].Recalculate = 0;

  assert( 0 <= NN[a].Vector );
  assert( NN[a].Vector < BookSize(CB) );
  assert( 0 <= NN[a].VectorDist );
}


/* ----------------------------------------------------------------- */


static void InsertDistancesToHeap(CODEBOOK* CB, HEAP* H, NEARESTTYPE NN[])
{
  int i;

  for( i = 0; i < BookSize(CB); i++ )
    {
    Heap_insert(H, &(NN[i]), NULL, &(NN[i].HeapIndex));
    }
}


/* ----------------------------------------------------------------- */


static void UpdateLazyDistances(TRAININGSET*  TS,
                                CODEBOOK*     CB,
                                PARTITIONING* P,
                                NEARESTTYPE   NN[],
                                HEAP*         H,
                                int*          a,
                                int*          b,
                                llong*        Distab)
{
  NEARESTTYPE* pair;
  int          nrecalculated = 0;

  assert( *a < *b );
  assert( NN[*a].Recalculate != 0 );

/* PrintMessage("H=%i CB=%i ", Heap_size(H), BookSize(CB)); */

  do
    {
    UpdateDistance(TS, CB, P, NN, *a);

    Heap_insert(H, &(NN[*a]), NULL, &(NN[*a].HeapIndex));
    pair = Heap_removeroot(H, NULL);

    *a = pair->WhoAmI;
    nrecalculated++;
    } while( NN[*a].Recalculate != 0 );

  *b = pair->Vector;
  *Distab = pair->VectorDist;

  Heap_remove(H, NN[*b].HeapIndex, NULL);

  if( *a > *b )
    {
    swapint(*a, *b);
    }

/*
PrintMessage("W,V,a,b=%i,%i,%i,%i -> %5.0f CB%i\n",
       pair->WhoAmI, pair->Vector, *a, *b, (double)*Distab, BookSize(CB));
*/
  if( Value(ShowProgress) >= 3 )
    {
    PrintMessage("UpdLazy:rec=%3i size=%4i Distab=%7.0f a=%4i\n",
           nrecalculated, BookSize(CB), (double)*Distab, *a);
    }

  /* We assume ... something. */
  assert( 0 <= *a );
  assert( *a < *b );
  assert( *b < BookSize(CB) );
}


/* ====================  Fast PNN by Equitz  ======================= */


/* ----------------------------------------------------------------- */


static void SearchPairInBucket(TRAININGSET*  TS,
                               CODEBOOK*     CB,
                               PARTITIONING* P,
                               NEARESTTYPE   NN[],
                               KDTREE*       tree,
                               int           bindex)
{
  int       i, j;
  llong     Dist;
  KDBUCKET* buc = &tree->bucket[bindex];
  int*      next = tree->next;

  assert( 1 <= buc->size );

  buc->dist = MAXLLONG;
  for( i = buc->first; next[i] != -1; i = next[i] )
    {
    for( j = next[i]; j != -1; j = next[j] )
      {
      Dist = MergeDistortion(TS, CB, P, NN, i, j);
      if( Dist < buc->dist )
        {
        buc->cand1 = i;
        buc->cand2 = j;
        buc->dist = Dist;
        }
      }
    }
}


/* ----------------------------------------------------------------- */


static void CalculateBucketDistances(TRAININGSET*  TS,
                                     CODEBOOK*     CB,
                                     PARTITIONING* P,
                                     NEARESTTYPE   NN[],
                                     KDTREE*       tree)
{
  STACK*  S = S_make();
  KDNODE* node;

  S_push(S, tree->root);
  while( ! S_empty(S) )
    {
    node = S_pop(S);
    if( node->bucket == -1 )
      {
      assert( node->right != NULL );
      assert( node->left != NULL );
      S_push(S, node->right);
      S_push(S, node->left);
      }
    else
      {
      SearchPairInBucket(TS, CB, P, NN, tree, node->bucket);
      }
    }

  S_free(S);
}


/* ----------------------------------------------------------------- */


static void SelectMergeBuckets(TRAININGSET*  TS,
                               CODEBOOK*     CB,
                               PARTITIONING* P,
                               NEARESTTYPE   NN[],
                               KDTREE*       tree)
{
  STACK*  S = S_make();
  KDNODE* node;
  int     i;

  tree->inorder = 0;
  S_push(S, tree->root);
  while( ! S_empty(S) )
    {
    node = S_pop(S);
    if( node->bucket == -1 )
      {
      assert( node->right != NULL );
      assert( node->left != NULL );
      S_push(S, node->right);
      S_push(S, node->left);
      }
    else
      {
      for( i = tree->inorder-1;
           0 <= i &&
           tree->bucket[node->bucket].dist <
           tree->bucket[tree->bucketorder[i]].dist;
           i-- )
        {
        tree->bucketorder[i+1] = tree->bucketorder[i];
        }
      tree->bucketorder[i+1] = node->bucket;
      tree->inorder++;
      }
    }

  S_free(S);
}


/* ================================================================= */


static llong CalculateAllDistances(TRAININGSET*  TS,
                                   CODEBOOK*     CB,
                                   PARTITIONING* P,
                                   NEARESTTYPE   NN[],
                                   int*          a,
                                   int*          b,
                                   llong*        Distab)
{
  int   i, j;
  llong Dist;
  llong totalerror = 0LL;

  for( i = 0; i < BookSize(CB); i++ )
    {
    NN[i].VectorDist = MAXLLONG;
    NN[i].Vector     = i;
    }

  *Distab = MAXLLONG;
  for( i = 0; i < BookSize(CB) - 1; i++ )
    {
    if( Value(ShowProgress) >= 4 )
      {
      PrintMessage("\nCalcAll:%4i: Dist=%7.0f Vect=%4i Distab=%7.0f a=%4i ",
                   i, (double)NN[i].VectorDist, NN[i].Vector, (double)*Distab, *a);
      }

    for( j = i + 1; j < BookSize(CB); j++ )
      {
      Dist = MergeDistortion(TS, CB, P, NN, i, j);
      if( Dist < NN[i].VectorDist )
        {
        NN[i].Vector     = j;
        NN[i].VectorDist = Dist;
        }
      if( Dist < NN[j].VectorDist )
        {
        NN[j].Vector     = i;
        NN[j].VectorDist = Dist;
        }
      }

    if( NN[i].VectorDist < *Distab )
      {
      *Distab = NN[i].VectorDist;
      *a = i;
      *b = NN[i].Vector;
      }

    totalerror += NN[i].LocalDist;

    if( Value(ShowProgress) >= 4 )
      {
      PrintMessage("\nCalcAll:%4i: Dist=%7.0f Vect=%4i Distab=%7.0f a=%4i ",
             i, (double)NN[i].VectorDist, NN[i].Vector, (double)*Distab, *a);
      }
    }

  if( *a > *b )
    {
    swapint(*a, *b);
    }

  if( Value(ShowProgress) >= 4 )
    {
    PrintMessage(" done\n");
    }

  assert( totalerror >= 0LL );

  return( totalerror );
}


/*--------------------------------------------------------------------*/


static llong CalculateSomeDistances(TRAININGSET*  TS,
                                    CODEBOOK*     CB,
                                    PARTITIONING* P,
                                    NEARESTTYPE   NN[],
                                    int*          a,
                                    int*          b,
                                    llong*        Distab)
{
  int    i, j;
  llong  Dist;
  int    nrecalculated = 0;
  llong  totalerror = 0LL;

  if( Value(GLAIterations) > 0 || Value(PartitionRemapping) )
    {
    CalculatePartitionErrors(TS, CB, P, NN);
    }

  *Distab = MAXLLONG;
  *a = 0; /* ?? */
  *b = 0; /* ?? */

  for( i = 0; i < BookSize(CB); i++ )
    {
    if( NN[i].Recalculate )
      {
      nrecalculated++;
      NN[i].VectorDist = MAXLLONG;
      for( j = 0; j < BookSize(CB); j++ )
        {
        if( i != j )
          {
          Dist = MergeDistortion(TS, CB, P, NN, i, j);
          if( Dist < NN[i].VectorDist )
            {
            NN[i].Vector     = j;
            NN[i].VectorDist = Dist;
            }
          if( Dist < NN[j].VectorDist )  /* Rare? */
            {
            NN[j].Vector     = i;
            NN[j].VectorDist = Dist;
            }
          }
        NN[i].Recalculate = 0;
        }
      if( Value(ShowProgress) >= 4 )
        {
        PrintMessage("CalcSome:rec=%3i i=%4i Dist%7.0f "
                     "Vect=%4i Distab=%7.0f a=%4i\n",
                     nrecalculated, i, (double)NN[i].VectorDist,
                     NN[i].Vector, (double)*Distab, *a);
        }
      }
    if( NN[i].VectorDist < *Distab )
      {
      *Distab = NN[i].VectorDist;
      *a = i;
      *b = NN[i].Vector;
      }

    assert( NN[i].LocalDist >= 0LL );

    totalerror += NN[i].LocalDist;
    }

  if( *a > *b )
    {
    swapint(*a, *b);
    }

  if( Value(ShowProgress) >= 3 )
    {
    PrintMessage("CalcSome:rec=%3i size=%4i Distab=%7.0f a=%4i\n",
                 nrecalculated, BookSize(CB), (double)*Distab, *a);
    }

  assert( totalerror >= 0LL );

  return( totalerror );
}


/*--------------------------------------------------------------------*/


void PairwiseNearestNeighbour(TRAININGSET*  TS,
                              CODEBOOK*     CB,
                              SASchedule*   SAS,
                              PARTITIONING* OrigP,
                              int           ResultCBSize)
{
  int           a, b;
  llong         Distab;
  PARTITIONING* P;
  NEARESTTYPE*  NN;
  KURITADATA    KD;
  HEAP          H;
  KDTREE        tree;
  llong         Perror = 0LL;

  if( BookSize(CB) <= ResultCBSize ) /* Nothing to do. */
    {
    return;
    }
  InitializePNN(TS, CB, OrigP, &P, &NN, &KD, &H, &tree);
  GenerateInitialPartition(TS, CB, P, NN);

  /* Calculate initial distances. */
  switch( Value(MergeMethod) )
    {
    case Kissing:
    case TrivialON3:
      {
      CalculateAllDistances(TS, CB, P, NN, &a, &b, &Distab);
      break;
      }
    case Kurita:
      {
      CalculateKuritaDistances(TS, CB, P, &KD, &a, &b, &Distab);
      break;
      }
    case LazyPNN:
      {
      CalculateAllDistances(TS, CB, P, NN, &a, &b, &Distab);
      InsertDistancesToHeap(CB, &H, NN);
      Heap_remove(&H, NN[a].HeapIndex, NULL);
      Heap_remove(&H, NN[b].HeapIndex, NULL);
      break;
      }
    case FastPNN:
      {
      MakeKDTree(CB, &tree);
      CalculateBucketDistances(TS, CB, P, NN, &tree);
      SelectMergeBuckets(TS, CB, P, NN, &tree);
      break;
      }
    default:
      break;
    }

  /* Merge clusters. */
  Merge(TS, CB, P, NN, &KD, &H, &tree, a, b, Distab, ResultCBSize);
  IterateGLA(TS, CB, SAS, P);

  /* Main PNN loop */
  while( BookSize(CB) > ResultCBSize )
    {
    /* Select two nearest clusters (a,b). */
    switch( Value(MergeMethod) )
      {
      case Kissing:
        {
        Perror = CalculateSomeDistances(TS, CB, P, NN, &a, &b, &Distab);
        break;
        }
      case TrivialON3:
        {
        Perror = CalculateAllDistances(TS, CB, P, NN, &a, &b, &Distab);
        break;
        }
      case Kurita:
        {
        UpdateKuritaDistances(TS, CB, P, &KD, &a, &b, &Distab);
        Perror = 0LL;
        break;
        }
      case LazyPNN:
        {
        UpdateLazyDistances(TS, CB, P, NN, &H, &a, &b, &Distab);
        Perror = 0LL;
        break;
        }
      case FastPNN:
        {
        CalculateBucketDistances(TS, CB, P, NN, &tree);
        SelectMergeBuckets(TS, CB, P, NN, &tree);
        break;
        }
      default:
        break;
      }
    Merge(TS, CB, P, NN, &KD, &H, &tree, a, b, Distab, ResultCBSize);
    IterateGLA(TS, CB, SAS, P);

    if( Value(ShowProgress) >= 2 )
      {
      PrintMessage("PNN CB:%5u %9.4f\n",
                   BookSize(CB), SE2MSE(TS, (double)Perror));
      fflush(stdout);
      }
    }

  ShutDownPNN(TS, CB, OrigP, P, NN, &KD, &H, &tree);
}


