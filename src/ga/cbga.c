/*-------------------------------------------------------------------*/
/* CBGA.C         Pasi Fränti & Olli Virmajoki                       */
/*                                                                   */
/* - Genetic Algorithm - "Next generation"                           */
/* - Current version supports:                                       */
/*    RANDOM, CENTROID DISTANCE, PNN-old, PNN-new, ADOPTION SWAP,    */
/*    IS, Adaptive IS                                                */
/*-------------------------------------------------------------------*/

#define ProgName       "CBGA"
#define VersionNumber  "Version 0.12"
#define LastUpdated    "09.08.04"

/* ----------------------------------------------------------------- */

#include <assert.h>
#include <float.h>
#include <time.h>
#include <values.h>

#define  FACTFILE  "cbga.fac"

#include "parametr.c"
#include "cb.h"
#include "file.h"
#include "memctrl.h"
#include "random.h"
#include "pnn.h"
#include "sa.h"
#include "solution.h"
#include "sortcb.h"
#include "textfile.h"

/*-------------------------- DOS special -----------------------------*/

#if defined(__MSDOS__)
#include <conio.h>
#endif

/*--------------------------- Constants ------------------------------*/

#define  MaxGenerations  100
#define  CBX(i)          ((i%2==0)? &S1->CB : &S2->CB)

/*------------------------ Global variables --------------------------*/

int      CrossSetSize;
int      Survivors;

/*---------------------------  M i s c	------------------------------*/

#define  fpositive(a)   ((a) > FLT_EPSILON)
#define  swap(a,b)      { int t = a; a = b; b = t; }
#define  round(a)		((llong) (a+0.5))

#if defined(max) || defined(min)
#undef max
#undef min
#endif
#define  max(a,b) ((a) > (b) ? (a) : (b))
#define  min(a,b) ((a) < (b) ? (a) : (b))
#define  ScaleBetween(a,b,c) max((a), min((b), (c)))

#define       fpositive(a)  ((a) > FLT_EPSILON)
#define       NumberOfIterations (Value(NumberOfGenerations) ? \
                                  Value(NumberOfGenerations) : \
                                  Max(NumberOfGenerations))
#define       NoImprovement (!fpositive(PrevError-Error[0]) \
                             && Value(NumberOfGenerations)==0 \
                             && g>2)
static void RemoveEmptyPartitions(TRAININGSET*  TS,
                                  CODEBOOK*     CB,
                                  PARTITIONING* P);


/* =================  T I M E   F U N C T I O N S  ================= */


static long SetWatch(long* start)
{
  time_t tmp;
  return( *start = (long)time(&tmp) );
}


/* ----------------------------------------------------------------- */


static void PrintTime(long start)
{
  time_t tmp;
  long t = (long)time(&tmp) - start;
  long h,m,s;

  s = t % 60;
  t = t / 60;
  m = t % 60;
  h = t / 60;

  switch( Value(QuietLevel) )
    {
    case 0:  break;
    case 1:
    case 2:
    case 3:
    case 4:  printf("  %li:%02li:%02li\n", h, m, s);
             break;
    default: printf("Time: %li:%02li:%02li\n", h, m, s);
             break;
    }
}


/*============================  P R I N T  ===============================*/


void PrintInfo(void)
{
  printf("\n%s %s %s\n", ProgName, VersionNumber, LastUpdated);
  printf("GA - next generation generates a codebook using Genetic Algorithm.\n");
  printf("Usage: %s [%coption] <training set> <codebook>\n", ProgName, OPTION_SYMBOL);
  printf("For example: %s bridge tmp\n", ProgName);
  printf("\n  Options:\n");
  PrintOptions();
  printf("\n%s %s %s\n", ProgName, VersionNumber, LastUpdated);
}


/*-------------------------------------------------------------------*/


static void PrintOperatingInfo(char* TSName, char* CBName)
{
  if( Value(QuietLevel) >= 2 )
    {
    printf("\n%s %s %s\n", ProgName, VersionNumber, LastUpdated);
    printf("Training Set:     %s\n", TSName);
    printf("Codebook:         %s\n", CBName);
    printf("\n");
    PrintSelectedOptions();
    printf("\n");
    }
}


/*-------------------------------------------------------------------*/


static void PrintProgress(CODEBOOK*   CB,
                          SASchedule* SAS,
                          int         gen,
                          double      PrevError,
                          double*     Error,
                          long        watch)
{
  #define COLUMNS (Value(GenerationSize)<6 ? Value(GenerationSize) : 6)
  double Diff,NewError;
  int i;

  PrevError = PrintableError(PrevError, CB);
  NewError  = PrintableError(Error[0],  CB);
  Diff = (PrevError - NewError);

  switch(  Value(QuietLevel) )
	{
	case 0:  return;
    case 1:  return;
    case 2:  printf("%3i: %9.4f \n", gen, NewError); break;
    case 3:  printf("%3i: %9.4f  ", gen, NewError); PrintTime(watch); break;
    case 4:  printf("%3i: ",gen);
			 for (i=0; i<COLUMNS; i++)
                 printf("%.4f ", PrintableError(Error[i],CB));
             PrintTime(watch);
			 break;
    default: printf("Iter=%3i   Error=%9.4f  Change=%-9.4f ",
			 gen, NewError, Diff);
			 if( SASInUse(SAS) )
			   {
			   printf(" Temp=%8.4f ", SAS->CurrentTemperature);
			   }
            PrintTime(watch);
			break;
	}
  fflush(stdout);
}


/*-------------------------------------------------------------------*/


static void PrintResult(long watch, double error, int iterations)
{
  switch( Value(QuietLevel) )
	 {
     case 0:
     case 1:
     case 2:  printf("%9.4f  (%i)", error, iterations); break;
	 default: {
                printf("Distortion: %-9.4f\n", error);
                printf("Iterations: %-3i\n", iterations);
			  } break;
	 }
  PrintTime(watch);
}


/* ================================================================= */


static void Initialize(TRAININGSET* TS,
                       SOLUTION*    Sold[],
                       SOLUTION*    Snew[],
                       SASchedule*  SAS)
{
  int i;
  int s = (Value(CrossingMethod)==8 ? MaxGenerations : Value(GenerationSize));

  if( Value(QuietLevel)>=4 )
      printf("Initializing for %i solutions.\n",Value(GenerationSize));
  if( Value(CodebookSize)>=BookSize(TS) )
      {
      printf("Codebook size (%i) > training set size (%i).\n",
              Value(CodebookSize),BookSize(TS));
      exit(-1);
      }
  for (i=0; i<s; i++)
     {
     Sold[i] = (SOLUTION*) allocate( sizeof(SOLUTION) );
     Snew[i] = (SOLUTION*) allocate( sizeof(SOLUTION) );
     CreateNewSolution(TS, Sold[i], Value(CodebookSize));
     CreateNewSolution(TS, Snew[i], 2*Value(CodebookSize));
     }

  InitializeSASchedule(Value(TemperatureDistribution),
                       Value(TemperatureFunction),
                       Value(TemperatureChange),
               (double)Value(TemperatureInitial) * (double)(TS->MaxValue) / 100.0,
                      (Value(VectorRandomizing) == RandomCentroid ||
                       Value(VectorRandomizing) == RandomBoth),
                      (Value(VectorRandomizing) == RandomTrainingSet ||
                       Value(VectorRandomizing) == RandomBoth),
                       SAS);
}


/* ----------------------------------------------------------------- */


static void ShutDown(TRAININGSET* TS, SOLUTION** Sold, SOLUTION** Snew)
{
  int i;

  FreeCodebook(TS);
  for (i=0; i<Value(GenerationSize); i++)
     {
	 FreeSolution(Sold[i]);
	 FreeSolution(Snew[i]);
     deallocate(Sold[i]);
     deallocate(Snew[i]);
     }
}


/*======================  SOLUTION HANDLING  ========================*/


static int FreqToIndex(TRAININGSET* TS, int sum)

{
  int i=0,l=0;

  while (l<sum)  l += VectorFreq(TS, i++);
  return (i-1);
}


/*-------------------------------------------------------------------*/


static void GenerateRandomVector(TRAININGSET* TS, CODEBOOK* CB, int vec)
{
  int     i;
  YESNO   found=NO;
  int     t;
  int     index;

  while(!found)
     {
	 found=YES;
     if( Value(RandomVector)==FREQ_WEIGHTED )
		{
        t=irand(1,TotalFreq(TS));
        index=FreqToIndex(TS, t);
		}
     else
		{
        index = irand(1,BookSize(TS)-1);
		}
	 for (i=0; i<BookSize(CB); i++)
		if (EqualVectors(Vector(TS,index), Vector(CB,i), VectorSize(CB)))
           {
		   found=NO; break;
           }
     }
  CopyVector(Vector(TS,index), Vector(CB,vec), VectorSize(CB));
  VectorFreq(CB,vec)=1;
}


/*-------------------------------------------------------------------*/


static void GenerateRandomCodebook(TRAININGSET* TS, CODEBOOK* CB)
{
  int i;

  for (i=0; i<BookSize(CB); i++)
      {
      GenerateRandomVector(TS,CB,i);
      }
}


/*-------------------------------------------------------------------*/


static void GenerateRandomPartitioning(TRAININGSET* TS, PARTITIONING* P)
{
  int i;
  int newcluster;

  /* Step 1: Put first M training vectors to different clusters */
  for (i=0; i<PartitionCount(P); i++)
	 {
     ChangePartition(TS, P, i, i);
	 }

  /* Step 2: Make random swaps between cluster */
  for (i=0; i<BookSize(TS); i++)
	 {
	 /* Change mapping only if the old partition remains non-empty */
	 if( UniqueVectors(P,Map(P,i)) > 1 )
		{
        newcluster = irand(0,PartitionCount(P)-1);
		ChangePartition(TS, P, newcluster, i);
		}
	 }
}


/*-------------------------------------------------------------------*/


static void AddNoiseToCodebook(CODEBOOK* CB, SASchedule* SAS)
{
  int i;

  if( !SASUseToCB(SAS) ) return;
  if( Value(QuietLevel)>=5 )
      printf("--> Randomizing codebook.\n");

  for( i=0; i<BookSize(CB); i++ )
    {
	RandomizeVectorBySA(SAS, &Node(CB, i), &Node(CB, i),
						VectorSize(CB), CB->MaxValue);
	}
}


/*-------------------------------------------------------------------*/


static void OutputCodebook(char* CBName, CODEBOOK* CB, TRAININGSET* TS)
{
  int ndups = DuplicatesInCodebook(CB);

  if( ndups>0 && Value(QuietLevel)>=2 )
    {
    printf("WARNING: Final CB contains %i duplicates.\n", ndups);
    }
  SortCodebook(CB, DATA_DESCENDING);
  if(TS->InputFormat==TXT) {
      SaveCB2TXT(CB, CBName, TS->MinMax, NULL, NULL);
  }
  else {
      WriteCodebook(CBName, CB, Value(OverWrite));
  }

}


/*--------------------------------------------------------------------*/


static void SortSolutions(TRAININGSET* TS, SOLUTION* Ss[], double Error[])
{
  int       i,j;
  SOLUTION* S;
  double    e;

  for (i=0; i<Value(GenerationSize); i++)
     for (j=0; j<i; j++)
        if (Error[i]<Error[j])
            {
            S=Ss[i]; Ss[i]=Ss[j]; Ss[j]=S;
            e=Error[i]; Error[i]=Error[j]; Error[j]=e;
            }
}


/*--------------------------------------------------------------------*/


static void CopySolutions(SOLUTION* Source[], SOLUTION* Dest[])
{
  int i;

  for (i=0; i<Value(GenerationSize); i++)
     {
     CopySolution(Source[i], Dest[i]);
     }
}


/*--------------------------------------------------------------------*/


void SortByCentroidDistance(TRAININGSET* TS, CODEBOOK* CB)
{
  int		 i,j;
  llong*	 CBdist;
  llong 	 t;
  VECTORTYPE v;
  VECTORTYPE centroid;

  centroid = CreateEmptyVector(VectorSize(CB));
  CBdist   = (llong*) allocate(BookSize(CB) * sizeof(llong));

  /* calculate distances between code vectors and TS centroid */
  CodebookCentroid(TS, centroid);
  for (i=0; i<BookSize(CB); i++)
	  CBdist[i]=VectorDist(Vector(CB,i), centroid, VectorSize(TS));

  /* sort codebook according to the distances */
  for (i=0; i<BookSize(CB); i++)
	  for (j=0; j<i; j++)
		if (CBdist[j]>CBdist[i])
            {
			t=CBdist[i]; CBdist[i]=CBdist[j]; CBdist[j]=t;
			v=Vector(CB,i); Vector(CB,i)=Vector(CB,j); Vector(CB,j)=v;
            }

  deallocate(CBdist);
  FreeVector(centroid);
}


/*===================== Iterative shrinking =========================*/


int FindSecondNearest(CODEBOOK* TS, CODEBOOK* CB, int i, int tabu1, int tabu2)
/* Finds the nearest excluding candidates tabu1 and tabu2.
   The function can therefore find either 2nd or 3rd nearest. */
{
  int	j;
  int   MinIndex = 0;
  double factor;
  llong MinD = MAXLLONG;
  llong dist;

  for( j=0; j<BookSize(CB); j++ )
	 {
	 if(j==tabu1 || j==tabu2) continue;
     dist = VectorDist(Vector(CB,j), Vector(TS, i), VectorSize(CB));
     factor = ((double)VectorFreq(TS,i) * (double)VectorFreq(CB,j)) /
        ((double)VectorFreq(TS,i) + (double)VectorFreq(CB,j));
     dist = round( factor * (double)dist);
     if( dist<MinD ) { MinD=dist; MinIndex=j; }
     }
  return( MinIndex );
}


/*-------------------------------------------------------------------*/


static void GenerateSecondaryPartitioning(TRAININGSET*	TS,
										  CODEBOOK* 	CB,
										  PARTITIONING* P1,
										  PARTITIONING* P2)
{
  int  i,c;

  for( i=0; i<BookSize(TS); i++ )
	 {
     c = FindSecondNearest( TS, CB, i, -1, Map(P1,i) );
	 ChangePartition(TS, P2, c, i);
     }
}


/*-------------------------------------------------------------------*/


static void GenerateInitialPartitions(TRAININGSET*  TS,
                                      CODEBOOK*     CB,
                                      PARTITIONING* P1,
                                      PARTITIONING* P2)
{
  if( BookSize(TS) == BookSize(CB) )
	 {
	 PutAllInOwnPartition(TS, P1);
	 GenerateSecondaryPartitioning(TS, CB, P1, P2);
	 }
  else /* We have an initial codebook. */
     {
     GenerateOptimalPartitioning(TS, CB, P1);
     GenerateSecondaryPartitioning(TS, CB, P1, P2);
	}
}


/*--------------------------------------------------------------------*/


static void FindClusterToBeRemoved(CODEBOOK* CB,
								   llong	 RemovalCost[],
								   int* 	 cluster)
{
  llong MinDist=MAXLLONG;
  int i;

  for( i=0; i<BookSize(CB); i++ )
	 {
	 if( RemovalCost[i] < MinDist )
		{
		(*cluster) = i;
		MinDist = RemovalCost[i];
		}
	 }
  if( Value(QuietLevel)>=6 ) printf("Cluster %i selected.\n",*cluster);
}


/*-------------------------------------------------------------------*/


static void Init(int* Size, int Freq[], CODEBOOK* CB)
{
    int i;
    for( i=0; i<BookSize(CB); i++)
        Freq[i]=0;
    *Size=0;
}


/*-------------------------------------------------------------------*/


static void Clear(int* Size, int Elem[], int Freq[])
{
   int i;
   for( i=0; i<*Size; i++)
       Freq[Elem[i]]=0;
   *Size=0;
}


/*-------------------------------------------------------------------*/


static void AddVector(llong Counter[],
                      int* vector, int dim, int ti, int freq)
{
    int k;
    for( k = 0; k < dim; k++ )
    {
        Counter[ti*(dim) + k] += vector[k]*freq;
    }
}


/*-------------------------------------------------------------------*/


static void CopVector(llong Counter[],
                      int* vector, int dim, int ti, int freq)
{
    int k;
    for( k = 0; k < dim; k++ )
    {
        Counter[ti*(dim) + k] = vector[k]*freq;
    }
}


/*-------------------------------------------------------------------*/


static void Add(int* Size, int Elem[], int Freq[], llong Counter[],
                int* vector, int dim, int ti, int freq)
{
    if(Freq[ti]==0)
    {
        Elem[*Size]=ti;
        *Size=*Size+1;
        CopVector(Counter, vector, dim, ti, freq);
        Freq[ti]=freq;
    }
    else
    {
        AddVector(Counter, vector, dim, ti, freq);
        Freq[ti]+=freq;
    }
}


/*-------------------------------------------------------------------*/


static void CalcCentroids(int ti, int Freq[], llong Counter[],
                          CODEBOOK* Cents, BOOKNODE* bn)
{
    int k;

    for( k = 0; k<VectorSize(Cents); k++)
    {
        bn->vector[k] = round((double)Counter[ti *VectorSize(Cents) +k]
            /(double)Freq[ti]);
    }
}


/*-------------------------------------------------------------------*/


static void UpdatePrimaryPartition(TRAININGSET* TS,
                                   PARTITIONING* P1, PARTITIONING* P2,
                                   int* Size, int Elem[], int Freq[],
                                   llong Counter[], int obsolete)
{
    int i, temp;
    Clear(Size, Elem, Freq);
    while ( UniqueVectors(P1, obsolete))
    {
        i = FirstVector(P1,obsolete);
        Add(Size, Elem, Freq, Counter, Vector(TS,i), VectorSize(TS),
            Map(P2,i), VectorFreq(TS,i));
        temp = Map(P2,i);
        ChangePartition(TS, P2, Map(P1,i), i);
        ChangePartition(TS, P1, temp, i);
    }
}


/*-------------------------------------------------------------------*/


static void CalcPartitionCentroids(CODEBOOK* CB,
                                   PARTITIONING* P1,
                                   int Size,
                                   int Elem[])
{
    int j;
    for ( j = 0; j <Size; j++)
    {
        PartitionCentroid(P1, Elem[j], &Node(CB, Elem[j]));
    }
}


/*-------------------------------------------------------------------*/


static int MustBeUpdated(PARTITIONING* P1,
                         PARTITIONING* P2,
                         int Freq[],
                         int i,
                         int obsolete,
                         int UpdateStyle)
{
    int must = NO;
    switch( UpdateStyle)
      {
      case MINIMUM:
          {
          if( Map(P2, i) == obsolete)
              must=YES;
          break;
          }
      case STANDARD:
          {
          if( (Map(P2,i) == obsolete) || ( Freq[Map(P1,i)] && Freq[Map(P2,i)] ) )
              must=YES;
          break;
          }
      case EXTENSIVE:
          {
              if( (Map(P2,i) == obsolete) || Freq[Map(P1,i)] || Freq[Map(P2,i)] )
              must=YES;
          break;
          }
      default:
        {
            must=NO;
            break;
        }
      }
    return must;
}


/*-------------------------------------------------------------------*/


static void UpdateSecondaryPartition(TRAININGSET* TS,
									 CODEBOOK* CB,
									 PARTITIONING* P1,
									 PARTITIONING* P2,
                                     int Size,
                                     int Freq[],
                                     int obsolete,
                                     int UpdateStyle,
                                     int Quiet)
{
  int  i,ti;

  for( i=0; i<BookSize(TS); i++ )
  {
      if( MustBeUpdated(P1, P2, Freq, i, obsolete, UpdateStyle))
        {
           ti = FindSecondNearest(TS, CB, i, obsolete, Map(P1,i) );
	       ChangePartition(TS, P2, ti, i);
        }
  }
}


/*-------------------------------------------------------------------*/


static void RemoveObsolete(TRAININGSET* TS,
                           CODEBOOK* CB,
                           PARTITIONING* P1,
                           PARTITIONING* P2,
                           int obsolete)
{
    int last = BookSize(CB)-1;
    CopyNode(&Node(CB, last), &Node(CB, obsolete), VectorSize(CB));
    BookSize(CB)--;
    JoinPartitions(TS, P1, obsolete, last);
    JoinPartitions(TS, P2, obsolete, last);
}


/*-------------------------------------------------------------------*/


static llong CalcSubRemCost(CODEBOOK* CB, int Freq[],
                            CODEBOOK* Cents, int ti)
    /* We merge the sub cluster (ti) and the cluster (ti) */
{
    double factor;
    llong dist;

    dist = VectorDist(Vector(Cents,ti), Vector(CB,ti), VectorSize(CB));
    factor = ((double)Freq[ti] * (double)VectorFreq(CB,ti)) /
        ((double)Freq[ti] + (double)VectorFreq(CB,ti));
    dist = round( factor * (double)dist);

    return (dist);
}


/*-------------------------------------------------------------------*/


static llong CalcRemCost(TRAININGSET* TS,
                         CODEBOOK* CB,
                         int i,
                         int j)
    /* We remove the training vector (i) from the cluster (j) */
{

    llong dist;

    dist = VectorDist(Vector(TS,i), Vector(CB,j), VectorSize(CB)) * (llong)VectorFreq(TS,i);

    return (dist);
}


/*-------------------------------------------------------------------*/


static llong CalcMergeCost(TRAININGSET* TS, CODEBOOK* CB, int i, int j)
    /* We merge the training vector (i) and the code vector (j) */
{
    llong dist;
    double factor;

	dist = VectorDist(Vector(TS,i), Vector(CB,j), VectorSize(TS));
    factor=((double)VectorFreq(TS,i) * (double)VectorFreq(CB,j)) /
         ((double)VectorFreq(TS,i) + (double)VectorFreq(CB,j));
    dist = round( factor*(double)dist);
    return (dist);
}


/*-------------------------------------------------------------------*/


static void CalcRemovalCosts(TRAININGSET* TS, CODEBOOK* CB,
                             PARTITIONING* P1, PARTITIONING* P2,
                             int* Size, int Elem[], int Freq[],
                             llong Counter[], CODEBOOK* Cents,
                             llong RemCost[], int RemovalCalc)
{
   int i, j, l;

   Init(Size, Freq, CB);
   for( j=0; j <BookSize(CB); j++)
   {
     Clear(Size, Elem, Freq);
     RemCost[j]=0;
     for ( i = FirstVector(P1, j); ! EndOfPartition(i); i = NextVector(P1, i) )
     {
        Add(Size, Elem, Freq, Counter, Vector(TS,i), VectorSize(CB),
            Map(P2,i), VectorFreq(TS,i));
        RemCost[j] -= CalcRemCost(TS, CB, i, j);
        if(RemovalCalc==SIMPLE)
        {
            RemCost[j] += CalcMergeCost(TS, CB, i, Map(P2,i));
        }
     }
     if(RemovalCalc==EXACT)
     {
     for( i = 0; i <*Size; i++)
        CalcCentroids(Elem[i], Freq, Counter, Cents, &Node(Cents,Elem[i]));
     for ( i = FirstVector(P1, j); ! EndOfPartition(i); i = NextVector(P1, i) )
     {
        RemCost[j] += CalcRemCost(TS, Cents, i, Map(P2,i));
     }
     for( l = 0; l <*Size; l++)
     {
        RemCost[j] += CalcSubRemCost(CB, Freq, Cents, Elem[l]);
     }
     }/* end RemovalCalc==EXACT */

   }
}


/*-------------------------------------------------------------------*/


static void CreateDataStructures(CODEBOOK* CB, llong**  RemCost, int** Elem,
                               int** Freq, llong** Counter, CODEBOOK* Cents)
{
  *RemCost = (llong*) allocate( BookSize(CB) * sizeof(llong) );
  *Elem = (int*) allocate( BookSize(CB) * sizeof(int) );
  *Freq = (int*) allocate( BookSize(CB) * sizeof(int) );
  *Counter = (llong*) allocate( BookSize(CB) * VectorSize(CB) * sizeof(llong) );
  CreateNewCodebook(Cents, BookSize(CB), CB);
}


/*-------------------------------------------------------------------*/


static void FreeDataStructures(llong*  RemCost,
                               int* Elem,
                               int* Freq,
                               llong* Counter,
                               CODEBOOK* Cents)
{
  deallocate(RemCost);
  deallocate(Elem);
  deallocate(Freq);
  deallocate(Counter);
  FreeCodebook(Cents);
}


/*-------------------------------------------------------------------*/


static void IterativeShrinking( TRAININGSET*   TS,
                                CODEBOOK*      CB,
                                PARTITIONING*  P1,
				                PARTITIONING*  P2,
                                int            ResultCBSize,
                                int            RemovalCalc,
                                int            UpdateStyle,
                                int            Quiet)
{
  llong*    RemCost;
  int       Size;
  int*      Elem;
  int*      Freq;
  llong*    Counter;
  CODEBOOK  Cents;
  int       obsolete;

  CreateDataStructures(CB, &RemCost, &Elem, &Freq, &Counter, &Cents);
  if(Value(PNNParameter)==YES)  RemoveEmptyPartitions(TS, CB, P1);
  GenerateInitialPartitions(TS, CB, P1, P2);

  while( BookSize(CB)>ResultCBSize)
  {
    CalcRemovalCosts(TS, CB, P1, P2, &Size, Elem, Freq, Counter, &Cents, RemCost, RemovalCalc);
    FindClusterToBeRemoved(CB, RemCost, &obsolete);
    UpdatePrimaryPartition(TS, P1, P2, &Size, Elem, Freq, Counter, obsolete);
    CalcPartitionCentroids(CB, P1, Size, Elem);
    UpdateSecondaryPartition(TS, CB, P1, P2, Size, Freq, obsolete, UpdateStyle, Quiet);
    RemoveObsolete(TS, CB, P1, P2, obsolete);
  }

  FreeDataStructures(RemCost, Elem, Freq, Counter, &Cents);
}


/*========================  PNN Routines ============================*/


typedef struct { llong  VectorDist;
                 int    Vector;
                 YESNO  Recalculate;
               } NEARESTTYPE;


/*-------------------------------------------------------------------*/


static llong MergeDistortion(CODEBOOK* CB, int p1, int p2)
{
  llong   dist;
  double  factor;

  if( VectorFreq(CB,p1)==0 || VectorFreq(CB,p2)==0 )
     {
     dist = 0LL;
     }
  else
     {
     dist   = VectorDist(Vector(CB,p1), Vector(CB,p2), VectorSize(CB));
     factor = ((double)VectorFreq(CB,p1) * (double)VectorFreq(CB,p2)) /
              ((double)VectorFreq(CB,p1) + (double)VectorFreq(CB,p2));
     dist   = round( factor * (double) dist);
     }
  return(dist);
}


/*-------------------------------------------------------------------*/


static void FindNearestCluster(CODEBOOK* CB, int a, int* b, llong* dist)
{
  llong Dist;
  int i;

  (*dist)=MAXLLONG; (*b)=0;
  for( i=0; i<BookSize(CB); i++ )
     {
     if( a!=i )
        {
        Dist = MergeDistortion(CB, a, i);
        if( Dist < (*dist) )
           {
           (*dist) = Dist;
           (*b)    = i;
           }
        }
     }
}


/*-------------------------------------------------------------------*/


static void CalculateAllDistances(CODEBOOK* CB, NEARESTTYPE NN[])
{
  int i;

  for( i=0; i<BookSize(CB); i++ )
    {
    FindNearestCluster(CB, i, &(NN[i].Vector), &(NN[i].VectorDist));
    }
}


/*-------------------------------------------------------------------*/


static void FindMergedPair(CODEBOOK* CB, NEARESTTYPE NN[], int* a, int* b)
{
  llong MinDist=MAXLLONG;
  int i;

  for( i=0; i<BookSize(CB); i++ )
    {
    if( NN[i].VectorDist < MinDist )
       {
       (*a) = i;
       (*b) = NN[i].Vector;
       MinDist = NN[i].VectorDist;
       }
    }
  if( (*a)>(*b) ) swap((*a),(*b));
}


/*-------------------------------------------------------------------*/


static void CentroidOfTwoVectors(CODEBOOK* CB, BOOKNODE* cent, int p1, int p2)
{
  int k;
  int freqsum;

  freqsum = VectorFreq(CB,p1) + VectorFreq(CB,p2);

  /* Calculate mean values for centroid */
  for( k=0; k<VectorSize(CB); k++)
    {
	if( freqsum > 0 )
	cent->vector[k] = round( (double)(VectorFreq(CB, p1)*VectorScalar(CB, p1, k) +
                                      VectorFreq(CB, p2)*VectorScalar(CB, p2, k))
                             / (double)freqsum );
	else cent->vector[k] = 0;
    }
  cent->freq = freqsum;
}


/* ----------------------------------------------------------------- */


static void RemoveLastVector(CODEBOOK* CB, int a, int b, NEARESTTYPE NN[])
{
  int last = BookSize(CB)-1;
  int i;

  if( b != last )
    {
    CopyNode(&Node(CB,last), &Node(CB,b), VectorSize(CB));
    NN[b] = NN[last];
    for( i=0; i<BookSize(CB); i++ )
       {
       if( NN[i].Vector==last )  NN[i].Vector = b;
       }
    }
}


/*--------------------------------------------------------------------*/


static void RemoveEmptyPartitions(TRAININGSET*  TS,
                                  CODEBOOK*     CB,
                                  PARTITIONING* P)
{
  int i,last;

  for( i=0; i<BookSize(CB); i++ )
     {
     if( VectorFreq(CB,i)==0  &&  BookSize(CB)>Value(CodebookSize) )
        {
        last = BookSize(CB)-1;
        CopyNode(&Node(CB,last), &Node(CB,i), VectorSize(CB));
        JoinPartitions(TS, P, i, last);
        BookSize(CB)--;
        i--;
        }
     }
}


/*--------------------------------------------------------------------*/


static void MarkClustersForRecalculation(CODEBOOK* CB,
                                         NEARESTTYPE NN[],
                                         int a, int b)
{
  int i;

  for( i=0; i<BookSize(CB); i++ )
    {
	NN[i].Recalculate = ( NN[i].Vector==a || NN[i].Vector==b ) ? YES : NO;
	}
  NN[a].Recalculate = YES;
}


/*--------------------------------------------------------------------*/


static void RecalculateDistances(CODEBOOK* CB, NEARESTTYPE NN[])
{
  int  i;
  int  nrec=0;

  for( i=0; i<BookSize(CB); i++ )
    {
    if( NN[i].Recalculate )
       {
       FindNearestCluster(CB, i, &(NN[i].Vector), &(NN[i].VectorDist));
       NN[i].Recalculate = NO;
       nrec++;
       }
    }
  if( Value(QuietLevel) >= 6 )
     {
     printf("Recalculated: %i / %i.\n", nrec, BookSize(CB));
     }
}


/*-------------------------------------------------------------------*/


static void MergeVectors(TRAININGSET*  TS,
                         CODEBOOK*     CB,
                         PARTITIONING* P,
                         NEARESTTYPE   NN[],
                         int           a,
                         int           b)
{
  MarkClustersForRecalculation(CB, NN, a, b);
  CentroidOfTwoVectors(CB, &Node(CB, a), a, b);
  JoinPartitions(TS, P, a, b);
  RemoveLastVector(CB, a, b, NN);
  BookSize(CB)--;
}


/*-------------------------------------------------------------------*/


static void PNNFast(TRAININGSET*   TS,
                    CODEBOOK*      CB,
                    PARTITIONING*  P,
                    int            ResultCBSize)
{
  NEARESTTYPE*  NN;
  int  a,b;

  if(Value(PNNParameter)==YES)  RemoveEmptyPartitions(TS, CB, P);
  NN = (NEARESTTYPE*) allocate(BookSize(CB) * sizeof(NEARESTTYPE));
  CalculateAllDistances(CB, NN);
  while( BookSize(CB) > ResultCBSize )
     {
     FindMergedPair(CB, NN, &a, &b);
     MergeVectors(TS, CB, P, NN, a, b);
     RecalculateDistances(CB, NN);
     }
  deallocate(NN);
}


/*====================  Initial generation  =========================*/


static void GenerateInitialSolution(TRAININGSET* TS,
                                    SOLUTION*    S,
                                    SASchedule*  SAS,
                                    double*      Error)
{
  AddGenerationMethod(&S->CB, ProgName);
  AddGenerationMethod(&S->CB, VersionNumber);
  BookSize(&S->CB) = Value(CodebookSize);
  PartitionCount(&S->P) = Value(CodebookSize);
  switch( Value(InitialGeneration) )
     {
     case RANDOM_CB:
           GenerateRandomCodebook(TS, &S->CB);
           GenerateOptimalPartitioning(TS, &S->CB, &S->P);
           break;
     case RANDOM_P:
           GenerateRandomPartitioning(TS, &S->P);
           GenerateOptimalCodebook(TS, &S->CB, &S->P);
           break;
     case RANDOM_P_NOISE:
           GenerateRandomPartitioning(TS, &S->P);
           GenerateOptimalCodebook(TS, &S->CB, &S->P);
           AddNoiseToCodebook(&S->CB, SAS);
           break;
     case RGLA:
           GenerateRandomCodebook(TS, &S->CB);
           GenerateOptimalPartitioning(TS, &S->CB, &S->P);
           IterateGLAForSolution(TS, S, Value(GLAIterations), MSE);
           break;
     case RGLA_NOISE:
           GenerateRandomCodebook(TS, &S->CB);
           GenerateOptimalPartitioning(TS, &S->CB, &S->P);
           IterateGLAForSolution(TS, S, Value(GLAIterations), MSE);
           AddNoiseToCodebook(&S->CB, SAS);
           break;
     case R_PNN_GLA:
           BookSize(&S->CB) = 2*Value(CodebookSize);
           PartitionCount(&S->P) = 2*Value(CodebookSize);
           GenerateRandomCodebook(TS, &S->CB);
           GenerateOptimalPartitioning(TS, &S->CB, &S->P);
	   PNNFast(TS, &S->CB, &S->P, Value(CodebookSize));
           S->optimality = OPT_NONE;
           IterateGLAForSolution(TS, S, Value(GLAIterations), MSE);
           AddNoiseToCodebook(&S->CB, SAS);
           break;
     default:
           GenerateRandomCodebook(TS, &S->CB);
           GenerateOptimalPartitioning(TS, &S->CB, &S->P);
           break;
     }

/*(*Error) = AverageErrorForSolution(TS, &S->CB, &S->P, MSE); */
  (*Error) = AverageErrorCBFast(TS, &S->CB, &S->P, MSE);
}


/*-------------------------------------------------------------------*/


static void GenerateInitialGeneration(TRAININGSET* TS,
                                      SOLUTION*    Ss[],
                                      SASchedule*  SAS,
                                      double       Error[])
{
  int i;

  if( BookSize(TS) < Value(CodebookSize) )
      {
      printf("ERROR: Training set size (%i) < requested codebook size (%i).\n",
             BookSize(TS), Value(CodebookSize));
      exit( -1 );
      }

  for (i=0; i<Value(GenerationSize); i++)
      {
      GenerateInitialSolution(TS, Ss[i], SAS, &Error[i]);
      }
}


/*=================  CROSSOVER operations for CB  ====================*/


static int FindNearestNeighborVector(CODEBOOK* CB, int v)
{
  llong MinError=MAXLLONG;
  llong e;
  int	MinIndex=0;
  int	i;

  for (i=0; i<BookSize(CB); i++)
	 {
	 e = VectorDist(Vector(CB,i), Vector(CB,v), VectorSize(CB));
	 if( (e<MinError) && (i!=v) ) { MinError = e; MinIndex = i; }
	 }
  return( MinIndex );
}


/*-------------------------------------------------------------------*/


static int FindFirstNonEmpty(CODEBOOK* CB, PARTITIONING* P, int v)
{
  int	i;

  for (i=0; i<BookSize(CB); i++)
	 {
     if(UniqueVectors(P,i) && i!=v) return(i);
	 }
  printf("WARNING: No non-empty clusters found.\n");
  return( 0 );
}


/*--------------------------------------------------------------------*/


static void CrossRandom(TRAININGSET* TS,
                        SOLUTION*    S1,
                        SOLUTION*    S2,
                        SOLUTION*    S3,
                        SASchedule*  SAS)
/* Generates Snew from S1 and S2 using RANDOM CROSSOVER method. */
{
  int      vsize=VectorSize(&S1->CB);
  int	   i,j,n;
  YESNO    found;
  double   error;

  for (i=0; i<Value(CodebookSize); i++)
    {
	found=NO;
	while(!found)
	   {
	   found=YES;
	   n = irand(0,Value(CodebookSize)-1);
	   for(j=0; j<i; j++)
          if(EqualVectors(Vector(CBX(i),n),Vector(&S3->CB,j),vsize))
             { found=NO; break; }
	   }
    CopyNode(&Node(CBX(i),n), &Node(&S3->CB,i), vsize);
	}
  error = GenerateOptimalPartitioning(TS, &S3->CB,&S3->P);
}


/*--------------------------------------------------------------------*/


static void ReplaceDuplicates(CODEBOOK* CB, CODEBOOK* CBreserve)
{
  int i,j;

  for ( i=0; i<BookSize(CB); i++ )
    for( j=i+1; j<BookSize(CB); j++ )
       {
       if( EqualVectors( Vector(CB,i), Vector(CB,j), VectorSize(CB) ))
          {
          GenerateRandomVector(CBreserve, CB, i);
          }
       }
}


/*--------------------------------------------------------------------*/


static void CrossCentroidDistance(TRAININGSET* TS,
                                  SOLUTION*    S1,
                                  SOLUTION*    S2,
                                  SOLUTION*    S3,
                                  SASchedule*  SAS)
/* Generates Snew from S1 and S2 using CENTROID DISTANCE method. */
{
  int   vsize = VectorSize(&S1->CB);
  int   bsize = BookSize(&S1->CB);
  int   i;
  int   r1 = irand(1,2);
  int   r2 = 3-r1;

  SortByCentroidDistance(TS, &S1->CB);
  SortByCentroidDistance(TS, &S2->CB);
  for (i=0; i<bsize/2; i++)
      CopyVector(Vector(CBX(r1),i), Vector(&S3->CB,i), vsize);
  for (i=bsize/2; i<bsize; i++)
      CopyVector(Vector(CBX(r2),i), Vector(&S3->CB,i), vsize);
  ReplaceDuplicates(&S3->CB,CBX(r1));
}


/*--------------------------------------------------------------------*/


static void CrossAdoptionSwap(TRAININGSET* TS,
							  SOLUTION*    S1,
							  SOLUTION*    S2,
							  SOLUTION*    S3,
							  long		   watch)
/* Generates Snew from S1 and S2 using Adoption Swap */
{
  llong  error;
  int    i, c, guess;
  int    vsize = VectorSize(&S1->CB);
  int    bsize = BookSize(&S1->CB);

  CopyCodebook(&S1->CB, &S3->CB);
  CopyPartitioning(&S1->P, &S3->P);
  for(i=0; i<bsize; i++)
     {
     guess = Map(&S3->P,FirstVector(&S2->P,i));
     c = FindNearestVector( &Node(&S2->CB,i), &S3->CB, &error, guess, EUCLIDEANSQ);
     CopyVector( Vector(&S2->CB,i), Vector(&S3->CB,c), vsize );
     LocalRepartitioningGeneral(TS, &S3->CB, &S3->P, c, EUCLIDEANSQ);
     RepartitionDueToNewVectorGeneral(TS, &S3->CB, &S3->P, c, EUCLIDEANSQ);
     GenerateOptimalCodebook(TS, &S3->CB, &S3->P);
     S3->optimality = OPT_CB;
     }
}


/*--------------------------------------------------------------------*/


static void CrossAdaptivePairwise(TRAININGSET* TS,
                                  SOLUTION*    S1,
                                  SOLUTION*    S2,
                                  SOLUTION*    S3,
                                  SASchedule*  SAS,
                                  long         watch)
/* Generates Snew from S1 and S2 using Adaptive Pairwise Swap */
{
  double e;
  llong  error;
  int    i, c;
  int    vsize = VectorSize(&S1->CB);
  int    bsize = BookSize(&S1->CB);

  CopyCodebook(&S1->CB, &S3->CB);
  CopyPartitioning(&S1->P, &S3->P);
  for(i=0; i<bsize; i++)
     if(irand(0,1)==0)
        {
        c = FindNearestVector( &Node(&S3->CB,i), &S2->CB, &error, i, EUCLIDEANSQ);
        CopyVector( Vector(&S2->CB,c), Vector(&S3->CB,i), vsize );
        }
  e = GenerateOptimalPartitioning(TS, &S3->CB, &S3->P);
  FillEmptyPartitions(TS, &S3->CB, &S3->P);
}


/*--------------------------------------------------------------------*/


static void CrossPNNold(TRAININGSET* TS,
                        SOLUTION*    S1,
                        SOLUTION*    S2,
                        SOLUTION*    S3,
                        SASchedule*  SAS,
                        long         watch)
/* Generates Snew from S1 and S2 using PNN - old */
{
  CODEBOOK  TStmp;
  double e;
  int	 i;
  int	 booksize = BookSize(&S1->CB);
  int	 vsize = VectorSize(&S1->CB);

  if( Value(QuietLevel)>=5 )
    { printf("--> Combining codebooks by PNN."); PrintTime(watch); }

  /* Combine codebooks: CB3 <- CB1 + CB2 */
  BookSize(&S3->CB) = 2*Value(CodebookSize);
  for(i=0; i<booksize; i++)
     CopyNode( &Node(&S1->CB,i), &Node(&S3->CB,i), vsize );
  for(i=0; i<booksize; i++)
     CopyNode( &Node(&S2->CB,i), &Node(&S3->CB,booksize+i), vsize );

  /* Reduce codebook size from 2M -> M using PNN algorithm */
  if( Value(QuietLevel)>=5 )
     { printf("--> Reducing codebook by PNN."); PrintTime(watch); }
  CreateNewCodebook(&TStmp, BookSize(&S3->CB), &S3->CB);
  CopyCodebook(&S3->CB, &TStmp);
  SetVersionString(&TStmp, CBFILE);
  SetPNNParameters(0,0,0,0,0,0,0,0,0);
  PairwiseNearestNeighbour(&TStmp, &S3->CB, SAS, NULL, Value(CodebookSize));
  if( !BookSize(&S3->CB)==Value(CodebookSize) )
     printf("!!!!!!!!!!! %i !!!!!!!!! \n", BookSize(&S3->CB));
  e = GenerateOptimalPartitioning(TS, &S3->CB, &S3->P);
  FreeCodebook(&TStmp);

  if( Value(QuietLevel)>=5 )
     { printf("--> Crossing completed."); PrintTime(watch); }
}


/*--------------------------------------------------------------------*/


static void CrossPNNnew(TRAININGSET* TS,
                        SOLUTION*    S1,
                        SOLUTION*    S2,
                        SOLUTION*    S3,
                        SASchedule*  SAS,
                        long         watch)
/* Generates Snew from S1 and S2 using PNN - new */
{
  llong  d1, d2;
  int    i, p1, p2;
  int    booksize = BookSize(&S1->CB);
  int    vsize = VectorSize(&S1->CB);

  if( Value(QuietLevel)>=5 )
    { printf("--> Combining codebooks by PNN."); PrintTime(watch); }

  /* Combine codebooks: CB3 <- CB1 + CB2 */
  BookSize(&S3->CB) = 2*Value(CodebookSize);
  for(i=0; i<booksize; i++)
     CopyNode( &Node(&S1->CB,i), &Node(&S3->CB,i), vsize );
  for(i=0; i<booksize; i++)
     CopyNode( &Node(&S2->CB,i), &Node(&S3->CB,booksize+i), vsize );

  if( Value(QuietLevel)>=5 )
     { printf("--> Combining partitionings."); PrintTime(watch); }

  /* Combine partitions: P3 <- Best partitions from P1 and P2 */
  CopyPartitioning(&S1->P, &S3->P);
  PartitionCount(&S3->P) = 2*Value(CodebookSize);
  for(i=0; i<BookSize(TS); i++)
    {
    p1 = Map(&S1->P, i);
    p2 = Map(&S2->P, i);
    d1 = VectorDist(Vector(TS,i), Vector(&S1->CB,p1), vsize);
    d2 = VectorDist(Vector(TS,i), Vector(&S2->CB,p2), vsize);
    if(d2<d1) ChangePartition(TS, &S3->P, booksize+p2, i);
    }

  /* Calculate optimal codebook (2M) for combined partitioning */
  GenerateOptimalCodebook(TS, &S3->CB, &S3->P);
  S3->optimality = OPT_CB;

  /* Reduce codebook size from 2M -> M using PNN algorithm */
  if( Value(QuietLevel)>=5 )
     { printf("--> Reducing codebook by PNN."); PrintTime(watch); }
  PNNFast(TS, &S3->CB, &S3->P, Value(CodebookSize));

  if( Value(QuietLevel)>=5 )
     { printf("--> Crossing completed."); PrintTime(watch); }
}


/*--------------------------------------------------------------------*/


static void CrossIS(TRAININGSET* TS,
                        SOLUTION*    S1,
                        SOLUTION*    S2,
                        SOLUTION*    S3,
                        SASchedule*  SAS,
                        long         watch)
/* Generates Snew from S1 and S2 using IS  */
{
  llong  d1, d2;
  int    i, p1, p2;
  int    booksize = BookSize(&S1->CB);
  int    vsize = VectorSize(&S1->CB);
  PARTITIONING  P2;
  CreateNewPartitioning(&P2, TS, 2*Value(CodebookSize));

  if( Value(QuietLevel)>=5 )
    { printf("--> Combining codebooks by IS."); PrintTime(watch); }

  /* Combine codebooks: CB3 <- CB1 + CB2 */
  BookSize(&S3->CB) = 2*Value(CodebookSize);
  for(i=0; i<booksize; i++)
     CopyNode( &Node(&S1->CB,i), &Node(&S3->CB,i), vsize );
  for(i=0; i<booksize; i++)
     CopyNode( &Node(&S2->CB,i), &Node(&S3->CB,booksize+i), vsize );

  if( Value(QuietLevel)>=5 )
     { printf("--> Combining partitionings."); PrintTime(watch); }

  /* Combine partitions: P3 <- Best partitions from P1 and P2 */
  CopyPartitioning(&S1->P, &S3->P);
  PartitionCount(&S3->P) = 2*Value(CodebookSize);
  for(i=0; i<BookSize(TS); i++)
    {
    p1 = Map(&S1->P, i);
    p2 = Map(&S2->P, i);
    d1 = VectorDist(Vector(TS,i), Vector(&S1->CB,p1), vsize);
    d2 = VectorDist(Vector(TS,i), Vector(&S2->CB,p2), vsize);
    if(d2<d1) ChangePartition(TS, &S3->P, booksize+p2, i);
    }

  /* Calculate optimal codebook (2M) for combined partitioning */
  GenerateOptimalCodebook(TS, &S3->CB, &S3->P);
  S3->optimality = OPT_CB;

  /* Reduce codebook size from 2M -> M using IS algorithm */
  if( Value(QuietLevel)>=5 )
     { printf("--> Reducing codebook by IS."); PrintTime(watch); }
  IterativeShrinking(TS, &S3->CB, &S3->P, &P2, Value(CodebookSize),
      Value(RemovalCostCalculation), Value(Update), Value(QuietLevel));

  if( Value(QuietLevel)>=5 )
     { printf("--> Crossing completed."); PrintTime(watch); }
  FreePartitioning(&P2);
}


/*--------------------------------------------------------------------*/


static void CrossSolutions(TRAININGSET* TS,
                           SOLUTION*    S1,
                           SOLUTION*    S2,
                           SOLUTION*    S3,
                           SASchedule*  SAS,
                           long         watch)
{
  switch( Value(CrossingMethod) )
    {
    case 0:  CopySolution(S1, S3); break;
    case 1:  CrossRandom(TS, S1, S2, S3, SAS); break;
    case 2:  CrossCentroidDistance(TS, S1, S2, S3, SAS); break;
    case 3:  /* LP */
    case 4:  CrossAdoptionSwap(TS, S1, S2, S3, watch); break;
    case 5:  CrossAdaptivePairwise(TS, S1, S2, S3, SAS, watch);  break;
    case 6:  CrossPNNold(TS, S1, S2, S3, SAS, watch);  break;
    case 7:  CrossPNNnew(TS, S1, S2, S3, SAS, watch);  break;
    case 8:  /* Adaptive IS, crossing is the same */
    case 9:  CrossIS(TS, S1, S2, S3, SAS, watch);  break;
    default: break;
    }
}


/*==================  Genetic operations for CB  =====================*/


static void CalculateSurvivors(void)
{
  int s;

  switch( Value(SelectionMethod) )
     {
     case ROULETTE:
             Survivors = 1;
             CrossSetSize = Value(GenerationSize);
             break;
     case ELITIST1:
     case ZIGZAG:
             s=1; while ( (s*(s-1)/2) < Value(GenerationSize) ) s++;
             Survivors = 1;
             CrossSetSize = s;
             break;
     case ELITIST2:
             s=1; while ( (s*(s+1)/2) < Value(GenerationSize) ) s++;
             CrossSetSize = s;
             Survivors = s;
             break;
     }
}


/*--------------------------------------------------------------------*/


static int SelectRandom(double Error[])
		/* Selects a random codebook by roulette wheel. */
{
  double total=0;
  double s;
  int	 i;

  for (i=0; i<CrossSetSize; i++)
	  total += 1/(1+Error[i]);

  s = total * frand();
  total = 0;

  for (i=0; i<CrossSetSize; i++)
	 {
     total += 1/(1+Error[i]);
	 if (total >= s) return(i);
	 }

  return(0);
}


/*--------------------------------------------------------------------*/


static void SelectNextPair(double Error[], int* s1, int* s2)
{
  /* ELITIST: (0,1) (0,2) ... (0,c-1)
			  (1,2) (1,3) ... (1,c-1)
			  ...
			  (c-2,c-1)
  */

  switch( Value(SelectionMethod) )
     {
     case ROULETTE:
	{
        (*s1) = SelectRandom(Error);
        (*s2) = SelectRandom(Error);
        while((*s1)==(*s2))  (*s2) = SelectRandom(Error);
        break;
	}
     case ELITIST1:
     case ELITIST2:
	{
	if( (++(*s2)) == CrossSetSize )  { (*s1)++; (*s2)=(*s1)+1; }
	if( (*s1) == CrossSetSize ) { printf("Crossover iterator overflow.\n"); exit(-1); }
	break;
        }
   case ZIGZAG:
        {
        do {
           if (((*s1)+(*s2))%2) { (*s1)++; (*s2)=max((*s2)-1,0); }
           else                 { (*s2)++; (*s1)=max((*s1)-1,0); }
           } while ((*s1) < (*s2));
		if( (*s1) == CrossSetSize ) { printf("Crossover iterator overflow.\n"); exit(-1); }
        break;
        }
   default:
		{
		(*s1)=0; (*s2)=0; break;
		}
   }
}


/*--------------------------------------------------------------------*/


static void Mutate(TRAININGSET* TS, SOLUTION* S, long watch)
{
  int i;

  if( Value(Mutations)==0 ) return;

  if( Value(QuietLevel)>=5 )
     { printf("Generating %i mutations.",Value(Mutations));
	   PrintTime(watch); }

  for (i=0; i<=Value(Mutations); i++)
	 {
     GenerateRandomVector( TS, &S->CB, irand(0,BookSize(&S->CB)-1) );
	 }
  S->optimality = OPT_NONE;
}


/*--------------------------------------------------------------------*/


static void IterateGLASelective(TRAININGSET* TS,
                                SOLUTION*    Ss[],
                                double       Error[])
{
  int i;
  int s=min(Value(SelectiveGLA),Value(GenerationSize));

  for(i=0; i<s; i++)
     {
     IterateGLAForSolution(TS, Ss[i], Value(GLAIterations), MSE);
     /*Error[i] = AverageErrorForSolution(TS, &Ss[i]->CB, &Ss[i]->P, MSE);*/
     Error[i] = AverageErrorCBFast(TS, &Ss[i]->CB, &Ss[i]->P, MSE);
     }
}


/*--------------------------------------------------------------------*/


static void GenerateNextGeneration(TRAININGSET* TS,
                                   SOLUTION*    Sold[],
                                   SOLUTION*    Snew[],
				   SASchedule*	SAS,
                                   double       Error[],
                                   int          generation,
                                   long         watch)
{
   int	i, gs;
   int	first=0;
   int	second=0;

  if( Value(QuietLevel)>=5 )
      printf("\nGenerating %i. generation.\n",generation);

  /* Generate GenSize-Survivors new solutions. Survivors already exist. */
  for(i=Survivors; i<Value(GenerationSize); i++)
     {
     if( Value(QuietLevel)>=5 ) printf("--> Generating %i. solution.\n",i);
     SelectNextPair(Error, &first, &second);
     CrossSolutions(TS, Sold[first], Sold[second], Snew[i], SAS, watch);
     Mutate(TS, Snew[i], watch);
     if( Value(GLAIterations) && !Value(SelectiveGLA) )
        {
        IterateGLAForSolution(TS, Snew[i], Value(GLAIterations), MSE);
        if( Value(QuietLevel)>=5 )
          { printf("--> Normal GLA completed."); PrintTime(watch); }
        }
     AddNoiseToCodebook(&Snew[i]->CB, SAS);
     /*Error[i] = AverageErrorForSolution(TS, &Snew[i]->CB, &Snew[i]->P, MSE);*/
     Error[i] = AverageErrorCBFast(TS, &Snew[i]->CB, &Snew[i]->P, MSE);
     if( Value(QuietLevel)>=4 )
        {
        printf("Crossing pair %i and %i.\n",first,second);
        if((generation==NumberOfIterations) && Value(QuietLevel)>=6 )
           {
		   printf("\nCrossover codebooks:\n");
		   PrintCodebook(&Sold[first]->CB);
		   PrintCodebook(&Sold[second]->CB);
		   printf("\nNEW codebook:\n");
		   PrintCodebook(&Snew[i]->CB);
		   printf("\n\n");
           }
        }
     }

  /* Adaptive IS: Increase population by adding new solution. */
  if( Value(CrossingMethod)==8 )
     {
     gs = Value(GenerationSize);
     SetValue(GenerationSize, min( 2*gs, MaxGenerations) );
     CalculateSurvivors();
     for(i=gs; i<Value(GenerationSize); i++)
        {
        GenerateInitialSolution(TS, Snew[i], SAS, &Error[i]);
        }
     }
}


/*===========================  M A I N  ==============================*/


int main(int argc, char* argv[])
{
  TRAININGSET   TS;
  SOLUTION*     Sold[MaxGenerations];
  SOLUTION*     Snew[MaxGenerations];
  SASchedule    SAS;
  char          TSName[64] = { '\0' };
  char          CBName[64] = { '\0' };
  double	Error[MaxGenerations];
  double        PrevError=0;
  long          watch;
  int			g;
  int           GoOn=YES;

  ParameterInfo paraminfo[3] = { { TSName, FormatNameTS, 0, INFILE  },
                                 { CBName, FormatNameCB, 0, OUTFILE } };

  ParseParameters(argc, argv, 2, paraminfo);
  initrandom(Value(RandomSeed));
  PrintOperatingInfo(TSName, CBName);
  ReadTrainingSet(TSName, &TS);
  Initialize(&TS, Sold, Snew, &SAS);
  CalculateSurvivors();
  SetWatch(&watch);

  GenerateInitialGeneration(&TS, Snew, &SAS, Error);
  SortSolutions(&TS, Snew, Error);
  PrintProgress(&Snew[0]->CB, &SAS, 0, PrevError, Error, watch);

  for(g=1; g<=NumberOfIterations; g++)
     {
     PrevError=Error[0];
     CopySolutions(Snew, Sold);
     GenerateNextGeneration(&TS, Sold, Snew, &SAS, Error, g, watch);
     SortSolutions(&TS, Snew, Error);
     if( Value(SelectiveGLA) )
        {
        IterateGLASelective(&TS, Snew, Error);
        SortSolutions(&TS, Snew, Error);
        if( Value(QuietLevel)>=5 )
          { printf("--> Selective GLA completed."); PrintTime(watch); }
        }
     /* if( Value() && NoImprovement ) IncreaseGenerationSize(); */
     PrintProgress(&Snew[0]->CB, &SAS, g, PrevError, Error, watch);
     DecreaseTemperature(&SAS);
	 #if defined(__MSDOS__)
	 if (kbhit()) GoOn=(getch()!=27);
	 #endif
     if(!GoOn || NoImprovement) break;
     }

  PrintResult(watch, PrintableError(Error[0], &Snew[0]->CB), g-1 );
  OutputCodebook(CBName, &Snew[0]->CB, &TS);
  ShutDown(&TS, Sold, Snew);
  checkmemory();
  return( 0 );
}

