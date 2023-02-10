/*-------------------------------------------------------------------*/
/* CBIS.C         Pasi Fränti & Olli Virmajoki & Timo Kaukoranta     */
/*                                                                   */
/* - Iterative shrinking                                             */
/*                                                                   */
/*-------------------------------------------------------------------*/

#define ProgName       "CBIS"
#define VersionNumber  "Version 0.25"
#define LastUpdated    "31.8.03 OV"

/* ----------------------------------------------------------------- */

#include <assert.h>
#include <float.h>
#include <time.h>
#include <values.h>

#define  FACTFILE  "cbis.fac"

#include "parametr.c"
#include "cb.h"
#include "file.h"
#include "interfc.h"
#include "memctrl.h"
/*#include "solution.h"*/ 
#include "sortcb.h"


/*---------------------------  M i s c	------------------------------*/

/*#define  STATISTICS 1*/
#ifdef STATISTICS
#define STATCODE(a)a
#else 
#define STATCODE(a)
#endif

STATCODE
(
llong ndistcalcs=0LL;
llong ndistcalcstotal=0LL;
llong nidistcalcstotal=0LL;
llong n2distcalcs=0LL;
llong n2distcalcstotal=0LL;
llong nSize=0LL;
llong nSizetotal=0LL;
llong niteration=0LL;
)

#define  swap(a,b)      { int t = a; a = b; b = t; }
#define  round(a)		((llong) (a+0.5))

#if defined(max) || defined(min)
#undef max
#undef min
#endif
#define  max(a,b) ((a) > (b) ? (a) : (b))
#define  min(a,b) ((a) < (b) ? (a) : (b))
#define  ScaleBetween(a,b,c) max((a), min((b), (c)))


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

  if(Value(QuietLevel)==0) return;
  s=t%60; t=t/60; m=t%60; h=t/60;
  printf("Time: %li:%02li:%02li\n", h, m, s);
}


/*============================  P R I N T  ===============================*/


void PrintInfo(void)
{
  printf("\n%s   %s   %s\n", ProgName, VersionNumber, LastUpdated);
  printf("Generates a codebook using Iterative Shrinking.\n");
  printf("Usage: %s [%coption] <training set> <codebook>\n", ProgName, OPTION_SYMBOL);
  printf("For example: %s bridge tmp\n", ProgName);
  printf("\n  Options:\n");
  PrintOptions();
  printf("\n%s   %s   %s\n", ProgName, VersionNumber, LastUpdated);
}


/*-------------------------------------------------------------------*/


static void PrintOperatingInfo(char* TSName, char* InitCBName,
                               char* CBName, char* PAName)
{
  if( Value(QuietLevel) >= 2 )
    {
    printf("\n%s   %s   %s\n", ProgName, VersionNumber, LastUpdated);
    printf("Training Set:     %s\n", TSName);
    if(InitCBName[0] != '\0')
    printf("Initial codebook: %s\n", InitCBName);
	 printf("Result codebook:  %s\n", CBName);
	 if(Value(SavePartition))
	 printf("Result partition: %s\n", PAName);
    printf("\n");
    PrintSelectedOptions();
    printf("\n");
    }
}


/*-------------------------------------------------------------------*/


static void PrintProgress(TRAININGSET*  TS,
                          CODEBOOK*     CB,
                          PARTITIONING* P,
                          int           Quiet,
                          long          watch)
{
    double Error=0.0;
    PARTITIONING PT;

    if(Quiet<2) return;
    STATCODE
        (
        ndistcalcstotal += ndistcalcs;
        n2distcalcstotal += n2distcalcs;
        nSizetotal += nSize;
        niteration += 1;
        )
    

    switch(Quiet)
	{ 
        case 2:  /*if(BookSize(CB) <=256)
                 {
                     CreateNewPartitioning(&PT, TS, BookSize(CB));
                     CopyPartitioning(P, &PT);
                     GenerateOptimalPartitioning(TS, CB, &PT);
                     Error = AverageErrorForSolution(TS, CB, &PT, MSE);
                     Error = PrintableError(Error, CB);
                     printf("Size=%3i: %9.4f \n",BookSize(CB), Error);
                     FreePartitioning(&PT);
                 }*/
                 /*    printf("Size=%3i: %9.4f ndistcalcs: %7.0f n2distcalcs: %7.0f \n", 
                     BookSize(CB), Error, (double)ndistcalcs, (double)n2distcalcs);
                 */
                 /*break;*/
        default:
             printf("Size=%3i: %9.4f  ",  BookSize(CB), Error);
             PrintTime(watch); break;
	}
    fflush(stdout);
    STATCODE
        (
        ndistcalcs=0LL;
        n2distcalcs=0LL;
        nSize=0LL;
        )
}


/*======================  SOLUTION HANDLING  ========================*/

 
static void InitializeSolution(char*		 InitCBName,
							   TRAININGSET*  TS,
							   CODEBOOK*	 CB,
                               PARTITIONING* P1,
                               PARTITIONING* P2)
{
  if( Value(CodebookSize)>=BookSize(TS) )
  {
  		ErrorMessage("Codebook size (%i) > training set size (%i).\n", Value(CodebookSize),BookSize(TS));
	  	ExitProcessing(-1);
  }

  if( InitCBName[0] != '\0' ) /* Initial CB is given. */
  {
     ReadCodebook(InitCBName, CB);
     AddGenerationMethod(CB, "Iterative shrinking (from CB0)");
  }
  else /* Create initial CB from TS */
  {
     CreateNewCodebook(CB, BookSize(TS), TS);
     CopyCodebook(TS, CB);
     SetVersionString(CB, CBFILE);
     AddGenerationMethod(CB, "Iterative shrinking");
  }

  CreateNewPartitioning(P1, TS, BookSize(CB));
  CreateNewPartitioning(P2, TS, BookSize(CB));
}


/*-------------------------------------------------------------------*/


static void SaveSolution(char*         CBName,
                         char*         PAName,
                         TRAININGSET*  TS,
                         CODEBOOK*     CB,
                         PARTITIONING* P)
{
  int ndups;

  ndups = DuplicatesInCodebook(CB);
  if( ndups && Value(QuietLevel)>=2 )
  {
     printf("WARNING: Final CB contains %i duplicates.\n", ndups);
  }
  SortCodebook(CB, DATA_DESCENDING);
  WriteCodebook(CBName, CB, Value(OverWrite));
  if(Value(SavePartition))
  {
     WritePartitioning(PAName, P, TS, Value(OverWrite));
  }
}


/*-------------------------------------------------------------------*/


static void ShutDown(TRAININGSET* TS, CODEBOOK* CB,
                     PARTITIONING* P1, PARTITIONING* P2)
{
  FreeCodebook(TS);
  FreeCodebook(CB);
  FreePartitioning(P1);
  FreePartitioning(P2);
}


/*===================== Iterative shrinking =========================*/


static void PrintPartitions(TRAININGSET*  TS,
                            CODEBOOK*     CB,
                            PARTITIONING* P1,
                            PARTITIONING* P2,
                            llong         RemovalCost[])
{
  int c,i;

  for( c=0; c<BookSize(CB); c++ )
	 {
     printf("%3i: P1=",c);
	 for( i=FirstVector(P1,c); !EndOfPartition(i); i=NextVector(P1,i) )
        printf("%2i(%3i,%3i)", i, VectorScalar(TS,i,0), VectorScalar(TS,i,1));
     printf("  P2=");
	 for( i=FirstVector(P1,c); !EndOfPartition(i); i=NextVector(P1,i) )
        printf("<%2i>",Map(P2,i));
     printf("   cost=%4li\n", (long) RemovalCost[c]);
	 }
  printf("\n");
}


/*-------------------------------------------------------------------*/


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
     STATCODE
         (
         n2distcalcs++;
         )
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
										  PARTITIONING* P2,
                                          long          watch)
{
  int  i,c;

  for( i=0; i<BookSize(TS); i++ )
	 {
     c = FindSecondNearest( TS, CB, i, -1, Map(P1,i) );
	 ChangePartition(TS, P2, c, i);
     }
  if( Value(QuietLevel)>=4 )
     {
     printf("\nSecondary partitions completed. ");
     PrintTime(watch);
     }
}


/*-------------------------------------------------------------------*/


static void GenerateInitialPartitions(TRAININGSET*  TS,
                                      CODEBOOK*     CB,
                                      PARTITIONING* P1,
                                      PARTITIONING* P2,
                                      long          watch)
{
  if( BookSize(TS) == BookSize(CB) )
	 {
	 PutAllInOwnPartition(TS, P1);
	 GenerateSecondaryPartitioning(TS, CB, P1, P2, watch);
	 }
  else /* We have an initial codebook. */
     {
     GenerateOptimalPartitioning(TS, CB, P1);
     GenerateSecondaryPartitioning(TS, CB, P1, P2, watch);
	}
}


/*--------------------------------------------------------------------*/


static void CalculateRemovalCosts(TRAININGSET*  TS,
                                  CODEBOOK*     CB,
                                  PARTITIONING* P1,
                                  PARTITIONING* P2,
                                  llong         RemovalCost[])
/*   Simple variant; (not in use)   */
{
  int	 i,p1,p2;
  llong d1, d2;
  double factor;

  for( i=0; i<BookSize(CB); i++ )
	 {
	 RemovalCost[i] = 0;
	 }
  for( i=0; i<BookSize(TS); i++ )
	 {
	 p1 = Map(P1,i);
     p2 = Map(P2,i);

     d1 = VectorFreq(TS,i) * VectorDist(Vector(TS,i), Vector(CB,p1), VectorSize(TS));
	 d2 = VectorDist(Vector(TS,i), Vector(CB,p2), VectorSize(TS));
     factor=((double)VectorFreq(TS,i) * (double)VectorFreq(CB,p2)) /
         ((double)VectorFreq(TS,i) + (double)VectorFreq(CB,p2));
     d2 = round( factor*(double)d2);
     RemovalCost[p1] += (d2-d1);
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
  if( Value(QuietLevel)>=4 ) printf("Cluster %i selected.\n",*cluster);
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
        bn->vector[k] = round((double)Counter[ti * VectorSize(Cents) + k]/(double)Freq[ti]);
    }
}


/*-------------------------------------------------------------------*/


static void UpdatePrimaryPartition(TRAININGSET* TS, 
                                   PARTITIONING* P1, PARTITIONING* P2, 
                                   int* Size, int Elem[], int Freq[], llong Counter[], 
                                   int obsolete)
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
    STATCODE
            (
                nSize=*Size;
            )
}


/*-------------------------------------------------------------------*/


static void CalcPartitionCentroids(CODEBOOK* CB, PARTITIONING* P1, int Size, int Elem[])
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
                                     int Quiet,
                                     long watch)
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


static void RemoveObsolete(TRAININGSET* TS, CODEBOOK* CB, PARTITIONING* P1, 
                           PARTITIONING* P2, int obsolete)
{
    int last = BookSize(CB)-1;
    CopyNode(&Node(CB, last), &Node(CB, obsolete), VectorSize(CB));
    BookSize(CB)--; 
    JoinPartitions(TS, P1, obsolete, last);
    JoinPartitions(TS, P2, obsolete, last); 
}


/*-------------------------------------------------------------------*/


static llong CalcSubRemCost(CODEBOOK* CB, int Freq[], CODEBOOK* Cents, int ti)
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
                             int* Size, int Elem[], int Freq[], llong Counter[], 
                             CODEBOOK* Cents, llong RemCost[], int RemovalCalc)
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
        STATCODE
            (
            ndistcalcs+=1;
            )
        if(RemovalCalc==SIMPLE)
        {
            RemCost[j] += CalcMergeCost(TS, CB, i, Map(P2,i));
            STATCODE
                (
                    ndistcalcs+=1;
                )
        }
     }
     if(RemovalCalc==EXACT)
     {
     for( i = 0; i <*Size; i++)
        CalcCentroids(Elem[i], Freq, Counter, Cents, &Node(Cents,Elem[i]));
     for ( i = FirstVector(P1, j); ! EndOfPartition(i); i = NextVector(P1, i) )
     {  
        RemCost[j] += CalcRemCost(TS, Cents, i, Map(P2,i));
        STATCODE
            (
            ndistcalcs+=1;
            )  
     }
     for( l = 0; l <*Size; l++)
     {
        RemCost[j] += CalcSubRemCost(CB, Freq, Cents, Elem[l]);
        STATCODE
            (
            ndistcalcs+=1;
            )
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


static void PrintResult(int Quiet)
{
   if(Quiet<2) return;
   STATCODE
       (
        printf("\nSizeAverage: %7.2f distcalcsAverage: %7.2f distcalcs2Average: %7.2f",
            (double)nSizetotal/niteration, (double)ndistcalcstotal/niteration, (double)n2distcalcstotal/niteration);
        printf("\nnidistcalcs: %7.0f ndistcalcs: %7.0f n2distcalcs: %7.0f ",
            (double)nidistcalcstotal, (double)ndistcalcstotal, (double)n2distcalcstotal);   
       )     
}


/*-------------------------------------------------------------------*/


static void FreeDataStructures(llong*  RemCost, int* Elem,
                               int* Freq, llong* Counter, CODEBOOK* Cents)
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
  long      watch;
  double    Error=0;
  int       obsolete;
  
  SetWatch(&watch);
  CreateDataStructures(CB, &RemCost, &Elem, &Freq, &Counter, &Cents);
  GenerateInitialPartitions(TS, CB, P1, P2, watch);
  PrintProgress(TS, CB, P1, Quiet, watch);
  STATCODE
      (
        niteration=0LL;
        nidistcalcstotal=n2distcalcstotal;
        n2distcalcstotal=0LL;
      )
  while( BookSize(CB)>ResultCBSize)
  {
    CalcRemovalCosts(TS, CB, P1, P2, &Size, Elem, Freq, Counter, &Cents, RemCost, RemovalCalc);
    FindClusterToBeRemoved(CB, RemCost, &obsolete);
    UpdatePrimaryPartition(TS, P1, P2, &Size, Elem, Freq, Counter, obsolete); 
    CalcPartitionCentroids(CB, P1, Size, Elem);
    UpdateSecondaryPartition(TS, CB, P1, P2, Size, Freq, obsolete, UpdateStyle, Quiet, watch);
    RemoveObsolete(TS, CB, P1, P2, obsolete);
    PrintProgress(TS, CB, P1, Quiet, watch);
  }

  GenerateOptimalPartitioning(TS, CB, P1);
  Error = AverageErrorForSolution(TS, CB, P1, MSE);
  if( Quiet >= 0 ) printf("%9.4f  ", Error);
  if( Quiet >= 1 ) PrintTime(watch);

  PrintResult(Quiet);
  FreeDataStructures(RemCost, Elem, Freq, Counter, &Cents);
}


/*===========================  M A I N  ==============================*/


int main(int argc, char* argv[])
{
  TRAININGSET   TS;
  CODEBOOK		CB;
  PARTITIONING  P1;
  PARTITIONING  P2;
  char			TSName[64]	   = { '\0' };
  char			CBName[64]	   = { '\0' };
  char			InitCBName[64] = { '\0' };
  char			PAName[64]	   = { '\0' };
  
  ParameterInfo paraminfo[3] = { { TSName,	   FormatNameTS, 0, INFILE	},
								 { InitCBName, FormatNameCB, 1, INFILE	},
								 { CBName,	   FormatNameCB, 0, OUTFILE } };

  ParseParameters(argc, argv, 3, paraminfo);
  PickFileName(CBName, PAName);
  CheckFileName(PAName, FormatNamePA);
  PrintOperatingInfo(TSName, InitCBName, CBName, PAName);
  ReadTrainingSet(TSName, &TS);
  InitializeSolution(InitCBName, &TS, &CB, &P1, &P2);
  switch( Value(ISMethod) )
  {
  case NOTHING: 
        printf("\nWe just do nothing!\n\n");
        break;
  case IS: /* Iterative Shrinking */
        IterativeShrinking(&TS, &CB, &P1, &P2, Value(CodebookSize), 
            Value(RemovalCostCalculation), Value(Update), Value(QuietLevel));
        break;
  }
  SaveSolution(CBName, PAName, &TS, &CB, &P1);
  ShutDown(&TS, &CB, &P1, &P2);
  return( 0 );
}

