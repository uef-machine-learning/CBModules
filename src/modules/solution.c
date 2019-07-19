/*-------------------------------------------------------------------*/
/* SOLUTION.C      Juha Kivijärvi                                    */
/*                                                                   */
/* Data structures and interface for solution abstraction            */
/*                                                                   */
/*-------------------------------------------------------------------*/

#define ProgName           "SOLUTION"
#define VersionNumber      "Version 0.11"
#define LastUpdated        "6.4.2001"

/* ----------------------------------------------------------------- */
/*      Changelog:                                                   */
/*                                                                   */
/*      0.11  added SortSolution                                     */
/*      0.10  added SolutionToBitString, BitStringToSolution and     */
/*            SizeOfSolution (JK)                                    */
/*                                                                   */
/* ----------------------------------------------------------------- */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <values.h>

#include "memctrl.h"
#include "file.h"
#include "interfc.h"
#include "cb.h"
#include "solution.h"
#include "owntypes.h"

#define  MAX_GLA   100000000

/* ===================   Solution interface   ====================== */


CBFILETYPE ReadSolution(TRAININGSET* TS, SOLUTION* S, char* FileName)
{
  CBFILETYPE ftype;
  
  ftype = DetermineCBFileTypeConsideringOrder(FileName, 23);
  
  switch (ftype)
    {
    case NOTFOUND:
      {
      ErrorMessage("ERROR: couldn't find the file: %s\n", FileName);
      ExitProcessing(-1);
      }
    case CBFILE:
      {
      ReadCodebook(FileName, &S->CB);
      CreateNewPartitioning(&S->P, TS, BookSize(&S->CB));
      GenerateOptimalPartitioning(TS, &S->CB, &S->P);
      S->optimality = OPT_PA;
      break;
      }
    case PAFILE:
      {
      ReadPartitioning(FileName, &S->P, TS);
      CreateNewCodebook(&S->CB, S->P.PartitionCount, TS);
      GenerateOptimalCodebook(TS, &S->CB, &S->P);
      S->optimality = OPT_CB;
      break;
      }
    default:
      {
      ErrorMessage("ERROR: invalid file type\n");
      ExitProcessing(-1);
      }    
    }

  return ftype;
}
    

/*---------------------------------------------------------------------*/


void CreateNewSolution(TRAININGSET* TS, SOLUTION* S, int booksize)
{
  CreateNewCodebook(&S->CB, booksize, TS);
  CreateNewPartitioning(&S->P, TS, booksize);
  S->optimality = OPT_NONE;
}
  

/*---------------------------------------------------------------------*/


void FreeSolution(SOLUTION* S)
{
  FreeCodebook(&S->CB);
  FreePartitioning(&S->P);
}  


/*---------------------------------------------------------------------*/


void IncreaseSolutionSize(SOLUTION* S, int newsize)
{
  IncreaseCodebookSize(&S->CB, newsize);
  IncreaseNumberOfPartitions(&S->P, newsize);
}

  
/*---------------------------------------------------------------------*/


void CopySolution(SOLUTION* sourceS, SOLUTION* destS)
{
  CopyCodebook(&sourceS->CB, &destS->CB);
  CopyPartitioning(&sourceS->P, &destS->P);
  destS->optimality = sourceS->optimality;
}
  

/*---------------------------------------------------------------------*/


void ChangeCodeVectorInSolution(TRAININGSET*   TS, 
                                SOLUTION*      S, 
                                int            item, 
                                VECTORTYPE     v,
                                ERRORFTYPE     errorf,
                                YESNO          dolocalrep)
{
/*  Puts vector v to place item. If dolocalrep, does
    local repartitioning. Changes in optimality:
       before        after
       OPT_NONE      OPT_NONE
       OPT_CB        OPT_NONE
       OPT_PA        (dolocalrep?) OPT_PA : OPT_NONE
       OPT_BOTH      (dolocalrep?) OPT_PA : OPT_NONE      */

  CopyVector(v, Vector(&S->CB, item), VectorSize(&S->CB));

  if (dolocalrep)
    {
    LocalRepartitioningGeneral(TS, &S->CB, &S->P, item, DistType(errorf));
    RepartitionDueToNewVectorGeneral(TS, &S->CB, &S->P, item,
                                                        DistType(errorf));
    if (S->optimality == OPT_BOTH ) S->optimality = OPT_PA;
    if (S->optimality == OPT_CB ) S->optimality = OPT_NONE;
    }
  else S->optimality = OPT_NONE;
}
      
    
/*---------------------------------------------------------------------*/


void ChangePartitionInSolution (TRAININGSET* TS,
                                SOLUTION*    S,
                                int          item,
                                int          Pindex,
                                ERRORFTYPE   errorf,
                                YESNO        updatecentroids)
{
/* Moves item to partition Pindex. If updatecentroids,
   updates partition centroids. Changes in optimality:
      before             after
      OPT_NONE           OPT_NONE
      OPT_CB             (updatecentroids?) OPT_CB : OPT_NONE
      OPT_PA             OPT_NONE
      OPT_BOTH           (updatecentroids?) OPT_CB : OPT_NONE   */

  int oldP = Map(&S->P, item);

  ChangePartition(TS, &S->P, Pindex, item);
  
  if (updatecentroids)
    {
    if (CCFreq(&S->P, oldP) > 0)
      PartitionCentroid(&S->P, oldP, &Node(&S->CB, oldP));
    if (CCFreq(&S->P, Pindex) > 0)
      PartitionCentroid(&S->P, Pindex, &Node(&S->CB, Pindex));
    if ( S->optimality == OPT_PA ) S->optimality = OPT_NONE;
    if ( S->optimality == OPT_BOTH ) S->optimality = OPT_CB;
    }
  else S->optimality = OPT_NONE;
}
                                 
                                
/*---------------------------------------------------------------------*/


void SolutionToBitString(SOLUTION* S, void* bs)
/* Saves all the data in a solution into a bit string (except
   for vector names).
   Space for bs should be allocated. May need changes when the
   definition of SOLUTION is changed. */
{
int i;
int ts;                    /* size(s) of a memory area */
void* bsp;                 /* pointer to free space */

bsp = bs;
/* Codebook */
ts = MaxVersionLength + 8 * sizeof(int) + MaxGenMethodLength;
memcpy(bsp, (void*) &(S->CB), ts);
bsp += ts;
ts = VectorSize(&(S->CB)) * sizeof(VECTORELEMENT);
for (i=0; i<S->CB.AllocatedSize; i++)
    {
    memcpy(bsp, (void*) Node(&(S->CB),i).vector, ts);
    bsp += ts;
    memcpy(bsp, (void*) &(Node(&(S->CB),i).vmean), 2 * sizeof(int));
    bsp += 2 * sizeof(int);
    }
memcpy(bsp, (void*) &(S->CB.AllocatedSize), sizeof(int));
bsp += sizeof(int);

/* Partitioning */
ts = MaxVersionLength + 2 * sizeof(int);
memcpy(bsp, (void*) &(S->P), ts);
bsp += ts;
ts = S->P.TSsize * sizeof(int);
memcpy(bsp, (void*) S->P.Map, ts);
bsp += ts;
memcpy(bsp, (void*) S->P.First, ts);
bsp += ts;
memcpy(bsp, (void*) S->P.Next, ts);
bsp += ts;
ts = S->P.Vsize * sizeof(llong);
for (i=0; i<S->P.PartitionCount; i++)
    {
    memcpy(bsp, (void*) S->P.CC[i].counter, ts);
    bsp += ts;
    memcpy(bsp, (void*) &(S->P.CC[i].freq), sizeof(int));
    bsp += sizeof(int);
    }
ts = S->P.PartitionCount * sizeof(int);
memcpy(bsp, (void*) S->P.Uniques, ts);
bsp += ts;
ts = 2 * sizeof(int) + MaxGenMethodLength;
memcpy(bsp, (void*) &(S->P.Vsize), ts);
bsp += ts;
    
/* Optimalitytype */
memcpy(bsp, (void*) &(S->optimality), sizeof(OPTIMALITYTYPE));
}


/*---------------------------------------------------------------------*/


void BitStringToSolution(SOLUTION* S, void* bs)
/* Puts the data from a bitstring to a solution. Space for S
   should be allocated. Needs changes when SolutionToBitString
   is changed. */
{
int i;
int ts;
void* bsp;

bsp = bs;
/* Codebook */
ts = MaxVersionLength + 8 * sizeof(int) + MaxGenMethodLength;
memcpy((void*) &(S->CB), bsp, ts);
bsp += ts;
ts = VectorSize(&(S->CB)) * sizeof(VECTORELEMENT);
for (i=0; i<S->CB.AllocatedSize; i++)
    {
    memcpy((void*) Node(&(S->CB),i).vector, bsp, ts);
    bsp += ts;
    memcpy((void*) &(Node(&(S->CB),i).vmean), bsp, 2 * sizeof(int));
    bsp += 2 * sizeof(int);
    }
memcpy((void*) &(S->CB.AllocatedSize), bsp, sizeof(int));
bsp += sizeof(int);

/* Partitioning */
ts = MaxVersionLength + 2 * sizeof(int);
memcpy((void*) &(S->P), bsp, ts);
bsp += ts;
ts = S->P.TSsize * sizeof(int);
memcpy((void*) S->P.Map, bsp, ts);
bsp += ts;
memcpy((void*) S->P.First, bsp, ts);
bsp += ts;
memcpy((void*) S->P.Next, bsp, ts);
bsp += ts;
ts = S->P.Vsize * sizeof(llong);
for (i=0; i<S->P.PartitionCount; i++)
    {
    memcpy((void*) S->P.CC[i].counter, bsp, ts);
    bsp += ts;
    memcpy((void*) &(S->P.CC[i].freq), bsp, sizeof(int));
    bsp += sizeof(int);
    }
ts = S->P.PartitionCount * sizeof(int);
memcpy((void*) S->P.Uniques, bsp, ts);
bsp += ts;
ts = 2 * sizeof(int) + MaxGenMethodLength;
memcpy((void*) &(S->P.Vsize), bsp, ts);
bsp += ts;

/* Optimalitytype */
memcpy((void*) &(S->optimality), bsp, sizeof(OPTIMALITYTYPE));

/* Set namepointers to NULLs */
for (i=0; i<S->CB.AllocatedSize; i++)
    VectorName(&(S->CB), i) = NULL;
}


/*---------------------------------------------------------------------*/


int SizeOfSolution(SOLUTION* S)
/* Returns the length of the bitstring of SolutionBitString. */
{
return ( 2 * MaxVersionLength +
(13 + S->CB.AllocatedSize + 3 * S->P.TSsize + 2 * S->P.PartitionCount)
    * sizeof(int) +
2 * MaxGenMethodLength +
( VectorSize(&(S->CB)) + 1 ) * sizeof(VECTORELEMENT) * S->CB.AllocatedSize +
S->P.Vsize * sizeof(llong) * S->P.PartitionCount +
sizeof(OPTIMALITYTYPE) );
}


/*-------------------  GLA related routines  --------------------------*/


void ChangeOptimality(TRAININGSET*   TS,
                      SOLUTION*      S,
                      OPTIMALITYTYPE optimality,
                      ERRORFTYPE     errorf)
{
  if (S->optimality == optimality) return;
  if (S->optimality == OPT_BOTH)   return;

  switch (optimality)
    {
    case OPT_NONE:
      break;
    case OPT_CB:
      {
      GenerateOptimalCodebook(TS, &S->CB, &S->P);
      S->optimality = OPT_CB;
      break;
      }
    case OPT_PA:
      {
      GenerateOptimalPartitioningGeneral(TS, &S->CB, &S->P, errorf);
      S->optimality = OPT_PA;
      break;
      }
    case OPT_BOTH:
      {
      IterateGLAForSolution(TS, S, MAX_GLA, errorf);
      S->optimality = OPT_BOTH;
      break;
      }
    }
}
           
    
/*---------------------------------------------------------------------*/


static void CheckChangedCodeVectors(SOLUTION*    S,          
                                    CODEBOOK*    oldCB,      
                                    int*         newcount,   
                                    int*         newindices) 
{
  int i;
  
  *newcount = 0;
                                     
  for( i = 0; i < BookSize(&S->CB); i++ )
    {
    if (!(EqualVectors(Vector(&S->CB,i), Vector(oldCB,i), VectorSize(&S->CB))))
      {
      newindices[*newcount] = i;
      (*newcount)++;
      }
    }
}


/*---------------------------------------------------------------------*/


static void NewSolutionByLocalRepartitioning(TRAININGSET* TS,
                                             SOLUTION*    S,
                                             int          newcount,
                                             int*         newindices,
                                             DISTANCETYPE disttype)
{     
  int     i, j, best;
  llong   olddist, newdist;
  YESNO   changed;
  
  for ( i=0; i<BookSize(TS); i++ )
    {
    best = -1;
    olddist = VectorDistance(Vector(&S->CB,Map(&S->P,i)),  
                             Vector(TS,i), 
                             VectorSize(TS),
                             MAXLLONG,
                             disttype);
    changed = NO;
    for (j=0; j<newcount; j++)
      if ( Map(&S->P,i) == newindices[j] )
        {
        changed = YES;
        break;
        }
    
    if (changed)
      {
      for ( j=0; j<BookSize(&S->CB); j++)
        {
        newdist = VectorDistance(Vector(&S->CB,j),  
                                 Vector(TS,i), 
                                 VectorSize(TS),
                                 olddist,
                                 disttype);

        if ( newdist < olddist )
          {
          best    = j;
          olddist = newdist;
          }
        }
      }
    else
      {
      for ( j=0; j<newcount; j++)
        {
        newdist = VectorDistance(Vector(&S->CB,newindices[j]),  
                                 Vector(TS,i), 
                                 VectorSize(TS),
                                 olddist,
                                 disttype);
        if ( newdist < olddist )
          {
          best    = newindices[j];
          olddist = newdist;
          }
        }
      }
    
    if ( best > -1 )
      ChangePartition(TS, &S->P, best, i);
    }
}      


/*---------------------------------------------------------------------*/


void IterateGLAForSolution(TRAININGSET* TS, 
                           SOLUTION*    S, 
                           int          count,
                           ERRORFTYPE   errorf)
{
  int      round;
  int*     newindices;
  int      newcount = 1;
  CODEBOOK cb;

  if (count < 1) return;
  
  switch (S->optimality)
    {
    case OPT_NONE:        /* start by generating an optimal partitioning */
    case OPT_CB:
      {
      ChangeOptimality(TS, S, OPT_PA, errorf);
      IterateGLAForSolution(TS, S, count-1, errorf);
      ChangeOptimality(TS, S, OPT_CB, errorf);
      FillEmptyPartitions(TS, &S->CB, &S->P);
      break;
      }
    case OPT_PA:
      {
      newindices = allocate (BookSize(&S->CB) * sizeof(int) );
      CreateNewCodebook(&cb, BookSize(&S->CB), TS);
      round = 0;

      while (round < count && newcount > 0)
        {
        CopyCodebook(&S->CB, &cb);
        GenerateOptimalCodebook(TS, &S->CB, &S->P);
        FillEmptyPartitions(TS, &S->CB, &S->P);
        CheckChangedCodeVectors(S, &cb, &newcount, newindices);
        NewSolutionByLocalRepartitioning(TS, S, newcount, newindices,
                                         DistType(errorf));
        round++;
        }
      deallocate(newindices);
      FreeCodebook(&cb);
      break;
      }
    case OPT_BOTH:
      break;
    }      
}          


/*--------------------------------------------------------------------*/


void SortSolution(TRAININGSET* TS, SOLUTION* S)
/* Sorts solution so that TS[0] belongs to cluster 0 etc. */
{
  int i;
  int *newindices;               /* cluster i will become newindices[i] */
  int latest = 0;
  SOLUTION newS;

  newindices = allocate(BookSize(&S->CB) * sizeof(int));
  for (i=0; i<BookSize(&S->CB); i++)
    newindices[i] = -1;

  for (i=0; i<BookSize(TS); i++)
    if (newindices[Map(&S->P, i)] < 0)
      newindices[Map(&S->P, i)] = latest++;

  if (latest < BookSize(&S->CB)) /* There were empty clusters */
    for (i=0; i<BookSize(&S->CB); i++)
      if (newindices[i] < 0) newindices[i] = latest++;

  CreateNewSolution(TS, &newS, BookSize(&S->CB));
  for (i=0; i<BookSize(&S->CB); i++)
    CopyNode(&Node(&S->CB, i), &Node(&(newS.CB), newindices[i]), 
	     VectorSize(&S->CB));
  for (i=0; i<BookSize(TS); i++)
    ChangePartition(TS, &(newS.P), newindices[Map(&S->P, i)], i);
  CopySolution(&newS, S);
  FreeSolution(&newS);

  deallocate(newindices);
}
