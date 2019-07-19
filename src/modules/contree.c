/*-------------------------------------------------------------------*/
/* CONTREE.C      Eugene Ageenko                                     */
/*                Pasi Fr„nti                                        */
/*                                                                   */
/* CONTEXT TREE (Variable size context) MODELLING                    */
/*                                                                   */
/* .CONTEXT TREE TYPES ARE SUPPORTED:                                */
/*   = STATIC                                                        */
/*   = SEMIADAPTIVE (cost of tree expansion is taken into conside-   */
/*                   ration when constructed)                        */
/*                                                                   */
/* .CONTEXT TREE BUILDING STRATEGIES ARE IMPLEMETED:                 */
/*   = TOP-DOWN, level by level construction/pruning                 */
/*   = BOTTOM-UP one pass pruning, (full tree constructed before)    */
/*   - HYBRID. (first BOTTOM-UP, then continue with TOP-DOWN)        */
/*                                                                   */
/* .TREE PRUNING STRATEGY                                            */
/*   - IF NO SAVING observed if contexts is divided then prune!      */
/*   - SAVING calculated as OLDBITS - NEWBITS - COST.                */
/*   - COST is 2 for S/A CT, and 8 for the S/A JBIG in additional    */
/*                                                                   */
/* .COST FUNCTION IS OPTIMIZED FOR:                                  */
/*   = STATIC modelling (for normal JBIG)                            */
/*     - saving  is calculated on the basis of the final statistics  */
/*     - artificial LEARNING COST is added and                       */
/*     - controlled by ADJUSTMENT parameter (0-150%)                 */
/*   = SEMI-ADAPTIVE modelling                                       */
/*     - like previous but the MODELLING cost is taken in the con-   */
/*       sideration, (for the S/A JBIG or EDM), no L.C. of course.   */
/*   = DYNAMIC modelling (for normal JBIG, like in normal JBIG)      */
/*     - saving  is calculated as the cumulative sum for each pixel  */
/*       on the basis of current (previous) statistics               */
/*                                                                   */
/* .SPECIAL TRICKS FOR THE BETTER CT CONSTRUCTION ARE:               */
/*  = ONE LEVEL DELAYED PRUNING to avoid LOCAL INOPTIMALITY problem  */
/*    (applicable for TOP-DOWN pruning, for STATIC or DYNAMIC m-ing) */
/*     L.I. problem - when get no improvement on the context, due    */
/*     to high learning cost, but next division can provide us with  */
/*     huge improvement)                                             */
/*  = GREEDY EXPANSION                                               */
/*    (applicable for TOP-DOWN pruning, for STATIC modelling)        */
/*    - first tree expansion w/o Learning Cost (adjustment is Zero)  */
/*    - then pruning tree BOTTOM-UP with the real Learning Cost      */
/*  = OPTIMAL MEMORY                                                 */
/*    (applicable for BOTTOM-UP pruning, for DYNAMIC modelling)      */
/*    - strategy is combinational of the sequential BOTTOM-UP        */
/*      techniques applied one by another to fit the given memory    */
/*                                                                   */
/* .OTHER FEATURES AND COMMENTS                                      */
/*   - only tree structure is written to the file ".CT", starting    */
/*     from the given (minimal) full tree depth                      */
/*   - HYBRID tree building strategy is outlined:                    */
/*       1)  Build tree to the given starting depth                  */
/*       2)  Prune it up to the given minimal depth                  */
/*       3)  Add and prune new levels to the tree                    */
/*           one-by-one (new pass over the image)                    */
/*                                                                   */
/*-------------------------------------------------------------------*/

#define ProgName        "CONTREE"
#define VersionNumber   "0.14"
#define LastUpdated     "21.11.97"

/* ------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>

#include "binimage.h"
#include "file.h"
#include "memctrl.h"
#include "error.h"
#include "contree.h"

/* -------------------------- Misc -----------------------------------*/

#if defined(max) || defined(min)
#undef max
#undef min
#endif
#define  max(a,b) ((a) > (b) ? (a) : (b))
#define  min(a,b) ((a) < (b) ? (a) : (b))
#define  ScaleBetween(a,b,c)   max((a), min((b), (c)))

#if defined(round)
#undef round
#endif
#define round(a)      ( (int) ((a) + 0.5) )

#define  NODEMODELLINGCOST 8
#define  NODEDIVIDINGCOST 2
#define  GREEDYSAVINGS 25

typedef enum {SILENT, QUIET, NORMAL, PROGRESS, STAT, FULL_STAT, DEBUG_STAT} QUIETMODE;
typedef enum {NOTPRUNED, PRUNED} PRUNETYPE;
typedef enum {STATIC=0, SEMIADAPTIVE=1, ADAPTIVE=2} QM_MODELLINGTYPE;
typedef enum {NO_CT=0, STATIC_CT=1, SEMIADAPTIVE_CT=2, ADAPTIVE_CT=3} CT_MODELLINGTYPE;
typedef enum {GEZERO=0, DELAYED1=1, GREEDY=2, OPTIMALMEMORY=3} COSTFUNCTIONTYPE;
typedef enum {TOP_DOWN, BOTTOM_UP} PRUNEDIRECTION;

long CN = 0; /* Context Number */
int  EstimatedCost = 0;
long NodesAdded = 0;             /* added at the stage */
long NodesPruned = 0;            /* pruned at all at the stage */
long NewNodesPruned = 0;         /* pruned on the last level at the stage */
long CurrentLastLevel = 0;       /* current last level */
double OldBitsCurrent = 0.0;
double NewBitsCurrent = 0.0;
double NewBitsCurrent0 = 0.0;
double NewBitsCurrent1 = 0.0;
double BitsSaved = 0.0;
double BitsOld = 0.0;

static CONTREEDATA* Data;

static int  TemplateX[31] = {30,-1,0,-1,+1,-2,0,-2,+2,-1,+1,-2,+2,-3,0,-3,+3,-1,+1,-3,+3,-2,+2,-4, 0,-4,+4,-1,+1,-3,+3 };
static int  TemplateY[31] = {30,0,-1,-1,-1,0,-2,-1,-1,-2,-2,-2,-2,0,-3,-1,-1,-3,-3,-2,-2,-3,-3, 0,-4,-1,-1,-4,-4,-3,-3 };

/*    JBIG extended template   */
/*           27 24 28          */
/*     29 21 17 14 18 22 30    */
/*     19 11  9  6 10 12 20    */
/*  25 15  7  3  2  4  8 16 26 */
/*  23 13  5  1  ?             */
/*                             */

static double ScewLearningCost[] =

{
0.0,         0.000288558, 0.001154464, 0.002598411, 0.004621561,
0.007225546, 0.010412479, 0.014184963, 0.018546105, 0.023499531,
0.029049406, 0.035200451, 0.041957978, 0.049327907, 0.057316811,
0.065931945, 0.075181295, 0.085109755, 0.096199585, 0.108615584,
0.122380702, 0.137520008, 0.15406088,  0.172033228, 0.191469756,
0.212406252, 0.235270162, 0.260507465, 0.288189241, 0.318393905,
0.351208215, 0.38718513,  0.426913162, 0.470552605, 0.518517223,
0.571794566, 0.630781486, 0.696666714, 0.770378905, 0.853510218,
0.947881055, 1.055928948, 1.181030827, 1.327860107, 1.503334641,
1.718158701, 1.990349949, 2.353752578, 2.884555243, 3.826547107,
5.0
};

/*

p(0.995) = 4.79501567

{  0.0,
   0.028569152, 0.056583528, 0.084064265, 0.111031312, 0.137503524,
   0.163498732, 0.189033824, 0.214124805, 0.238786860, 0.263034406,
   0.286881148, 0.310340121, 0.333423734, 0.356143810, 0.378511623,
   0.400537930, 0.429428502, 0.472175804, 0.514299035, 0.555816155,
   0.596744360, 0.637100124, 0.676899239, 0.716156852, 0.754887502,
   0.812213971, 0.868791053, 0.924638087, 0.979773675, 1.034215715,
   1.105903347, 1.176711354, 1.246661061, 1.327268666, 1.412636233,
   1.501806312, 1.601878759, 1.709028002, 1.824942931, 1.953900254,
   2.098816123, 2.261725930, 2.446288999, 2.660989929, 2.915917047,
   3.230407652, 3.638206997, 4.215685252, 5.0        , 5.0};
*/

/*======================= Q M T - Interface ===========================*/


YESNO ReadCTHeader(FILE* f, int* fulldepth, long* contexts)
/* return Minimal Full Depth or 0 */
{
  int Id1, Id2;

  Id1 = getc(f);
  Id2 = getc(f);
  if( Id1 != 'C' || Id2 != 'T' )
  {
     ErrorMsg ("CT identification flag expected. Possible not a CONTEXT-TREE-FILE.\n");
     return(NO);
     /* not .QT file or broken */
  }
  fscanf(f, "%i %li", fulldepth, contexts);
  /* Skip white space; only one character can appear! */
  fgetc(f);
  return(YES);
}


/*-------------------------------------------------------------------*/


YESNO WriteCTHeader(FILE* f, int fulldepth, long contexts)
{
  putc('C', f);
  putc('T', f);
  putc(10, f);
  fprintf(f, "%i %li", fulldepth, contexts);  putc(10, f);
  return (YES);
}


/* ========================== MODEL TREE =========================== */


MODELTREE ContextAtLevel (MODELTREE M, IMAGE* image, int x, int y, int level)
{
   int i;

   for (i=0;i<level;i++)
   {
      if (M->Status == FINAL) return (NULL);
      /* we assume here that tree is full otherwise node has FINAL status */
      if (GetImageBitPixel(image, x+TemplateX[i+1], y+TemplateY[i+1] ) == WHITE)
         M = M->WhiteChild;
      else
         M = M->BlackChild;
   }
   return (M);
}


/* -------------------------- BUILDING ----------------------------- */


MODELTREE NewModelTreeNode()
{
  MODELTREE M;

  M = allocate(sizeof(MODELTREE_TYPE));
  M ->Status = NEWLEAF;
  M ->WhiteChild = NULL;
  M ->BlackChild = NULL;
  M ->CumDynEntropy = 0.0;
  if (Data->QMModelling == ADAPTIVE)
  {
     M ->Whites = 1L;
     M ->Total = 2L;
  }
  else
  {
     M ->Whites = 0L;
     M ->Total = 0L;
  }

  return(M);
}


/* ----------------------------------------------------------------- */

void ConstructLevel(MODELTREE M, int level)
/* level - new level number, root is zero */
{

  if (M->Status == FINAL) return;
  if (level==1)
  {
     M ->WhiteChild = NewModelTreeNode();
     M ->BlackChild = NewModelTreeNode();
     M ->Status = DIVIDED;
     NodesAdded++;
  }
  else
  {
     /* Are the IF statements necessary? The first IF ..=FINAL should be enough */
     if (M->WhiteChild != NULL) ConstructLevel(M->WhiteChild, level-1);
     if (M->BlackChild != NULL) ConstructLevel(M->BlackChild, level-1);
  }

}


/* ----------------------------------------------------------------- */


MODELTREE ConstructModelTree(int depth)
{
  MODELTREE M;

  M = NewModelTreeNode();
  if (depth > 0)
  {
     M ->WhiteChild = ConstructModelTree(depth-1);
     M ->BlackChild = ConstructModelTree(depth-1);
     M ->Status = DIVIDED;
  }

  return(M);
}


/* ------------------------ CALCULATING ---------------------------- */


int EstimateCost(int CTModelling, int QMModelling)
{
  int cost=0;

  if (QMModelling == SEMIADAPTIVE) cost += NODEMODELLINGCOST;
  if (CTModelling == SEMIADAPTIVE_CT) cost += NODEDIVIDINGCOST;
  return (cost);

}


/* ----------------------------------------------------------------- */


double Entropy(double p)
{
  return( p==0 ? 0 : p==1 ? 0 : - (p*log(p) + (1-p)*log(1-p)) / log(2) );
}


/* ----------------------------------------------------------------- */


double LearningCost(long Whites, long Total)
{
  double p,cost;
  int lookup;

  if (Total == 0) return (0.0);
  p = (double) Whites / Total;
  if (p<0.5) p = 1 - p;
  lookup = round((p - 0.5) * 100);

  cost = ScewLearningCost[lookup];
  if (lookup == 49)
     { if ( ( round((p - 0.990) * 200) ) == 1 ) cost = 4.79501567; }

  cost = max(cost,Data->MinLearnCost);
  cost = min(cost,Data->MaxLearnCost);

  return (cost);
}

/* ----------------------------------------------------------------- */


double EstimateLearningCost(MODELTREE M)
{
  double cost;
  if (Data->ConstLearnCost != 0) return (Data->ConstLearnCost);
  cost = LearningCost(M->WhiteChild->Whites,M->WhiteChild->Total) +
         LearningCost(M->BlackChild->Whites,M->BlackChild->Total) -
         LearningCost(M->Whites,M->Total);
  return (cost);
}

/* ----------------------------------------------------------------- */


double CalculateSavingIfDivided(MODELTREE M)
{
   long total1, total0;
   double pw, pw1, pw0, saving;

   pw = M->Total ? (double) M->Whites / M->Total : 0;
   OldBitsCurrent = Entropy(pw) * M->Total;

   total0 = M->WhiteChild->Total;
   total1 = M->BlackChild->Total;
   pw0 = total0 ? ((double) M->WhiteChild->Whites) / total0 : 0;
   pw1 = total1 ? ((double) M->BlackChild->Whites) / total1 : 0;
   NewBitsCurrent0 = Entropy(pw0) * total0;
   NewBitsCurrent1 = Entropy(pw1) * total1;
   NewBitsCurrent = NewBitsCurrent0 + NewBitsCurrent1;

   saving = OldBitsCurrent - NewBitsCurrent - EstimatedCost;
   return(saving);

}


/* ----------------------------------------------------------------- */


void CalculateTree(MODELTREE MT, IMAGE* image, int mindepth, int depth)
/* calculate data for level from mindepth to depth */

{

  MODELTREE M;
  int y,x,i,Pel,Pel_i,ImageSizeY,ImageSizeX;
  long amount;

  ImageSizeY=image->ImageSizeY; ImageSizeX=image->ImageSizeX;
  for (y=1;y<=ImageSizeY;y++)
      for (x=1;x<=ImageSizeX;x++)
      {
         M = ContextAtLevel(MT, image, x, y, mindepth);
         if (M==NULL) continue;
         Pel = GetImageBitPixel(image, x, y);

         for (i=mindepth;i<depth;i++)
         {
            if (M->Status == FINAL) { M=NULL; break; }
            Pel_i = GetImageBitPixel(image, x+TemplateX[i+1], y+TemplateY[i+1] );

            if (Data->QMModelling == ADAPTIVE)
            {
               amount = (Pel == WHITE ? M->Whites : M->Total - M->Whites);
               if (M->Total != 0L && amount != 0L) M->CumDynEntropy -=  log((double)amount / M->Total);
            }

            if (Pel == WHITE)  { M->Whites++; M->Total++; }
            else                 M->Total++;

            if (Pel_i == WHITE)  { M = M->WhiteChild; }
            else                 { M = M->BlackChild; }
         }

         if (M!=NULL)
         {
            if (Data->QMModelling == ADAPTIVE)
            {
               amount = (Pel == WHITE ? M->Whites : M->Total - M->Whites);
               if (M->Total != 0L && amount != 0L) M->CumDynEntropy -= log((double)amount / M->Total);
            }
            if (Pel == WHITE)  { M->Whites++; M->Total++; }
            else                 M->Total++;
         }

      }

}


/* ----------------------------------------------------------------- */


void CalculateTotalBitsStatic(MODELTREE M, int level)
{
   double pw;

   /* We find leaves not belove given level and calculate static entropy */
   if (level <= 0 || M->WhiteChild == NULL || M->BlackChild == NULL)
   {
      pw = M->Total ? (double) M->Whites / M->Total : 0;
      BitsOld += Entropy(pw) * M->Total;
   }
   else
   {
      CalculateTotalBitsStatic(M->WhiteChild, level-1);
      CalculateTotalBitsStatic(M->BlackChild, level-1);
   }
}


/* ----------------------------------------------------------------- */


void CalculateTotalBitsDynamic(MODELTREE M, int level)
{
   /* We find leaves not belove given level and calculate dynamic entropy */
   if (level <= 0 || M->WhiteChild == NULL || M->BlackChild == NULL)
   {
      BitsOld += M->CumDynEntropy / log(2);
   }
   else
   {
      CalculateTotalBitsDynamic(M->WhiteChild, level-1);
      CalculateTotalBitsDynamic(M->BlackChild, level-1);
   }
}


/* -------------------------- PRUNING ------------------------------ */


void PrintNodeStatistic(MODELTREE M, double S, double LC, char* path,
                        PRUNEDIRECTION direction, char* decision)
{
   if ((Data->QuietLevel>=FULL_STAT) ||
       (Data->QuietLevel>=STAT && direction== TOP_DOWN && decision[0]!='-') ||
       (Data->QuietLevel>=STAT && direction== BOTTOM_UP && decision[0]!='+'))
   {
      printf ("%-25s ",path);
      printf ("%12.2f -> %12.2f, save= %10.2f ", OldBitsCurrent, NewBitsCurrent, S);
      printf ("(%s)\n",decision);
      if (Data->QuietLevel>=DEBUG_STAT)
      {
         printf("%li/%li (%.2f) -> %li/%li (%.2f), %li/%li (%.2f), LC=%.2f\n\n",
         M->Whites,M->Total, OldBitsCurrent,
         M->WhiteChild->Whites,M->WhiteChild->Total, NewBitsCurrent0,
         M->BlackChild->Whites,M->BlackChild->Total, NewBitsCurrent1, LC);
      }
   }
}


/* ----------------------------------------------------------------- */


PRUNETYPE PruneNode (MODELTREE M, MODELTREE Parent,
                     NODESTATUS status_if_pruned, char* path)

{
   double S,LC=0.0;

   if (Data->QMModelling == ADAPTIVE)
   {
      OldBitsCurrent  = M->CumDynEntropy / log(2);
      NewBitsCurrent = (M->WhiteChild->CumDynEntropy +
                        M->BlackChild->CumDynEntropy) / log(2);

      if (Data->QuietLevel>=DEBUG_STAT)
      {
         NewBitsCurrent0 = M->WhiteChild->CumDynEntropy / log(2);
         NewBitsCurrent1 = M->BlackChild->CumDynEntropy / log(2);
      }
      S = OldBitsCurrent - NewBitsCurrent - EstimatedCost;
   }
   else
   {
      S = CalculateSavingIfDivided(M);
      if (Data->QMModelling == STATIC)
      {
         LC=EstimateLearningCost(M);
         S-=LC*(Data->Growth);
      }
   }

   if (S > 0)
   {
      M->Status = DIVIDED; /* by default */
      M->WhiteChild->Status = OLDLEAF;
      M->BlackChild->Status = OLDLEAF;
      BitsSaved += S;
      PrintNodeStatistic(M, S, LC, path,
        (status_if_pruned==FINAL ? TOP_DOWN : BOTTOM_UP), "+");

      /* we can check the parent status here and turn it back to the DIVIDED */
      if ( (Data->CostFunction == DELAYED1 &&
            status_if_pruned == FINAL) && /* i.e. we go down */
           Parent->Status == TMP_DIVIDED) Parent->Status = DIVIDED;

      return(NOTPRUNED);

   }
   else
   {
      if ( (Data->CostFunction == DELAYED1 &&
            status_if_pruned == FINAL) &&    /* i.e. we go down */
           CurrentLastLevel < Data->MaxDepth &&
           Parent->Status != TMP_DIVIDED)   /* parent isn't TMP_DIVIDED */
      {
         M->Status = TMP_DIVIDED;
         M->WhiteChild->Status = OLDLEAF;
         M->BlackChild->Status = OLDLEAF;
         /* BitsSaved += S; */  /* ??? what to do here? we ca not count */
         PrintNodeStatistic(M, S, LC, path,
            (status_if_pruned==FINAL ? TOP_DOWN : BOTTOM_UP), "?");
         return(NOTPRUNED);
      }
      else
      {
         M->Status = status_if_pruned;
         NewNodesPruned++; NodesPruned++;
         PrintNodeStatistic(M, S, LC, path,
            (status_if_pruned==FINAL ? TOP_DOWN : BOTTOM_UP), "-");
         deallocate(M->WhiteChild); M->WhiteChild = NULL;
         deallocate(M->BlackChild); M->BlackChild = NULL;
         return(PRUNED);
      }
   }

}


/* ----------------------------------------------------------------- */


PRUNETYPE PruneLevel(MODELTREE M, MODELTREE Parent, int level, char* path)
/* level to prune (parents) */
{
  char newpath[35];
  int wc_pruned,bc_pruned;

  if (M->Status == FINAL) return (NOTPRUNED);

  if (level > 0)
  {
     strcpy(newpath,path);
     strcat(newpath,"0");
     wc_pruned = PruneLevel(M->WhiteChild, M, level-1, newpath);
     strcpy(newpath,path);
     strcat(newpath,"1");
     bc_pruned = PruneLevel(M->BlackChild, M, level-1, newpath);

     if (wc_pruned == PRUNED && bc_pruned == PRUNED && M->Status == TMP_DIVIDED)
     {
        /* both cildren were pruned, then this node, which was candidate */
        /* for pruning (negative saving), must be pruned too             */
        M->Status = FINAL;
        NodesPruned++;
        PrintNodeStatistic(M, 0, 0, path, TOP_DOWN, "--");
        deallocate(M->WhiteChild); M->WhiteChild = NULL;
        deallocate(M->BlackChild); M->BlackChild = NULL;
        /* Oops! Here we must correct BitsSaved, but we can not :(       */
        /* So ... real BitsSaved will be slightly less due to ovedue.    */
     }
     return(NOTPRUNED);
  }
  else
  {
     return(PruneNode(M,Parent,FINAL,path));
  }

}


/* ----------------------------------------------------------------- */


void PostProcessTree(MODELTREE M, int level)
/* convert all LEAVES which are before STARTDEPTH to FINAL */
{

  if (level > 0)
  {
     if (M->WhiteChild == NULL || M->BlackChild == NULL)
     {
        M->Status = FINAL;
     }
     else
     {
        PostProcessTree(M->WhiteChild, level-1);
        PostProcessTree(M->BlackChild, level-1);
     }
  }

}


/* ----------------------------------------------------------------- */


void PruneTree(MODELTREE M, MODELTREE Parent, int depth, int levels, char* path)
/* bottom-up recursive pruning, *levels* to prune up with given */
{
  char newpath[35];

  if (levels<0) return; /* process is skipped */

  if (M->Status == FINAL) return;

  /* first go down */

  if (depth > 0)
  {
     strcpy(newpath,path);
     strcat(newpath,"0");
     PruneTree(M->WhiteChild, M, depth-1, levels, newpath);
     strcpy(newpath,path);
     strcat(newpath,"1");
     PruneTree(M->BlackChild, M, depth-1, levels, newpath);
  }

  /* when go back, prune the node, bottom nodes are just pruned before */

  if (depth <= levels)
  {

     if ( (M->WhiteChild->Status != NEWLEAF && M->WhiteChild->Status != OLDLEAF) ||
          (M->BlackChild->Status != NEWLEAF && M->BlackChild->Status != OLDLEAF) )
     /* only nodes, where both children are NEWLEAF can be pruned */
     {
        if (M->WhiteChild->Status == NEWLEAF || M->WhiteChild->Status == OLDLEAF) M->WhiteChild->Status = FINAL;
        if (M->BlackChild->Status == NEWLEAF || M->BlackChild->Status == OLDLEAF) M->BlackChild->Status = FINAL;
        return;
     }
     else /* actually pruning */
     {
        PruneNode(M,Parent,OLDLEAF,path);
     }

  }

}


/* --------------------- DEALLOCATING ------------------------------ */


void DeallocateModelTree(MODELTREE T)
{
   if (T->WhiteChild != NULL) DeallocateModelTree(T->WhiteChild);
   if (T->BlackChild != NULL) DeallocateModelTree(T->BlackChild);
   deallocate (T);
}


/* ====================== CONTEXT TREE ============================= */


long GetVariableContext(CONTREE C, IMAGE* image, int x, int y)
{
  int i = 1;

  while (C->Status == DIVIDED || C->Status == TMP_DIVIDED)
  {
     if (GetImageBitPixel(image, x+TemplateX[i], y+TemplateY[i] ) == WHITE)
        C = C->WhiteChild;
     else
        C = C->BlackChild;
     i++;
  }

  return (C->ContNumber);

}


/* ------------------------ BUILDING ------------------------------- */


CONTREE NewContextTreeNode()
{
   CONTREE C;

   C = allocate(sizeof(CONTREE_TYPE));
   C->ContNumber = 0L;
   C->WhiteChild = NULL;
   C->BlackChild = NULL;
   C->Status = FINAL;
   return (C);
}


/* ------------------------------------------------------------------*/


CONTREE ConstructContextTree(MODELTREE M)
{
   CONTREE C;

   C = NewContextTreeNode();

   if (M->Status == DIVIDED || M->Status == TMP_DIVIDED)
   {
      C->Status = DIVIDED;
      C-> WhiteChild = ConstructContextTree(M->WhiteChild);
      C-> BlackChild = ConstructContextTree(M->BlackChild);
   }
   else /* == FINAL */
   {
      C->Status = FINAL;
      C->ContNumber = CN;
      CN++;
   }

   return (C);
}


/*-------------------- STORE / READ ---------------------------------*/


YESNO  StoreContextTree_(CONTREE C, BITSTREAM* bt, int fulldepth)
{

  if (fulldepth > 0)
  {
     assert(C->Status != FINAL);
     StoreContextTree_(C->WhiteChild, bt, fulldepth-1);
     StoreContextTree_(C->BlackChild, bt, fulldepth-1);
  }
  else if (C->Status == FINAL)
  {
     OutputBit(bt,0);
  }
  else
  {
     OutputBit(bt,1);
     StoreContextTree_(C->WhiteChild, bt, 0);
     StoreContextTree_(C->BlackChild, bt, 0);
  }

  return(YES);
}


YESNO  StoreContextTree(CONTREE C, FILE* f, int fulldepth)
{
  BITSTREAM bs;

  InitializeBitStream(&bs, f);
  StoreContextTree_(C, &bs, fulldepth);
  FlushOutput(&bs);
  return (YES);
}


/*-------------------------------------------------------------------*/


CONTREE ReadContextTree_(BITSTREAM* bt, int fulldepth)
{

  CONTREE C;

  C = NewContextTreeNode();
  if (fulldepth > 0)
  {
     C->Status = DIVIDED;
     C->WhiteChild=ReadContextTree_(bt, fulldepth-1);
     C->BlackChild=ReadContextTree_(bt, fulldepth-1);
  }
  else if (InputBit(bt)==1)
  {
     C->Status = DIVIDED;
     C->WhiteChild=ReadContextTree_(bt, 0);
     C->BlackChild=ReadContextTree_(bt, 0);
  }
  else /* bit = 0 */
  {
     C->Status = FINAL;
     C->ContNumber = CN;
     CN++;
  }

  return (C);
}


CONTREE  ReadContextTree(FILE* f, int fulldepth, long* contextsread)
{
  CONTREE C;
  BITSTREAM bs;
  CN=0L;
  InitializeBitStream(&bs, f);
  C=ReadContextTree_(&bs, fulldepth);
  FlushInput(&bs);
  *contextsread = CN;
  return (C);
}


/* --------------------- DEALLOCATING ------------------------------ */


void DeallocateContextTree(CONTREE T)
{
   if (T->WhiteChild != NULL) DeallocateContextTree(T->WhiteChild);
   if (T->BlackChild != NULL) DeallocateContextTree(T->BlackChild);
   deallocate (T);
}


/* --------------------- DRAWING ----------------------------------- */


void PrintContextTree(CONTREE T)
{



}


/* ================================================================= */
/*                      PRIMARY FUNCTION                             */
/* ================================================================= */


CONTREE MakeContextTree (IMAGE* image, CONTREEDATA* UserData, int QL)
{
  CONTREE  CT;
  MODELTREE  MT;
  int k,LCfactor = 0;
  long improvement,newimprovement;

  Data = UserData;
  EstimatedCost = EstimateCost(Data->CTModelling, Data->QMModelling);
  BitsSaved = 0;
  BitsOld = 0;
  NodesAdded = 0;
  NodesPruned = 0;
  NewNodesPruned = 0;
  CurrentLastLevel = Data->StartDepth;
  Data->NodesTried = 0;
  Data->NodesExpanded = 0;
  Data->ContextsTotal = 0;
  newimprovement = 1;

  if (Data->QMModelling == STATIC && Data->CostFunction == GREEDY)
  {
     LCfactor = Data->Growth;
     Data->Growth = 0;
  }
  /* GREEDY aproach affects only STATIC sheme, it means no LCost DOWN */
  /* with Learning Cost UP                                            */


  if (QL>2) printf ("Construncting the model tree (depth=%i). ",Data->StartDepth);
  MT = ConstructModelTree(Data->StartDepth);      /* full tree with startdepth  */

  if (QL>2) printf ("Calculating from %i to %i\n", Data->MinDepth, Data->StartDepth);
  CalculateTree(MT, image, Data->MinDepth, Data->StartDepth);  /* all values for last levels */

  /*
  if (Data->QMModelling == ADAPTIVE)
        CalculateTotalBitsDynamic(MT,Data->MinDepth);
  else  CalculateTotalBitsStatic(MT,Data->MinDepth);
  */

  if (Data->QMModelling != STATIC || Data->CostFunction != GREEDY)
  {
    if (QL>2) printf ("Pruning tree bottom-up: ");
    if (QL>3) printf ("\n");
    PruneTree(MT,NULL,Data->StartDepth-1,Data->StartDepth-1-Data->MinDepth,"");
    PostProcessTree(MT, Data->StartDepth);
    if (QL>2) printf ("%li nodes were pruned.\n\n",NodesPruned);

    Data->NodesExpanded -= NodesPruned;
  }

  k = Data->StartDepth;
  do {

     k++;
     CurrentLastLevel = k;
     if (k>Data->MaxDepth) break; /* we exceed the template */
     NodesAdded = 0;
     NodesPruned = 0;
     NewNodesPruned = 0;
     BitsSaved = 0;

     if (QL>2) printf ("Constructing level %i: ",k);
     ConstructLevel(MT, k);
     Data->NodesTried += NodesAdded;
     if (QL>2) printf ("%li nodes were expanded. ",NodesAdded);
     if (QL>2) printf ("Calculating.\n");
/*     CalculateModels(MT, k);  */
     ImageRewind(image);

     CalculateTree(MT, image, k, k);
     newimprovement = NodesAdded;
     improvement = NodesAdded;
     Data->NodesExpanded += NodesAdded;

     /* if GREEDY then not necessary to coarse-prune the last level */
     if (Data->CostFunction == GREEDY && k==Data->MaxDepth)
        {if (QL>3) printf ("\n"); continue;}
     /* it will breaks automaticly */

     if (QL>2) printf ("Pruning: ");
     if (QL>3) printf ("\n");
     PruneLevel(MT, NULL, k-1, "");
     newimprovement -= NewNodesPruned;
     improvement -= NodesPruned;
     Data->NodesExpanded -= NodesPruned;
     if (QL>2) printf ("%li new nodes pruned, %li nodes expanded, %.1f bit saved\n\n",NodesPruned,newimprovement,BitsSaved);

  } while (newimprovement > 0); /* UNTIL no improvement */

  /* correction of the MaxDepth */
  Data->MaxDepth = k-1;

  /* correction of the MinDepth might be implemented in the PruneTree() */

  if (Data->QMModelling == STATIC && Data->CostFunction == GREEDY)
  {
    Data->Growth = LCfactor;

    if (QL>2) printf ("Fine pruning tree bottom-up: ");
    if (QL>3) printf ("\n");
    ImageRewind(image);

    PruneTree(MT,NULL,Data->MaxDepth-1,Data->MaxDepth-1-Data->MinDepth,"");
    PostProcessTree(MT, Data->MaxDepth);
    if (QL>2) printf ("%li nodes were pruned.\n\n",NodesPruned);

    Data->NodesExpanded -= NodesPruned;
  }

  Data->FinalDepth    = k-1;
  Data->ContextsTotal = pow(2,Data->StartDepth) + Data->NodesExpanded;

  CN = 0L;
  CT = ConstructContextTree(MT);
  if (CN != Data->ContextsTotal) {ErrorMsg("CONTREE: Error constructing the Context Tree.\n %li contexts were expected, %li contexts were constructed.\n",Data->ContextsTotal,CN); CT=NULL;}
  DeallocateModelTree(MT);
  return (CT);
}


/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */

/* Now we can CompressImage(CT, k); */
/* Compression (to EDM)

CompressImage(CT, image):

FOR y=1 TO ImageSizeY
        FOR x=1 TO ImageSizeX
!!!              c = GetVariableContext (CT, image, x, y);
                        CompressPixel( Pixel[x,y], c);

*/

