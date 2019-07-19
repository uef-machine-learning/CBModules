#if ! defined(__SEARCH_H)
#define __SEARCH_H

/* ----------------------------------------------------------------- */

#include "split.h"

/*------------------------ Type definitions --------------------------*/

#define  FormatNameSS    "ss"

/*--------------------  Vector-mean search  --------------------------*/

int  FindNearestUsingVectorMean(BOOKNODE* v, CODEBOOK* CB, int* tested);
void GenerateOptimalPartitionUsingVmeans(TRAININGSET*  TS,
                                         CODEBOOK*     CB,
                                         PARTITIONING* P);

/*--------------------  PCA-tree interface  --------------------------*/

SEARCHTREE* GenerateSearchTree(CODEBOOK* CB);
void        FreeSearchTree(SEARCHTREE* ST);
void        PrintSearchTree(SEARCHTREE* ST, int VectorSize);
int         SearchFromTree(SEARCHTREE* ST, BOOKNODE* node, int VectorSize);
double      RunGLATreeSearch(TRAININGSET*  TS,
                             CODEBOOK*     CB,
                             PARTITIONING* P,
                             int           iterations);

/*--------------------------------------------------------------------*/

#endif /* __SEARCH_H */
