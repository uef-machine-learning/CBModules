#if ! defined(__SPLIT_H)
#define __SPLIT_H

/* ----------------------------------------------------------------- */

#include "cb.h"
#include "sa.h"

/*--------------------  Search tree data structure  ------------------*/
                            
struct STNODESTRUCT { VECTORTYPE  Centroid;
                      YESNO       Leaf;
                      int         ClusterIndex;  /* Only for leafs */
                      struct STNODESTRUCT* Left;
                      struct STNODESTRUCT* Right; };
typedef struct STNODESTRUCT   SEARCHTREE;

/* -------------------  Search tree routines  ----------------------- */

SEARCHTREE* GenerateSearchTreeForCodebook(CODEBOOK* CB);

/* ----------------------------------------------------------------- */

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
                        int GLAIterations);

void Split(TRAININGSET*  TS,
           CODEBOOK*     CB,
           SASchedule*   SAS,
           PARTITIONING* P,
           int           ResultCBSize);

/* ----------------------------------------------------------------- */

#endif /* __SPLIT_H */

