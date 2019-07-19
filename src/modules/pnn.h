#if ! defined(__PNN_H)
#define __PNN_H

/* ----------------------------------------------------------------- */

#include "cb.h"
#include "sa.h"

/* ----------------------------------------------------------------- */

void SetPNNParameters(int ShowProgress,
                      int MergeErrorType,
                      int CentroidCalculation,
                      int MergeMethod,
                      int PartitionRemapping,
                      int GLAIterations,
                      int MergePercent,
                      int MaxBucket,
                      int MinBucket);
void PairwiseNearestNeighbour(TRAININGSET*  TS,
                              CODEBOOK*     CB,
                              SASchedule*   SAS,
                              PARTITIONING* P,
                              int           ResultCBSize);

/* ----------------------------------------------------------------- */

#endif /* __PNN_H */

