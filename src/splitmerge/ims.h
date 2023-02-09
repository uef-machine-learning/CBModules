#if ! defined(__IMS_H)
#define __IMS_H

/* ----------------------------------------------------------------- */

#include "cb.h"
#include "sa.h"

/* ----------------------------------------------------------------- */

void SetIMSParameters(int ShowProgress,
                      int CodebookSize,
                      int AlternatingOrder,
                      int AlternatingAmount,
                      int MinCodebookSize,
                      int MaxCodebookSize,
                      int StepSizeChange,
                      int Iterations,
                      int IMSGLAIterations,
                      int RequireDistortionDecrease);
void IMS(TRAININGSET*  TS,
         CODEBOOK*     CB,
         SASchedule*   SAS,
         PARTITIONING* P,
         int*          iteration);

/* ----------------------------------------------------------------- */

#endif /* __IMS_H */


