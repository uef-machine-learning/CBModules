#if ! defined(__GLA_H)
#define __GLA_H

/*--------------------  Basic type definitions -----------------------*/

#include "cb.h"
#include "sa.h"


void SetGLAParameters(int ShowProgress,
                      int GLAMethod,
                      int UsePDS,
                      int PDSInitialGuess,
                      int CodevectorCalculation,
                      ERRORFTYPE errorf);

void GLA(TRAININGSET*  TS,
         CODEBOOK*     CB,
         SASchedule*   SAS,         /* NULL allowed */
         PARTITIONING* OrigP,       /* NULL allowed */
         int           IterLimit,
         double*       InitError,
         double*       FinalError,
         int*          Iterations);

/* ----------------------------------------------------------------- */

#endif /* __GLA_H */

