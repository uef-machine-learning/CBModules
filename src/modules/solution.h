#if ! defined(__SOLUTION_H)
#define __SOLUTION_H

#include "owntypes.h"
#include "cb.h"

/*---------------------------- Definitions ---------------------------*/

typedef  enum { OPT_NONE=0, OPT_CB=1, OPT_PA=2, OPT_BOTH=3 } OPTIMALITYTYPE;

/*---------------------  Solution data structure  --------------------*/

typedef  struct { CODEBOOK       CB;
                  PARTITIONING   P;  
                  OPTIMALITYTYPE optimality; } SOLUTION; 

/*--------------------  Solution interface  --------------------------*/

CBFILETYPE ReadSolution(TRAININGSET* TS, SOLUTION* S, char* FileName);
void CreateNewSolution(TRAININGSET* TS, SOLUTION* S, int booksize);
void FreeSolution(SOLUTION* S);
void IncreaseSolutionSize(SOLUTION* S, int newsize);
void CopySolution(SOLUTION* sourceS, SOLUTION* destS);

void ChangeCodeVectorInSolution(TRAININGSET* TS, 
                                SOLUTION*    S, 
                                int          item, 
                                VECTORTYPE   v,
                                ERRORFTYPE   errorf,
                                YESNO        dolocalrep);
void ChangePartitionInSolution (TRAININGSET* TS,
                                SOLUTION*    S,
                                int          item,
                                int          Pindex,
                                ERRORFTYPE   errorf,
                                YESNO        updatecentroids);

void SolutionToBitString(SOLUTION* S, void* bs);
void BitStringToSolution(SOLUTION* S, void* bs);
int  SizeOfSolution(SOLUTION* S);
void SortSolution(TRAININGSET* TS, SOLUTION* S);

/*-------------------  GLA & error related routines  --------------------*/

void ChangeOptimality(TRAININGSET*   TS,
                      SOLUTION*      S,
                      OPTIMALITYTYPE optimality,
                      ERRORFTYPE     errorf);

void IterateGLAForSolution(TRAININGSET* TS, 
                           SOLUTION*    S, 
                           int          count,
                           ERRORFTYPE   errorf);

#define SolutionError(TS, S, errorf) \
        AverageErrorForSolution(TS, &((S)->CB), &((S)->P), errorf)


/*--------------------------------------------------------------------*/

#endif /* __SOLUTION_H */
