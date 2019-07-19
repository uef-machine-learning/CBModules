#if ! defined(__SA_H)
#define __SA_H

#include "cb.h"

/* ----------------------------------------------------------------- */

typedef  struct { int    Distribution;
                  int    DecreaseFunction;
                  double Change;
                  double CurrentTemperature;
                  int    ApplyToCodevectors;
                  int    ApplyToTrainingvectors;
                } SASchedule;

/* ----------------------------------------------------------------- */

void        InitializeSASchedule(int          Distribution,
                                 int          DecreaseFunction,
                                 double       change,
                                 double       initial,
                                 int          ToCodevectors,
                                 int          ToTrainingvectors,
                                 SASchedule* t);
void        DecreaseTemperature(SASchedule* t);
BOOKNODE*   RandomizeVectorBySA(SASchedule* t,
                                BOOKNODE*   source,
                                BOOKNODE*   dest,
                                int         Vsize,
                                int         maxvalue);


#define SASInUse(SAS)   ((SAS) == NULL ? 0 : ((SAS)->ApplyToCodevectors || (SAS)->ApplyToTrainingvectors) )
#define SASUseToCB(SAS) ((SAS) == NULL ? 0 : (SAS)->ApplyToCodevectors)
#define SASUseToTS(SAS) ((SAS) == NULL ? 0 : (SAS)->ApplyToTrainingvectors)
#define SASTemp(SAS)    ((SAS) == NULL ? 0 : (SAS)->CurrentTemperature)
#define SASEffective(SAS,Vsize,maxvalue)                   \
          (SASInUse(SAS) ? (SASTemp(SAS) >= 0.5) : 0)

/* Future implementation
#define SASEffective(SAS,Vsize,maxvalue)                   \
          (SASInUse(SAS) ? (maxvalue) == 1 ?               \
                           (SASTemp(SAS)*(Vsize) >= 0.5) : \
                           (SASTemp(SAS)         >= 0.5)   \
                         : 0)
*/

/* ----------------------------------------------------------------- */

#endif /* __SA_H */

