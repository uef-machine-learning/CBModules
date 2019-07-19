#if ! defined(__CONTREE_H)
#define __CONTREE_H

#include "owntypes.h"

#define  FormatNameCT  "ct"

typedef enum { DIVIDED, TMP_DIVIDED, FINAL, OLDLEAF, NEWLEAF } NODESTATUS;

struct MODELTREE_NODE { long                     Whites;
                        long                     Total;
                        double                   CumDynEntropy;
                        struct MODELTREE_NODE*   WhiteChild;
                        struct MODELTREE_NODE*   BlackChild;
                        NODESTATUS               Status;
                      };

struct CONTREE_NODE   { long                     ContNumber;
                        struct CONTREE_NODE*     WhiteChild;
                        struct CONTREE_NODE*     BlackChild;
                        NODESTATUS               Status;
                      };

typedef struct MODELTREE_NODE   MODELTREE_TYPE;
typedef MODELTREE_TYPE*  MODELTREE;

typedef struct CONTREE_NODE  CONTREE_TYPE;
typedef CONTREE_TYPE*        CONTREE;

typedef struct { int MinDepth;
                 int StartDepth;
                 int MaxDepth;
                 int FinalDepth;
                 float MinLearnCost;
                 float MaxLearnCost;
                 float ConstLearnCost;
                 int CTModelling;
                 int QMModelling;
                 float Growth;
                 int CostFunction;
                 int QuietLevel;
                 int MemoryGiven;
                 long NodesTried;
                 long NodesExpanded;
                 long ContextsTotal;
               } CONTREEDATA;

/* ------------------ export functions -------------------- */

CONTREE MakeContextTree (IMAGE* image, CONTREEDATA* Data, int QuietLevel);
void DeallocateContextTree(CONTREE T);

YESNO ReadCTHeader(FILE* f, int* fulldepth, long* contexts);
YESNO WriteCTHeader(FILE* f, int fulldepth, long contexts);

CONTREE  ReadContextTree(FILE* f, int fulldepth, long* contextsread);
YESNO  StoreContextTree(CONTREE C, FILE* f, int fulldepth);

long GetVariableContext(CONTREE C, IMAGE* image, int x, int y);


#endif /* __CONTREE_H */
