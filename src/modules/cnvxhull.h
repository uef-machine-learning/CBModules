#ifndef __CNVXHULL_H
#define __CNVXHULL_H

#include "cb.h"

/* ----------------------------------------------------------------- */

struct CONVEXVERTEXSTRUCT
         {
         int                        index; /* Index of the training vector */
         struct CONVEXVERTEXSTRUCT* next;  /* Next element */
         };

typedef struct CONVEXVERTEXSTRUCT  CONVEXVERTEX;
typedef struct
          {
          int           size;
          CONVEXVERTEX* first;
          } CONVEXHULL;
typedef struct
          {
          int         size;
          CONVEXHULL* convexhulls;
          } CONVEXHULLSET;

/* ----------------------------------------------------------------- */

#define CV_index(CV)         ((CV)->index)
#define CV_next(CV)          ((CV)->next)
#define CH_size(CH)          ((CH)->size)
#define CH_first(CH)         ((CH)->first)
#define CHS_size(CHS)        ((CHS)->size)
#define CHS_hull(CHS,i)      (&((CHS)->convexhulls[i]))
#define CHS_hullsize(CHS,i)  (CH_size(CHS_hull((CHS),(i))))
#define CHS_hullfirst(CHS,i) (CH_first(CHS_hull((CHS),(i))))


/* ----------------------------------------------------------------- */

void ConstructConvexHulls(TRAININGSET*    TS,
                          PARTITIONING*   P,
                          CONVEXHULLSET** CHS);
void FreeConvexHulls(CONVEXHULLSET* CHS);

/* ----------------------------------------------------------------- */

#endif /* __CNVXHULL_H */

