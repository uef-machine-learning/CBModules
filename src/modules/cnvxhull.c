
/*-------------------------------------------------------------------*/
/* CNVXHULL.C    Timo Kaukoranta                                     */
/*               Rami Kaila                                          */
/*                                                                   */
/* - 2-dimensional convex hull                                       */
/*                                                                   */
/*                                                                   */
/*-------------------------------------------------------------------*/

#define  ProgName       "CNVXHULL"
#define  VersionNumber  "version 0.02a"
#define  LastUpdated    "27.10.99"        /* iak */

/* ----------------------------------------------------------------- */

#include <math.h>

#include "cnvxhull.h"
#include "memctrl.h"

/* ----------------------------------------------------------------- */

#define iSquare(a) ( ((double)a)*((double)a) )
#define dSquare(a) ( (a)*(a) )
#define vlength(x1,y1,x2,y2) sqrt( iSquare((x1)-(x2)) + iSquare((y1)-(y2)) )


/* Procedures for Convex Hull ===============================================*/

static CONVEXVERTEX* AddVertex(CONVEXVERTEX* last, int TSindex)
{
  /* Allocate a new element for the linked list */
  CV_next(last) = (CONVEXVERTEX*)allocate(sizeof(CONVEXVERTEX));
  CV_index(CV_next(last)) = TSindex;
  CV_next(CV_next(last)) = NULL;

  return( CV_next(last) );
}

/* ----------------------------------------------------------------- */

static void FindFirstAndSecondConvexVertex(TRAININGSET*  TS,
                                           PARTITIONING* P,
                                           int           Pindex,
                                           int*          first,
                                           int*          second)
{
  int TSindex;
  int x, y;
  int firstx, firsty;
  int secondx, secondy;
  int dx, dy;
  int seconddx, seconddy;

  /* Search two first convex vertex:
     The first one is the leftmost vector.
     The second one has greatest slope from the first one.
  */

  /* Search the leftmost vector (smallest x component). */
  firstx = firsty = TS->MaxValue + 1;
  for( TSindex = FirstVector(P, Pindex);
       !EndOfPartition(TSindex);
       TSindex = NextVector(P, TSindex) )
    {
    x = VectorScalar(TS, TSindex, 0);
    y = VectorScalar(TS, TSindex, 1);
    if( x < firstx || (x == firstx && y < firsty) )
      {
      *first = TSindex;
      firstx = x;
      firsty = y;
      }
    }

  /* Search the vector V for which the slope
     (Vy - firsty)/(Vx -firstx) is greatest.
  */
  TSindex = FirstVector(P, Pindex);
  if( TSindex == *first )
    {
    TSindex = NextVector(P, TSindex);
    }
  *second = TSindex;
  secondx = VectorScalar(TS, TSindex, 0);
  secondy = VectorScalar(TS, TSindex, 1);
  seconddx = secondx - firstx;
  seconddy = secondy - firsty;
  TSindex = NextVector(P, TSindex);

  while( !EndOfPartition(TSindex) )
    {
    if( TSindex != *first )
      {
      x = VectorScalar(TS, TSindex, 0);
      y = VectorScalar(TS, TSindex, 1);

      if( seconddx == 0 )
        {
        if( firstx == x && firsty < y )
          {
          *second = TSindex;
          secondx = x;
          secondy = y;
          seconddx = secondx - firstx;
          seconddy = secondy - firsty;
          }
        }
      else /* seconddx > 0 */
        {
        dx = x - firstx;
        dy = y - firsty;
        if( dx == 0 ||
            ((double)dy/(double)dx > (double)seconddy/(double)seconddx) )
          {
          *second = TSindex;
          secondx = x;
          secondy = y;
          seconddx = secondx - firstx;
          seconddy = secondy - firsty;
          }
        }
      }

    TSindex = NextVector(P, TSindex);
    }
}

/*---------------------------------------------------------------------------*/

static void FindNextConvexVertex(TRAININGSET*  TS,
                                 PARTITIONING* P,
                                 int           Pindex,
                                 int           previous,
                                 int           current,
                                 int*          next)
{
  int    TSindex;
  double MinCosAngle, CosAngle;
  double a, b, c, MaxDist;
  int    prevx, prevy, currx, curry, x, y;

  MinCosAngle = 1;
  MaxDist = 0.0;
  prevx = VectorScalar(TS, previous, 0);
  prevy = VectorScalar(TS, previous, 1);
  currx = VectorScalar(TS, current, 0);
  curry = VectorScalar(TS, current, 1);

  TSindex = FirstVector(P, Pindex);
  while( !EndOfPartition(TSindex) )
    {
    if( TSindex != current )
      {
      x = VectorScalar(TS, TSindex, 0);
      y = VectorScalar(TS, TSindex, 1);

      a = vlength(prevx, prevy, x, y);
      b = vlength(prevx, prevy, currx, curry);
      c = vlength(currx, curry, x, y);

      CosAngle = -(dSquare(a) - dSquare(b) - dSquare(c)) / (2.0 * b * c);
      if( CosAngle < MinCosAngle ||
          (CosAngle == MinCosAngle && c > MaxDist) )
        {
        MinCosAngle = CosAngle;
        MaxDist = c;
        *next = TSindex;
        }
      }
    TSindex = NextVector(P, TSindex);
    }
}

/*---------------------------------------------------------------------------*/

static void ConstructConvexHull(TRAININGSET*  TS,
                                PARTITIONING* P,
                                CONVEXHULL*   CH,
                                int           Pindex)
{
  CONVEXVERTEX* last = CH_first(CH);
  int start, previous, current, next;

  switch( UniqueVectors(P, Pindex) )
    {
    case 0: /* There is no training vectors */
      {
      CH_size(CH) = 0;
      break;
      }
    case 1: /* There is only one training vector */
      {
      last = AddVertex(last, FirstVector(P, Pindex));
      CH_size(CH) = 1;
      break;
      }
    case 2: /* There is only two training vectors */
      {
      last = AddVertex(last, FirstVector(P, Pindex));
      last = AddVertex(last, NextVector(P, FirstVector(P, Pindex)));
      CH_size(CH) = 2;
      break;
      }
    default: /* There is more than two training vectors */
      {
      FindFirstAndSecondConvexVertex(TS, P, Pindex, &previous, &current);
      last = AddVertex(last, previous);
      last = AddVertex(last, current);
      CH_size(CH) = 2;

      /* Find the training vector which has greatest angle between
         two previous training vectors and has longest distance (=c)
         againts it using cosine clause
         a2 = b2 + c2 - 2 * a * b * cos alfa
         -> Alfa = -(a2 - b2 - c2) / (2 * b * c) and
            Alfa <= 1
         a = (fx,fy),(ax,ay) b = (cx,cy),(fx,fy) c = (cx,cy),(ax,ay)
      */

      /* do this until the hull is drawn to starting place */
      start = previous;
      do
        {
        FindNextConvexVertex(TS, P, Pindex, previous, current, &next);
        last = AddVertex(last, next);
        CH_size(CH)++;
        previous = current;
        current = next;
        } while( current != start );
      break;
      }
    }
}


/*---------------------------------------------------------------------------*/

void ConstructConvexHulls(TRAININGSET*    TS,
                          PARTITIONING*   P,
                          CONVEXHULLSET** CHS)
{
  int Pindex;
  CONVEXVERTEX* dummy;

  /* Allocate memory for set header */
  *CHS = (CONVEXHULLSET*)allocate(sizeof(CONVEXHULLSET));
  CHS_size(*CHS) = PartitionCount(P);

  /* Allocate memory for hull headers */
  (*CHS)->convexhulls = (CONVEXHULL*)allocate(CHS_size(*CHS) * sizeof(CONVEXHULL));
  for( Pindex = 0; Pindex < CHS_size(*CHS); Pindex++ )
    {
    CHS_hullsize(*CHS, Pindex) = 0;
    /* Allocate dummy element for the linked list */
    CHS_hullfirst(*CHS, Pindex) = (CONVEXVERTEX*)allocate(sizeof(CONVEXVERTEX));
    CV_index(CHS_hullfirst(*CHS, Pindex)) = -1; /* Unnecessary */
    CV_next(CHS_hullfirst(*CHS, Pindex)) = NULL;
    }

  /* For each partition construct a convex hull */
  for( Pindex = 0; Pindex < PartitionCount(P); Pindex++ )
    {
    ConstructConvexHull(TS, P, CHS_hull(*CHS, Pindex), Pindex);

    /* Remove the dummy vertex. 'ConstructConvexHull' handles
       the size field.
    */
    dummy = CHS_hullfirst(*CHS, Pindex);
    CHS_hullfirst(*CHS, Pindex) = CV_next(dummy);
    deallocate(dummy);
    }
}

/* ----------------------------------------------------------------- */

void FreeConvexHulls(CONVEXHULLSET* CHS)
{
  int           Pindex;
  CONVEXVERTEX* CV;

  for( Pindex = 0; Pindex < CHS_size(CHS); Pindex++ )
    {
    while( CHS_hullfirst(CHS, Pindex) != NULL )
      {
      CV = CV_next(CHS_hullfirst(CHS, Pindex));
      deallocate(CHS_hullfirst(CHS, Pindex));
      CHS_hullfirst(CHS, Pindex) = CV;
      }
    }

  deallocate(CHS->convexhulls);
  deallocate(CHS);
}

