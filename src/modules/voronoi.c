#include <assert.h>

#include <math.h>

#include <stdio.h>



#include "memctrl.h"

#include "sort.h"

#include "stack.h"

#include "voronoi.h"



/* ----------------------------------------------------------------- */



#define EPSILON    0.0000001

#define iszero(x)  (-EPSILON < (x) && (x) < EPSILON)



#define sqr(a)     ((a)*(a))



#ifndef M_PI

#define M_PI 3.1415926535897932384626433832795

#endif





/* ================================================================= */



typedef struct

          {

          int  size;

          int  allocatedsize;

          int* set;

          } IndexSet;



/* ----------------------------------------------------------------- */



static IndexSet* IS_make(int size)

{

  IndexSet* IS;



  assert( 0 <= size );



  IS = allocate(sizeof(IndexSet));

  IS->size = 0;

  IS->allocatedsize = size;

  IS->set = allocate(IS->allocatedsize * sizeof(int));



  return( IS );

}



/* ----------------------------------------------------------------- */



static void IS_free(IndexSet* IS)

{

  deallocate(IS->set);

  deallocate(IS);

}



/* ----------------------------------------------------------------- */



static int IS_size(IndexSet* IS)

{

  return( IS->size );

}



/* ----------------------------------------------------------------- */



static void IS_add(IndexSet* IS, int item)

{

  if( IS->size == IS->allocatedsize )

    {

    int* tmp = IS->set;

    int  i;



    IS->allocatedsize = 2 * IS->allocatedsize + 1;

    IS->set = allocate(IS->allocatedsize * sizeof(int));

    for( i = 0; i < IS->size; i++ )

      {

      IS->set[i] = tmp[i];

      }

    deallocate(tmp);

    }



  IS->set[IS->size] = item;

  IS->size++;

}



/* ----------------------------------------------------------------- */



static void IS_set(IndexSet* IS, int index, int item)

{

  assert( index < IS->size );



  IS->set[index] = item;

}



/* ----------------------------------------------------------------- */



static int IS_index(IndexSet* IS, int index)

{

  assert( index < IS->size );



  return( IS->set[index] );

}



/* ----------------------------------------------------------------- */



static void IS_clear(IndexSet* IS)

{

  IS->size = 0;

}



/* ================================================================= */



static void CO_copy(Coordinate* dest, Coordinate* source)

{

  dest->x = source->x;

  dest->y = source->y;

}



/* ----------------------------------------------------------------- */



static YESNO CO_equal(Coordinate* a, Coordinate* b)

{

  return( iszero(a->x - b->x) && iszero(a->y - b->y) ? YES : NO );

}



/* ================================================================= */



static void UV_clear(UnitVector* UV)

{

  UV->dx = UV->dy = 0.0;

}



/* ----------------------------------------------------------------- */



static void UV_copy(UnitVector* dest, UnitVector* source)

{

  dest->dx = source->dx;

  dest->dy = source->dy;

}



/* ----------------------------------------------------------------- */



static double VectorLength(double deltax, double deltay)

{

  return( sqrt(deltax * deltax + deltay * deltay) );

}



/* ----------------------------------------------------------------- */



static void UV_set(UnitVector* UV, double deltax, double deltay)

{

  double len = VectorLength(deltax, deltay);



  if( iszero(len) )

    {

    UV->dx = 0.0;

    UV->dy = 0.0;

    }

  else

    {

    UV->dx = deltax / len;

    UV->dy = deltay / len;

    }

}



/* ----------------------------------------------------------------- */



static YESNO UV_equal(UnitVector* a, UnitVector* b)

{

  return( iszero(a->dx - b->dx) && iszero(a->dy - b->dy));

}



/* ----------------------------------------------------------------- */



static YESNO UV_opposite(UnitVector* a, UnitVector* b)

{

  return( iszero(a->dx + b->dx) && iszero(a->dy + b->dy));

}



/* ----------------------------------------------------------------- */



static double UV_angle(UnitVector* a)

{

  if( 0 <= a->dx )

    {

    if( 0 <= a->dy )

      {

      return( asin(a->dy) );

      }

    else

      {

      return( 2.0 * M_PI - asin(-a->dy) );

      }

    }

  else

    {

    if( 0 <= a->dy )

      {

      return( M_PI - asin(a->dy) );

      }

    else

      {

      return( M_PI + asin(-a->dy) );

      }

    }

}



/* ----------------------------------------------------------------- */



static double UV_dotproduct(UnitVector* a, UnitVector* b)

{

  return( a->dx * b->dx + a->dy * b->dy );

}



/* ================================================================= */



static void PrintEdge(Edge* e, int i)

{

  printf("E%2i:%c vx=%2i,%2i nx=%2i,%2i cv=%2i,%2i dir=%7.3f,%7.3f\n",

         i,

         (E_status(e) == DIRECTION ? 'D' : E_status(e) == CONNECTION ? 'C' : 'B'),

         E_vertex(e,0), E_vertex(e,1),

         E_next(e,0), E_next(e,1),

         E_code(e,0), E_code(e,1),

         E_dirx(e), E_diry(e));

}



/* ----------------------------------------------------------------- */



static void PrintEdges(VoronoiGraph* VG)

{

  int i;



  for( i = 0; i < VG_nEdge(VG); i++ )

    {

    PrintEdge(&VG->edge[i], i);

    }

}



/* ----------------------------------------------------------------- */



static void PrintVertex(Vertex* v, int i)

{

  printf("V%2i:edge=%2i coord=%7.3f,%7.3f\n",

         i, VX_edge(v), VX_x(v), VX_y(v));

}



/* ----------------------------------------------------------------- */



static void PrintVertices(VoronoiGraph* VG)

{

  int i;



  for( i = 0; i < VG_nVertex(VG); i++ )

    {

    PrintVertex(&VG->vertex[i], i);

    }

}



/* ----------------------------------------------------------------- */



static void PrintGraph(VoronoiGraph* VG)

{

  PrintVertices(VG);

  PrintEdges(VG);

}



/* ================================================================= */



int VG_EdgeEnd(VoronoiGraph* VG, int vertex, int edge)

{

  if( VG->edge[edge].vertex[0] == vertex )

    {

    return( 0 );

    }

  else

    {

    if( VG->edge[edge].vertex[1] == vertex )

      {

      return( 1 );

      }

    else

      {

      assert( 0 );

      return( -1 );

      }

    }

}



/* ----------------------------------------------------------------- */



int VG_EdgeNext(VoronoiGraph* VG, int vertex, int edge)

{

  int end = VG_EdgeEnd(VG, vertex, edge);



  return( VG->edge[edge].next[end] );

}



/* ----------------------------------------------------------------- */



static int PrevEdge(VoronoiGraph* VG, int vertex, int edge)

{

  int i, prev;



  assert( vertex != ENDOFLINK );



  prev = ENDOFLINK;

  i = VG_Edge(VG, vertex);

  while( i != ENDOFLINK && i != edge )

    {

    prev = i;

    i = VG_EdgeNext(VG, vertex, i);

    }



  return( i == edge ? prev : ENDOFLINK );

}



/* ----------------------------------------------------------------- */



static int LastEdge(VoronoiGraph* VG, int vertex)

{

  int i, prev;



  prev = i = VG_Edge(VG, vertex);

  while( i != ENDOFLINK )

    {

    prev = i;

    i = VG_EdgeNext(VG, vertex, i);

    }

  return( prev );

}



/* ================================================================= */



VoronoiGraph* VG_make(CODEBOOK* CB)

{

  VoronoiGraph* VG;



  VG = allocate(sizeof(VoronoiGraph));

  VG->nvertex = 0;

  VG->nedge   = 0;



  /* VoronoiVektoriin varataan nelinkertainen tila CODEBOOkiin

     tilaan n„hden, koska voronoipisteit„ on enint„„n kaksi

     kertaa CODEBOOKin tila ja koska voronoivektoriin lis„t„„n my”s

     alueen rajapisteet niin tarvitaan viel„ kaksinkertainen m„„r„

     pisteit„. Esim. CODEBOOkissa 1 piste => VoronoiVektorissa 4

     pistett„, jotka ovat alueen kulmapisteit„. Varsinaisia

     voronoipisteit„ ei esinny. T„m„ esim on huonoin tapaus.*/



  VG->nallocvertex = 4*BookSize(CB);

  VG->nallocedge   = 4*BookSize(CB);

  VG->vertex = allocate(VG->nallocvertex * sizeof(Vertex));

  VG->edge   = allocate(VG->nallocedge * sizeof(Edge));



  return( VG );

} 



/* ----------------------------------------------------------------- */



void VG_kill(VoronoiGraph* VG)

{

  assert( VG != NULL );



  deallocate(VG->vertex);

  deallocate(VG->edge);

  deallocate(VG);

}



/* ----------------------------------------------------------------- */



int VG_addVertex(VoronoiGraph* VG, Coordinate* c)

{

  int position = VG_nVertex(VG);



  CO_copy(VG_VertexCoordinate(VG, position), c);

  VG->vertex[position].edge = ENDOFLINK;

  VG_nVertex(VG)++;



  return( position );

}



/* ----------------------------------------------------------------- */



int VG_addEdge(VoronoiGraph* VG, int vertex1, int vertex2)

{

  int   position = VG_nEdge(VG);

  Edge* E = &(VG->edge[position]);

  int last1, last2;



  E->status = DIRECTION;

  E->vertex[0] = vertex1;

  E->vertex[1] = vertex2;



  if( vertex1 != ENDOFLINK )

    {

    last1 = LastEdge(VG, vertex1);

    if( last1 != ENDOFLINK )

      {

      VG->edge[last1].next[VG_EdgeEnd(VG, vertex1, last1)] = position;

      }

    else

      {

      VG->vertex[vertex1].edge = position;

      }

    }



  if( vertex2 != ENDOFLINK )

    {

    last2 = LastEdge(VG, vertex2);

    if( last2 != ENDOFLINK )

      {

      VG->edge[last2].next[VG_EdgeEnd(VG, vertex2, last2)] = position;

      }

    else

      {

      VG->vertex[vertex2].edge = position;

      }

    }



  E->next[0] = ENDOFLINK;

  E->next[1] = ENDOFLINK;



  UV_clear(&(E->direction));



  E->codevector[0] = -1;

  E->codevector[1] = -1;



  VG_nEdge(VG)++;



  return( position );

}



/* ----------------------------------------------------------------- */



void VG_changeEdgeVertex(VoronoiGraph* VG, int edge, int end, int vertex)

{

  int curvertex = VG_Vertex(VG, edge, end);

  int tmpedge;

  int tmpend;

  int last;



  if( curvertex != ENDOFLINK )

    {

    tmpedge = PrevEdge(VG, curvertex, edge);



    if( tmpedge == ENDOFLINK )

      {

      assert( VG->vertex[curvertex].edge == edge );

      VG_Edge(VG, curvertex) = VG_EdgeNext(VG, curvertex, edge);

      }

    else

      {

      tmpend = VG_EdgeEnd(VG, curvertex, tmpedge);

      VG->edge[tmpedge].next[tmpend] = VG_EdgeNext(VG, curvertex, edge);

      }

    }



  if( vertex != ENDOFLINK )

    {

    last = LastEdge(VG, vertex);

    if( last != ENDOFLINK )

      {

      VG->edge[last].next[VG_EdgeEnd(VG, vertex, last)] = edge;

      }

    else

      {

      VG_Edge(VG, vertex) = edge;

      }

    }



  VG->edge[edge].next[end] = ENDOFLINK;

  VG->edge[edge].vertex[end] = vertex;

}



/* ----------------------------------------------------------------- */



void VG_removeVertex(VoronoiGraph* VG, int vertex)

{

  int i, j;

  int end;

  int last = VG_nVertex(VG) - 1;



  /* Set all edges connected to 'vertex' to point ENDOFLINK. */

  i = VG_Edge(VG, vertex);

  while( i != ENDOFLINK )

    {

    j = i;

    i = VG_EdgeNext(VG, vertex, i);



    end = VG_EdgeEnd(VG, vertex, j);

    VG_changeEdgeVertex(VG, j, end, ENDOFLINK);

    }



  if( vertex < last )

    {

    /* Move last vertex to 'vertex' and update edges. */



    VG->vertex[vertex].coord = VG->vertex[last].coord;



    i = VG_Edge(VG, last);

    while( i != ENDOFLINK )

      {

      j = i;

      i = VG_EdgeNext(VG, last, i);



      end = VG_EdgeEnd(VG, last, j);

      VG_changeEdgeVertex(VG, j, end, vertex);

      }

    }



  VG_nVertex(VG)--;

}



/* ----------------------------------------------------------------- */



void VG_removeEdge(VoronoiGraph* VG, int edge)

{

   int vertex1 = VG_Vertex(VG, edge, 0);

   int vertex2 = VG_Vertex(VG, edge, 1);

   int preve;

   int end;

   int last;



   /* Remove the next links. */

   if( vertex1 != ENDOFLINK )

     {

     preve = PrevEdge(VG, vertex1, edge);

     if( preve != ENDOFLINK )

       {

       end = VG_EdgeEnd(VG, vertex1, preve);

       VG->edge[preve].next[end] = VG_EdgeNext(VG, vertex1, edge);

       }

     else

       {

       VG->vertex[vertex1].edge = VG_EdgeNext(VG, vertex1, edge);

       }

     }



   if( vertex2 != ENDOFLINK )

     {

     preve = PrevEdge(VG, vertex2, edge);

     if( preve != ENDOFLINK )

       {

       end = VG_EdgeEnd(VG, vertex2, preve);

       VG->edge[preve].next[end] = VG_EdgeNext(VG, vertex2, edge);

       }

     else

       {

       VG->vertex[vertex2].edge = VG_EdgeNext(VG, vertex2, edge);

       }

     }



   /* Move the last edge to the place of 'edge'. */

   last = VG_nEdge(VG) - 1;



   if( last != edge )

     {

     vertex1 = VG_Vertex(VG, last, 0);

     vertex2 = VG_Vertex(VG, last, 1);



     if( vertex1 != ENDOFLINK )

       {

       preve = PrevEdge(VG, vertex1, last);

       if( preve != ENDOFLINK )

         {

         end = VG_EdgeEnd(VG, vertex1, preve);

         VG->edge[preve].next[end] = edge;

         }

       else

         {

         VG->vertex[vertex1].edge = edge;

         }

       }



     if( vertex2 != ENDOFLINK )

       {

       preve = PrevEdge(VG, vertex2, last);

       if( preve != ENDOFLINK )

         {

         end = VG_EdgeEnd(VG, vertex2, preve);

         VG->edge[preve].next[end] = edge;

         }

       else

         {

         VG->vertex[vertex2].edge = edge;

         }

       }

     VG->edge[edge] = VG->edge[last];

     }



   VG->nedge--;

}



/* ----------------------------------------------------------------- */



void VG_setEdgeStatus(VoronoiGraph* VG, int edge, EdgeStatus status)

{

  VG->edge[edge].status = status;

}



/* ----------------------------------------------------------------- */



void VG_setEdgeDir(VoronoiGraph* VG, int edge, UnitVector* dir)

{

  UV_copy(&(VG->edge[edge].direction), dir);

}



/* ----------------------------------------------------------------- */



void VG_setEdgeCodes(VoronoiGraph* VG, int edge, int code1, int code2)

{

  VG->edge[edge].codevector[0] = code1;

  VG->edge[edge].codevector[1] = code2;

}



/* ----------------------------------------------------------------- */



void VG_turnEdgeAround(VoronoiGraph* VG, int edge)

{

  Edge* E = &VG->edge[edge];

  int tmp;



  tmp            = E_vertex(E, 0);

  E_vertex(E, 0) = E_vertex(E, 1);

  E_vertex(E, 1) = tmp;



  tmp            = E_next(E, 0);

  E_next(E, 0)   = E_next(E, 1);

  E_next(E, 1)   = tmp;



  tmp            = E_code(E, 0);

  E_code(E, 0)   = E_code(E, 1);

  E_code(E, 1)   = tmp;



  E_dirx(E) = -E_dirx(E);

  E_diry(E) = -E_diry(E);

}



/* ----------------------------------------------------------------- */



int VG_Vertex(VoronoiGraph* VG, int edge, int end)

{

  return( VG->edge[edge].vertex[end] );

}



/* ----------------------------------------------------------------- */



EdgeStatus VG_EdgeStatus(VoronoiGraph* VG, int edge)

{

  return( VG->edge[edge].status );

}



/* ----------------------------------------------------------------- */



UnitVector* VG_EdgeDir(VoronoiGraph* VG, int edge)

{

  return( &(VG->edge[edge].direction) );

}



/* ----------------------------------------------------------------- */



int VG_EdgeCode(VoronoiGraph* VG, int edge, int end)

{

  return( VG->edge[edge].codevector[end] );

}



/* ================================================================= */



typedef struct

          {

          int code1;

          int code2;

          int prevvor;

          int code3;

          } Triple;



/* ----------------------------------------------------------------- */



static void TripleToStack(STACK* S,

                          int    code1,

                          int    code2,

                          int    prevvor,

                          int    code3)

{

  Triple* item;



/*   printf("TTS %2i %2i %2i %2i\n", code1, code2, prevvor, code3); */



  item = allocate(sizeof(Triple));

  if( code1 < code2 )

    {

    item->code1   = code1;

    item->code2   = code2;

    }

  else

    {

    item->code1   = code2;

    item->code2   = code1;

    }

  item->prevvor = prevvor;

  item->code3   = code3;

  S_push(S, (void*)item);

}



/* ----------------------------------------------------------------- */



static void TripleFromStack(STACK* S,

                            int*   code1,

                            int*   code2,

                            int*   prevvor,

                            int*   code3)

{

  Triple* item;



  item = (Triple*)S_pop(S);

  *code1   = item->code1;

  *code2   = item->code2;

  *prevvor = item->prevvor;

  *code3   = item->code3;

  deallocate(item);



/* printf("TFS %2i %2i %2i %2i\n", *code1, *code2, *prevvor, *code3); */

}



/* ----------------------------------------------------------------- */



static void StatesToStack(STACK* S, IndexSet* Formers, int prevvor)

{

  int i;



  for( i = 1; i < IS_size(Formers) - 1; i++ )

    {

    /* Last paremeter???? */

    TripleToStack(S, IS_index(Formers, i), IS_index(Formers, i+1), prevvor, 0);

    }

}



/* ================================================================= */



static int cmpCodeVectors(const void* a, const void* b, const void* info)

/* Ascending order assumed */

{

  CODEBOOK* CB = (CODEBOOK*)info;

  int       i  = *(int*)a;

  int       j  = *(int*)b;



  if( VectorScalar(CB, i, 0) < VectorScalar(CB, j, 0) )

    {

    return( 1 );

    }

  else

    {

    if( VectorScalar(CB, i, 0) == VectorScalar(CB, j, 0) &&

        VectorScalar(CB, i, 1) <  VectorScalar(CB, j, 1) )

      {

      return( 1 );

      }

    else

      {

      return( 0 );

      }

    }

}



/* ----------------------------------------------------------------- */



static void GenerateOrder(CODEBOOK* CB, int** order)

{

  int i;



  *order = (int*)allocate(BookSize(CB) * sizeof(int));

  for(i = 0; i < BookSize(CB); i++)

    {

    (*order)[i] = i;

    }



  QuickSort(*order, BookSize(CB), sizeof(int), (void*)CB, cmpCodeVectors);

}



/* ----------------------------------------------------------------- */



static void FreeOrder(int* order)

{

  deallocate(order);

}



/* ================================================================= */



static YESNO OnSameLine(CODEBOOK* CB,

                        int*      order,

                        int       code1,

                        int       code2,

                        int       code3)

  {

  UnitVector UV12, UV13;



  UV_set(&UV12,

         VectorScalar(CB,order[code2],0) - VectorScalar(CB,order[code1],0),

         VectorScalar(CB,order[code2],1) - VectorScalar(CB,order[code1],1));

  UV_set(&UV13,

         VectorScalar(CB,order[code3],0) - VectorScalar(CB,order[code1],0),

         VectorScalar(CB,order[code3],1) - VectorScalar(CB,order[code1],1));



  if( UV_equal(&UV12, &UV13) || UV_opposite(&UV12, &UV13) )

     {

     return YES;

     }

   else

     {

     return NO;

     }

  }



/* ----------------------------------------------------------------- */



static void CalculateVertex(CODEBOOK*   CB,

                            int*        order,

                            int         code1,

                            int         code2,

                            int         code3,

                            Coordinate* coord)

  {

  /* Lasketaan kahden eri pisteparin muodostamalle janalle keskinormaalit.

     N„iden normaalien leikkauspisteess„ on voronoipiste */

  Coordinate cent12, cent23, cent13;

  VECTORTYPE v1 = Vector(CB, order[code1]);

  VECTORTYPE v2 = Vector(CB, order[code2]);

  VECTORTYPE v3 = Vector(CB, order[code3]);

  double     slope12 = 1.0;

  double     slope23 = 1.0;



  cent12.x = (v1[0]+v2[0]) / 2.0;

  cent23.x = (v2[0]+v3[0]) / 2.0;

  cent13.x = (v1[0]+v3[0]) / 2.0;



  cent12.y = (v1[1]+v2[1]) / 2.0;

  cent23.y = (v2[1]+v3[1]) / 2.0;

  cent13.y = (v1[1]+v3[1]) / 2.0;



  if( v2[0] != v1[0] )

    {

    slope12 = (v2[1]-v1[1]) / (double)(v2[0]-v1[0]);

    }

  if(v3[0] != v2[0])

    {

    slope23 = (v3[1]-v2[1]) / (double)(v3[0]-v2[0]);

    }



  /* Riipuen siit„ mink„ suuntainen keskinormaali on, niin saadaan

     eri tapauksissa laskettu eri tavalla voronoipiste (periaate

     on sama kaikissa, mutta toisten tapausten muodostamista voidaan

     kiert„„ */



  if( iszero(slope12) && ( v3[0] == v2[0] ) )

    {

    coord->x = cent12.x;

    coord->y = cent23.y;

    }

  else if( ( v2[0] == v1[0] ) && iszero(slope23) )

    {

    coord->x = cent23.x;

    coord->y = cent12.y;

    }

  else if( iszero(slope12) )

    {

    coord->x = cent12.x;

    coord->y = -(coord->x-cent23.x)/slope23 + cent23.y;

    }

  else if( iszero(slope23) )

    {

    coord->x = cent23.x;

    coord->y = -(coord->x-cent12.x)/slope12 + cent12.y;

    }

  else if(v2[0] == v1[0] )

    {

    coord->y = cent12.y;

    coord->x = -(coord->y-cent23.y-cent23.x/slope23)*slope23;

    }

  else if(v3[0] == v2[0])

    {

    coord->y = cent23.y;

    coord->x = -(coord->y-cent12.y-cent12.x/slope12)*slope12;

    }

  else

    {

    assert( !iszero(slope12-slope23) );

    assert( !iszero(slope12) );



    coord->x = ( slope12 * slope23 * (-cent12.y + cent23.y)

                -slope23 * cent12.x + slope12 * cent23.x )

               / (slope12 - slope23);

    coord->y = -(coord->x-cent12.x) / slope12 + cent12.y;

    }

  }





/* ----------------------------------------------------------------- */



static void SortFormers(CODEBOOK*   CB,

                        int*        order,

                        Coordinate* coord,

                        IndexSet*   Formers)

{

  double*    angle;

  double     base;

  int        i, j;

  UnitVector uv;

  double     tmpd;

  int        tmpi;



  angle = allocate(IS_size(Formers) * sizeof(double));



  UV_set(&uv,

         VectorScalar(CB, order[IS_index(Formers,0)], 0) - coord->x,

         VectorScalar(CB, order[IS_index(Formers,0)], 1) - coord->y);

  base = angle[0] = UV_angle(&uv);

  angle[0] -= base;



  for( i = 1; i < IS_size(Formers); i++ )

    {

    UV_set(&uv,

           VectorScalar(CB, order[IS_index(Formers,i)], 0) - coord->x,

           VectorScalar(CB, order[IS_index(Formers,i)], 1) - coord->y);

    angle[i] = UV_angle(&uv) - base;

    if( angle[i] < 0.0 )

      {

      angle[i] += 2.0 * M_PI;

      }

    for( j = i; 0 < j && angle[j] < angle[j-1]; j-- )

      {

      tmpd       = angle[j];

      angle[j]   = angle[j-1];

      angle[j-1] = tmpd;

      tmpi = IS_index(Formers, j);

      IS_set(Formers, j, IS_index(Formers, j-1));

      IS_set(Formers, j-1, tmpi);

      }

    }



  deallocate(angle);

}



/* ----------------------------------------------------------------- */



static double VoronoiDistance(Coordinate* coord, VECTORTYPE v)

{

  /* Lasketaan luodun voronoipisteen coord et„isyys CODEBOOKin

     pisteeseen v. Neli”juurta ei k„ytet„. */

  double dx = coord->x - v[0];

  double dy = coord->y - v[1];



  return( dx * dx + dy * dy );

}



/*-------------------------------------------------------------------*/



static int VertexExist(VoronoiGraph* VG, Coordinate* coord)

{

  int i;



  for( i = 0; i < VG_nVertex(VG); i++ )

    {

    if( CO_equal(VG_VertexCoordinate(VG, i), coord) )

      {

      return( i );

      }

    }

  return( ENDOFLINK );

}



/*-------------------------------------------------------------------*/

/*

static int EdgeExist(VoronoiGraph* VG, int vertex, UnitVector* UV)

{

  int i;



  for( i = 0; i < VG_nEdge(VG); i++ )

    {

    if( CO_equal(VG_VertexCoordinate(VG, i), coord) )

      {

      return( YES );

      }

    }

  return( NO );

}

*/

/*-------------------------------------------------------------------*/



static YESNO isVoronoiVertex(VoronoiGraph* VG,

                             CODEBOOK*     CB,

                             int*          order,

                             int           code1,

                             int           code2,

                             int           code3,

                             Coordinate*   coord,

                             IndexSet*     Formers)

{

  /* Lasketaan kolmikon avulla voronoiehdokas coord.

     Tarkastetaan onko jokin kolmikkoon kuulumatton CODEBOOKin

     vektori l„hemp„n„ voronoiehdokasta coord kuin kolmikon vektorit.

     Jos n„in on, niin coord ei ole voronoipiste. */

  double mindist, dist, xdist;

  int    i;



  if( code1 == code2 || code1 == code3 || code2 == code3 )

    {

    return( NO );

    }



  CalculateVertex(CB, order, code1, code2, code3, coord);



  if( VertexExist(VG, coord) != ENDOFLINK )

    {

    return( NO );

    }



  mindist = VoronoiDistance(coord, Vector(CB,order[code1]));

  IS_clear(Formers);

  IS_add(Formers, code1);

  IS_add(Formers, code2);

  IS_add(Formers, code3);



/* printf("iVV %2i %2i %2i %f\n", code1, code2, code3, mindist); */



  /* Koodivektorien j„rjest„minen mahdollistaa oikosulun k„yt”n. */

  for( i = code1 + 1; i < BookSize(CB); i++ )

    {

    if( i != code2 && i != code3 )

      {

      xdist = sqr(VectorScalar(CB, order[i], 0) - coord->x);

      if( xdist < mindist || iszero(xdist - mindist) )

        {

        dist = xdist + sqr(VectorScalar(CB, order[i], 1) - coord->y);

        if( dist < mindist )

          {

          return( NO );

          }

        else

          {

          if( iszero(mindist - dist) )

            {

            IS_add(Formers, i);

            }

          }

        }

      else

        {

        /* Short cut! */

        i = BookSize(CB);

        }

      }

    }



  for( i = code1 - 1; 0 <= i; i-- )

    {

    if( i != code2 && i != code3 )

      {

      xdist = sqr(VectorScalar(CB, order[i], 0) - coord->x);

      if( xdist < mindist || iszero(xdist - mindist) )

        {

        dist = xdist + sqr(VectorScalar(CB, order[i], 1) - coord->y);

        if( dist < mindist )

          {

          return( NO );

          }

        else

          {

          if( iszero(mindist - dist) )

            {

            IS_add(Formers, i);

            }

          }

        }

      else

        {

        /* Short cut! */

        i = -1;

        }

      }

    }



  return( YES );

}



/*-------------------------------------------------------------------*/



static void DirectionVector(UnitVector* UV,

                            VECTORTYPE  v1,

                            VECTORTYPE  v2,

                            VECTORTYPE  v3,

                            Coordinate* coord)

{

  /* Luodaan voronoipisteelle siit„ l„htev„ suunta (yksikk”vektori),

     joka kulkee pisteiden v1 ja v2 v„list„ ja poisp„in v3:sta. */



  UnitVector dir12;

  UnitVector dir13;

  UnitVector dir12left;

  double     angle;



  UV_set(&dir12, v2[0] - v1[0] , v2[1] - v1[1]);

  UV_set(&dir13, v3[0] - v1[0] , v3[1] - v1[1]);

  UV_set(&dir12left, -dir12.dy, dir12.dx);



  angle = UV_dotproduct(&dir12left, &dir13);



  if( 0 < angle ) /* v3 is on the left side of the line v1-v2. */

    {

    UV_set(UV, -dir12left.dx, -dir12left.dy);

    }

  else

    {

    UV_set(UV, dir12left.dx, dir12left.dy);

    }

}



/*-------------------------------------------------------------------*/



static void AddNewVoronoiEdge(VoronoiGraph* VG,

                              CODEBOOK*     CB,

                              int*          order,

                              int           codeA,

                              int           codeB,

                              int           codeC,

                              Coordinate*   coord,

                              int           vertex)

{

  int        e;

  UnitVector dir;



/*   printf("ANVE %2i,%2i,%2i %2i\n", codeA, codeB, codeC, vertex); */



  e = VG_addEdge(VG, vertex, ENDOFLINK);

  VG_setEdgeStatus(VG, e, DIRECTION);

  DirectionVector(&dir,

                  Vector(CB, order[codeA]),

                  Vector(CB, order[codeB]),

                  Vector(CB, order[codeC]),

                  coord);

  VG_setEdgeDir(VG, e, &dir);

  VG_setEdgeCodes(VG, e, order[codeA], order[codeB]);

}



/* ----------------------------------------------------------------- */



static void Linkage(VoronoiGraph* VG, int lastvertex, int prevvor)

{

  /* Luodaan linkki viimeksi lis„tyn pisteen ja sen pisteen (VorInd)

     v„lille, josta p„„stiin t„h„n pisteeseen. */

  int i;

  UnitVector uv;



  UV_set(&uv,

         VG_VertexCoordinate(VG, lastvertex)->x -

         VG_VertexCoordinate(VG, prevvor)->x,

         VG_VertexCoordinate(VG, lastvertex)->y -

         VG_VertexCoordinate(VG, prevvor)->y);



  assert( lastvertex != prevvor );

  for( i = VG_Edge(VG, prevvor);

       i != ENDOFLINK;

       i = VG_EdgeNext(VG, prevvor, i) )

    {

    if( VG_EdgeStatus(VG, i) == DIRECTION )

      {

      if( UV_equal(VG_EdgeDir(VG, i), &uv) )

        {

        VG_setEdgeStatus(VG, i, CONNECTION);

        VG_changeEdgeVertex(VG, i, 1, lastvertex);

        return;

        }

      }

    }



  printf("ERROR:lv=%2i pv=%2i\n", lastvertex, prevvor);

  PrintGraph(VG);

  assert( 0 );

}



/*-------------------------------------------------------------------*/



static int AddVoronoiVertex(VoronoiGraph* VG,

                            CODEBOOK*     CB,

                            int*          order,

                            int           code1,

                            int           code2,

                            int           code3,

                            Coordinate*   coord,

                            IndexSet*     Formers,

                            int           prevvor)

{

  int vertex;

  int i;

  int first, last;



  SortFormers(CB, order, coord, Formers);



  assert( code1 == IS_index(Formers, 0) );

  assert( code2 == IS_index(Formers, 1) ||

          code2 == IS_index(Formers, IS_size(Formers)-1) );





  if( VG_nVertex(VG) > 0 )

    {

    vertex = VG_addVertex(VG, coord);

    Linkage(VG, vertex, prevvor);



    /* Edge between code1 and code2 exists already. */

    if( code2 == IS_index(Formers, 1) )

      {

      first = 1;

      last  = IS_size(Formers)-1;

      }

    else

      {

      first = 0;

      last  = IS_size(Formers)-2;

      }

    }

  else

    {

    vertex = VG_addVertex(VG, coord);

    first = 0;

    last  = IS_size(Formers)-1;

    }



  for( i = first ; i <= last; i++ )

    {

    AddNewVoronoiEdge(VG, CB, order,

                      IS_index(Formers, i),

                      IS_index(Formers, (i+1) % IS_size(Formers)),

                      IS_index(Formers, (i+2) % IS_size(Formers)),

                      coord, vertex);

    }



/*   printf("After AVV\n"); */

/*   PrintGraph(VG); */



  return( vertex );

}





/* ----------------------------------------------------------------- */



static void GenerateEdge(VoronoiGraph* VG,

                         CODEBOOK*     CB,

                         int*          order,

                         int           code1,

                         int           code2,

                         int           minimi,

                         int           maksimi)

{

  /* Luodaan kahden erillisen pisteen v„lille niit„ erottava viiva.*/

  VECTORTYPE v1 = Vector(CB, order[code1]);

  VECTORTYPE v2 = Vector(CB, order[code2]);

  Coordinate coord1, coord2, cent;

  double     dir[2];

  double     tmpSX1, tmpSY1, tmpSX2, tmpSY2;

  UnitVector direction;

  int        vertex1, vertex2, edge;



  cent.x = (v1[0]+v2[0]) / 2.0;

  cent.y = (v1[1]+v2[1]) / 2.0;

  dir[0] = v1[1] - v2[1];

  dir[1] = v1[0] - v2[0];



  /* Suuntavektorin ensimm„isen komponentin tulee osoittaa oikealle */

  if( dir[0] <= 0 )

    {

    dir[0] = -dir[0];

    }

  else

    {

    dir[1] = -dir[1];

    }



  if( iszero(dir[0]) )

    {

    coord1.x = coord2.x = cent.x;

    coord1.y = minimi;

    coord2.y = maksimi;

		}

  else

    {

    if( iszero(dir[1]) )

      {

      coord1.y =  coord2.y = cent.y;

      coord1.x = minimi;

      coord2.x = maksimi;

      }

    else

      {

      tmpSX1 = (maksimi - cent.x) / dir[0];

      tmpSX2 = (cent.x - minimi) / dir[0];

      if( dir[1] > 0 )

        {

        tmpSY1 = (maksimi - cent.y) / dir[1];

        tmpSY2 = (cent.y - minimi) / dir[1];

        }

      else

        {

        tmpSY2 = (maksimi - cent.y) / -dir[1];

        tmpSY1 = (cent.y - minimi) / -dir[1];

        }

      if( tmpSX1 < tmpSY1 )

        {

        coord1.x = maksimi;

        coord1.y = cent.y + tmpSX1 * dir[1];

        }

      else

        {

        coord1.x = cent.x + tmpSY1 * dir[0];

        coord1.y = ( dir[1] < 0 ? minimi :  maksimi );

        }

      if( tmpSX2 < tmpSY2 )

        {

        coord2.x = minimi;

        coord2.y = cent.y - tmpSX2 * dir[1];

        }

      else

        {

        coord2.x = cent.x - tmpSY2 * dir[0];

        coord2.y =( dir[1] < 0 ? maksimi : minimi );

        }

      }

    }



  vertex1 = VG_addVertex(VG, &coord1);

  vertex2 = VG_addVertex(VG, &coord2);

  edge = VG_addEdge(VG, vertex1, vertex2);



  VG_setEdgeStatus(VG, edge, CONNECTION);

  UV_set(&direction, dir[0], dir[1]);

  VG_setEdgeDir(VG, edge, &direction);

  VG_setEdgeCodes(VG, edge, order[code1], order[code2]);

}



/* ----------------------------------------------------------------- */



static void AllCodevectorsInLine(VoronoiGraph* VG, CODEBOOK* CB, int* order)

  {

  /* Jos kaikki pisteet ovat linjassa, niin muodostetaan eri„vien

     vierekk„isten pisteiden keskipisteen kautta kulkeva voronoiviiva.

     Huom. Pisteet ovat kasvavassa j„rjestyksess„. */

  int i;



  for(i = 1; i < BookSize(CB); i++)

    {

    GenerateEdge(VG, CB, order, i-1, i, CB->MinValue, CB->MaxValue);

    }

}



/* ----------------------------------------------------------------- */



static void NextState(IndexSet* Formers, int* code1, int* code2, int* code3)

{

  if( *code2 == IS_index(Formers, 1) )

    {

    *code2 = IS_index(Formers, IS_size(Formers)-1);

    }

  else

    {

    *code2 = IS_index(Formers, 1);

    }

  *code3 = *code1;

}



/* ----------------------------------------------------------------- */



static void GenerateGraph(VoronoiGraph* VG, CODEBOOK* CB, int* order)

{

  int        prevvor, lastvertex;

  int        code1 = 0;

  int        code2 = 1;

  int        code3 = 1;

  Coordinate coord;

  YESNO      tripleok;

  YESNO      vertexok;

  STACK*     S;

  IndexSet*  Formers;



  assert( BookSize(CB) > 2 );



  S = S_make();

  Formers = IS_make(6); /* Basic need is 3, so there is little extra. */



  do

    {

    do

      {

      code3++;

      } while( code3 < BookSize(CB) &&

               OnSameLine(CB, order, code1, code2, code3) );



    vertexok = code3 < BookSize(CB)

               ? isVoronoiVertex(VG, CB, order,

                                 code1, code2, code3,

                                 &coord, Formers)

               : NO;

    } while( code3 < BookSize(CB) && !vertexok );



  if( vertexok )

    {

    prevvor = AddVoronoiVertex(VG, CB, order,

                               code1, code2, code3,

                               &coord, Formers, 0);

    if( code2 == IS_index(Formers, 1) )

      {

      TripleToStack(S, IS_index(Formers, 0), IS_index(Formers, 1), prevvor, 0);

      StatesToStack(S, Formers, prevvor);

      }

    else

      {

      StatesToStack(S, Formers, prevvor);

      TripleToStack(S, IS_index(Formers, 0), IS_index(Formers, IS_size(Formers)-1), prevvor, 0);

      }

    NextState(Formers, &code1, &code2, &code3);



    tripleok = YES;

    while( tripleok )

      {

      do

        {

        do

          {

          code3++;

          } while( code3 < BookSize(CB) &&

                   OnSameLine(CB, order, code1, code2, code3) );

        vertexok = code3 < BookSize(CB)

                   ? isVoronoiVertex(VG, CB, order,

                                     code1, code2, code3,

                                     &coord, Formers)

                   : NO;

        } while( code3 < BookSize(CB) && !vertexok );



      if( vertexok )

        {

        lastvertex = AddVoronoiVertex(VG, CB, order,

                                      code1, code2, code3,

                                      &coord, Formers, prevvor);

        StatesToStack(S, Formers, lastvertex);

        prevvor = lastvertex;



        NextState(Formers, &code1, &code2, &code3);

        tripleok = YES;

        }

      else

        {

        if( !S_empty(S) )

          {

          TripleFromStack(S, &code1, &code2, &prevvor, &code3);

          code3 = code1;

          tripleok = YES;

          }

        else

          {

          tripleok = NO;

          }

        }

      }

    }

  else

    {

    AllCodevectorsInLine(VG, CB, order);

    }



  S_free(S);

  IS_free(Formers);

}



/* ================================================================= */



static void CheckGraph(VoronoiGraph* VG)

{

  int i, j;



  for( i = 0; i < VG_nEdge(VG) - 1; i++ )

     {

     if( VG_EdgeStatus(VG, i) == DIRECTION )

       {

        for( j = i + 1; j < VG_nEdge(VG); j++ )

           {

           if( VG_EdgeStatus(VG, j) == DIRECTION )

             {

/*               if( UV_opposite(VG_EdgeDir(VG, i), VG_EdgeDir(VG, j)) ) */

              if( (VG_EdgeCode(VG, i, 0) == VG_EdgeCode(VG, j, 0) &&

                   VG_EdgeCode(VG, i, 1) == VG_EdgeCode(VG, j, 1)) ||

                  (VG_EdgeCode(VG, i, 0) == VG_EdgeCode(VG, j, 1) &&

                   VG_EdgeCode(VG, i, 1) == VG_EdgeCode(VG, j, 0)) )

                {

                VG_setEdgeStatus(VG, i, CONNECTION);

                VG_changeEdgeVertex(VG, i, 1, VG_Vertex(VG, j, 0));

                VG_removeEdge(VG, j);

                }

             }

           }

       }

     }

}



/*-------------------------------------------------------------------*/



static YESNO isInside(Coordinate* coord, int minvalue, int maxvalue)

{

  return( minvalue <= coord->x && coord->x <= maxvalue &&

          minvalue <= coord->y && coord->y <= maxvalue );

}



/*-------------------------------------------------------------------*/



static YESNO isStrictlyInside(Coordinate* coord, int minvalue, int maxvalue)

{

  return( minvalue < coord->x && coord->x < maxvalue &&

          minvalue < coord->y && coord->y < maxvalue );

}



/*-------------------------------------------------------------------*/



static YESNO isOnBorder(Coordinate* coord, int minvalue, int maxvalue)

{

  /* Equality for 'double' is uncertain. */

  return( ( (minvalue == coord->x || coord->x == maxvalue) &&

            (minvalue <= coord->y && coord->y <= maxvalue) ) ||

          ( (minvalue <= coord->x && coord->x <= maxvalue) &&

            (minvalue == coord->y || coord->y == maxvalue) ) );

}



/*-------------------------------------------------------------------*/



static void CheckEdges(VoronoiGraph* VG, int vertex)

{

  int i;



  for( i = VG_Edge(VG, vertex);

       i != ENDOFLINK;

       i = VG_EdgeNext(VG, vertex, i) )

    {



    }

}



/*-------------------------------------------------------------------*/



static void SearchBorderVertex(VoronoiGraph* VG,

                               int           edge,

                               int           minvalue,

                               int           maxvalue)

{

  /* Jos piste sijaitsee alueen sis„ll„, niin etsit„„n alueen rajalta

     piste, johon linkki asetetaan osoittamaan */

  Edge*      E = &VG->edge[edge];

  Vertex*    V = &VG->vertex[VG_Vertex(VG, edge, 0)];

  Coordinate coord;

  double     apuSX, apuSY;

  int        vertex;



  if( iszero(E_dirx(E)) )

    {

    coord.x = VX_x(V);

    coord.y = ( E_diry(E) > 0 ? maxvalue : minvalue);

		}

  else if( iszero(E_diry(E)) )

    {

    coord.y = VX_y(V);

    coord.x = ( E_dirx(E) > 0 ? maxvalue : minvalue );

		}

  else

    {

    apuSX = ( E_dirx(E) > 0

              ? (maxvalue - VX_x(V)) / E_dirx(E)

              : (VX_x(V) - minvalue) / -E_dirx(E) );



    apuSY = ( E_diry(E) > 0

              ? (maxvalue - VX_y(V)) / E_diry(E)

              : (VX_y(V) - minvalue) / -E_diry(E) );



    if( apuSX < apuSY )

      {

      coord.x = ( E_dirx(E) > 0 ? maxvalue : minvalue );

      coord.y = VX_y(V) + apuSX * E_diry(E);

      }

    else

      {

      coord.x = VX_x(V) + apuSY * E_dirx(E);

      coord.y = ( E_diry(E) < 0 ? minvalue : maxvalue );

      }

    }



  vertex = VertexExist(VG, &coord);

  if( vertex == ENDOFLINK )

    {

    vertex = VG_addVertex(VG, &coord);

    }

  VG_setEdgeStatus(VG, edge, CONNECTION);

  VG_changeEdgeVertex(VG, edge, 1, vertex);

}



/* ----------------------------------------------------------------- */



static YESNO CalculateIntersectionPoint(Coordinate* base1,

                                        UnitVector* dir1,

                                        Coordinate* base2,

                                        UnitVector* dir2,

                                        Coordinate* intersection)

{

  double coeff1, coeff2;

  double const1, const2;



  if( iszero(dir1->dx) )

    {

    if( iszero(dir2->dx) )

      {

      return( NO );

      }

    else

      {

      coeff2 = dir2->dy / dir2->dx;

      const2 = base2->y - coeff2 * base2->x;

      intersection->x = base1->x;

      intersection->y = coeff2 * base1->x + const2;

      }

    }

  else

    {

    if( iszero(dir2->dx) )

      {

      coeff1 = dir1->dy / dir1->dx;

      const1 = base1->y - coeff1 * base1->x;

      intersection->x = base2->x;

      intersection->y = coeff1 * base2->x + const1;

      }

    else

      {

      coeff1 = dir1->dy / dir1->dx;

      coeff2 = dir2->dy / dir2->dx;

      const1 = base1->y - coeff1 * base1->x;

      const2 = base2->y - coeff2 * base2->x;



      assert( !iszero(coeff1 - coeff2) );



      intersection->x = (const2 - const1) / (coeff1 - coeff2);

      intersection->y = (coeff1 * const2 - coeff2 * const1) / (coeff1 - coeff2);

      }

    }



  return( YES );

}



/* ----------------------------------------------------------------- */



static void CutEdge(VoronoiGraph* VG,

                    int           edge,

                    double        startx,

                    double        starty,

                    double        endx,

                    double        endy,

                    YESNO*        removed)

{

  Coordinate cutbase;

  UnitVector cutdir;

  UnitVector cutdirleft;

  UnitVector dir1, dir2;

  int        vertex1, vertex2;

  double     angle1, angle2;

  Coordinate intersection;

  int        newvertex;



  UV_set(&cutdir, endx - startx, endy - starty);

  UV_set(&cutdirleft, -cutdir.dy, cutdir.dx);

  vertex1 = VG_Vertex(VG, edge, 0);

  vertex2 = VG_Vertex(VG, edge, 1);

  UV_set(&dir1,

         VG_VertexCoordinate(VG, vertex1)->x - startx,

         VG_VertexCoordinate(VG, vertex1)->y - starty);

  UV_set(&dir2,

         VG_VertexCoordinate(VG, vertex2)->x - startx,

         VG_VertexCoordinate(VG, vertex2)->y - starty);

  angle1 = UV_dotproduct(&cutdirleft, &dir1);

  angle2 = UV_dotproduct(&cutdirleft, &dir2);



  if( 0 <= angle1 && 0 <= angle2 ) /* Both on the left side */

    {

    return;

    }

  else

    {

    if( angle1 <= 0 && angle2 <= 0 ) /* Both on the right side */

      {

      *removed = YES;

      return;

      }

    else

      {

      if( angle2 < angle1 )

        {

        VG_turnEdgeAround(VG, edge);

        vertex1 = vertex2;

        }

      cutbase.x = startx;

      cutbase.y = starty;

      CalculateIntersectionPoint(VG_VertexCoordinate(VG, vertex1),

                                 VG_EdgeDir(VG, edge),

                                 &cutbase,

                                 &cutdir,

                                 &intersection);

      newvertex = VertexExist(VG, &intersection);

      if( newvertex == ENDOFLINK )

        {

        newvertex = VG_addVertex(VG, &intersection);

        }

      VG_changeEdgeVertex(VG, edge, 0, newvertex);

      }

    }

}



/* ----------------------------------------------------------------- */



static void LimitEdge(VoronoiGraph* VG,

                      int           edge,

                      int           minvalue,

                      int           maxvalue,

                      YESNO*        removed)

{

  CutEdge(VG, edge, 0, 0, maxvalue, 0, removed);

  if( !*removed )

    {

    CutEdge(VG, edge, maxvalue, 0, maxvalue, maxvalue, removed);

    if( !*removed )

      {

      CutEdge(VG, edge, maxvalue, maxvalue, 0, maxvalue, removed);

        {

        if( !*removed )

          {

          CutEdge(VG, edge, 0, maxvalue, 0, 0, removed);

          }

        }

      }

    }

}



/* ----------------------------------------------------------------- */



static void MoveVertexToLine(VoronoiGraph* VG,

                             int           edge,

                             double        startx,

                             double        starty,

                             double        endx,

                             double        endy,

                             YESNO*        removed)

{

  Coordinate cutbase;

  UnitVector cutdir;

  UnitVector cutdirleft;

  UnitVector vertexdir;

  int        vertex;

  double     dirangle, vertexangle;

  Coordinate intersection;

  int        newvertex;



  UV_set(&cutdir, endx - startx, endy - starty);

  UV_set(&cutdirleft, -cutdir.dy, cutdir.dx);



  vertex = VG_Vertex(VG, edge, 0);

  UV_set(&vertexdir,

         VG_VertexCoordinate(VG, vertex)->x - startx,

         VG_VertexCoordinate(VG, vertex)->y - starty);

  vertexangle = UV_dotproduct(&cutdirleft, &vertexdir);



  if( vertexangle <= 0 )  /* Vertex is on the right side of the line. */

    {

    dirangle = UV_dotproduct(&cutdirleft, VG_EdgeDir(VG, edge));



    if( dirangle <= 0 )   /* Direction goes to the right. */

      {

      *removed = YES;

      }

    else

      {

      cutbase.x = startx;

      cutbase.y = starty;

      CalculateIntersectionPoint(VG_VertexCoordinate(VG, vertex),

                                 VG_EdgeDir(VG, edge),

                                 &cutbase,

                                 &cutdir,

                                 &intersection);

      newvertex = VertexExist(VG, &intersection);

      if( newvertex == ENDOFLINK )

        {

        newvertex = VG_addVertex(VG, &intersection);

        }

      VG_changeEdgeVertex(VG, edge, 0, newvertex);

      }

    }

}



/* ----------------------------------------------------------------- */



static void MoveVertexToBorder(VoronoiGraph* VG,

                               int           edge,

                               int           minvalue,

                               int           maxvalue,

                               YESNO*        removed)

{

  MoveVertexToLine(VG, edge, 0, 0, maxvalue, 0, removed);

  if( !*removed )

    {

    MoveVertexToLine(VG, edge, maxvalue, 0, maxvalue, maxvalue, removed);

    if( !*removed )

      {

      MoveVertexToLine(VG, edge, maxvalue, maxvalue, 0, maxvalue, removed);

        {

        if( !*removed )

          {

          MoveVertexToLine(VG, edge, 0, maxvalue, 0, 0, removed);

          }

        }

      }

    }

}



/* ----------------------------------------------------------------- */



static void RemoveExtraVertices(VoronoiGraph* VG)

{

  int i;



  i = 0;

  while( i < VG_nVertex(VG) )

    {

    if( VG_Edge(VG, i) == ENDOFLINK )

      {

      VG_removeVertex(VG, i);

      }

    else

      {

      i++;

      }

    }

}



/*-------------------------------------------------------------------*/



static void ModifyGraph(VoronoiGraph* VG, int minvalue, int maxvalue)

{

  /* Poistetaan ne voronoipisteet, jotka eiv„t ole vektoriavaruudessa.

     Vektoriavaruuden pisteiden DIRECTION kaaret asetetaan osoittamaan

     avaruuden reunalla olevaan keinotekoiseen pisteeseen. */



  int   i;

  int   vertex1, vertex2;

  YESNO removed;

  YESNO inside1, inside2;



  i = 0;

  while( i < VG_nEdge(VG) )

    {

    removed = NO;

    switch( VG_EdgeStatus(VG, i) )

      {

      case DIRECTION:

        {

        vertex1 = VG_Vertex(VG, i, 0);

        if( isStrictlyInside(VG_VertexCoordinate(VG, vertex1), minvalue, maxvalue) )

          {

          SearchBorderVertex(VG, i, minvalue, maxvalue);

          }

        else

          {

          MoveVertexToBorder(VG, i, minvalue, maxvalue, &removed);

          if( !removed )

            {

            SearchBorderVertex(VG, i, minvalue, maxvalue);

            }

          }

        break;

        }

      case CONNECTION:

        {

        vertex1 = VG_Vertex(VG, i, 0);

        vertex2 = VG_Vertex(VG, i, 1);

        inside1 = isInside(VG_VertexCoordinate(VG, vertex1), minvalue, maxvalue);

        inside2 = isInside(VG_VertexCoordinate(VG, vertex2), minvalue, maxvalue);

        if( !inside1 || !inside2 )

          {

          if( isStrictlyInside(VG_VertexCoordinate(VG, vertex1),

                               minvalue,

                               maxvalue) )

            {

            /* inside2 == NO */

            VG_changeEdgeVertex(VG, i, 1, ENDOFLINK);

            SearchBorderVertex(VG, i, minvalue, maxvalue);

            }

          else

            {

            if( isStrictlyInside(VG_VertexCoordinate(VG, vertex2),

                                 minvalue,

                                 maxvalue) )

              {

              /* inside1 == NO */

              VG_turnEdgeAround(VG, i);

              VG_changeEdgeVertex(VG, i, 1, ENDOFLINK);

              SearchBorderVertex(VG, i, minvalue, maxvalue);

              }

            else

              {

              LimitEdge(VG, i, minvalue, maxvalue, &removed);

              }

            }

          }

        break;

        }

      case BORDER:

      default:

        {

        break;

        }

      }

    if( removed )

      {

      VG_removeEdge(VG, i);

      }

    else

      {

      i++;

      }

    }



  RemoveExtraVertices(VG);

}



/* ================================================================= */



void VG_generate(VoronoiGraph* VG, CODEBOOK* CB)

{

  int* order;



  if( BookSize(CB) < 2 )

    {

    /* No edges. */

    return;

    }

  else

    {

    if( BookSize(CB) == 2 )

      {

      /* Only one edge. */

      GenerateEdge(VG, CB, order, 0, 1, CB->MinValue, CB->MaxValue);

      }

    else

      {

      GenerateOrder(CB, &order);



      GenerateGraph(VG, CB, order);



/*       printf("After GG\n"); */

/*       PrintGraph(VG); */



      CheckGraph(VG);



/*       printf("After CG\n"); */

/*       PrintGraph(VG); */



      ModifyGraph(VG, CB->MinValue, CB->MaxValue);



/*       printf("After MG\n"); */

/*       PrintGraph(VG); */



      FreeOrder(order);

      }

    }

}



