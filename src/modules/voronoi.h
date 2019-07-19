#ifndef __VORO_H
#define __VORO_H


#include "cb.h"


typedef enum { DIRECTION, CONNECTION, BORDER } EdgeStatus;
#define ENDOFLINK  -1


/* ----------------------------------------------------------------- */

typedef struct
          {
          double x, y;
          } Coordinate;

/* ----------------------------------------------------------------- */

typedef struct
          {
          double dx, dy;
          } UnitVector;

/* ----------------------------------------------------------------- */

typedef struct
          {
          Coordinate coord;
          int        edge;
          } Vertex;

/* ----------------------------------------------------------------- */

typedef struct
          {
          EdgeStatus status;
          int        vertex[2];
          int        next[2];
          UnitVector direction;
          int        codevector[2];
          } Edge;

/* ----------------------------------------------------------------- */

typedef struct
          {
          Vertex* vertex;       /* Array of vertices */
          Edge*   edge;         /* Array of edges */
          int     nvertex;      /* Number of vertices */
          int     nedge;        /* Number of edges */
          int     nallocvertex; /* Number of allocated space for vertices */
          int     nallocedge;   /* Number of allocated space for edges */
          } VoronoiGraph;

/* ----------------------------------------------------------------- */

#define VX_x(VX)    ((VX)->coord.x)
#define VX_y(VX)    ((VX)->coord.y)
#define VX_edge(VX) ((VX)->edge)

#define E_status(E)      ((E)->status)
#define E_vertex(E,end)  ((E)->vertex[end])
#define E_next(E,end)    ((E)->next[end])
#define E_dirx(E)        ((E)->direction.dx)
#define E_diry(E)        ((E)->direction.dy)
#define E_code(E,end)    ((E)->codevector[end])

/* ----------------------------------------------------------------- */

VoronoiGraph* VG_make(CODEBOOK* CB);
void          VG_kill(VoronoiGraph* VG);
void          VG_generate(VoronoiGraph* VG, CODEBOOK* CB);

int           VG_addVertex(VoronoiGraph* VG, Coordinate* c);
int           VG_addEdge(VoronoiGraph* VG, int vertex1, int vertex2);
void          VG_removeVertex(VoronoiGraph* VG, int vertex);
void          VG_removeEdge(VoronoiGraph* VG, int edge);


void          VG_setEdgeStatus(VoronoiGraph* VG, int edge, EdgeStatus status);
void          VG_setEdgeDir(VoronoiGraph* VG, int edge, UnitVector* dir);
void          VG_setEdgeCodes(VoronoiGraph* VG, int edge, int code1, int code2);
void          VG_setEdgeVertex(VoronoiGraph* VG, int edge, int end, int vertex);
int           VG_Vertex(VoronoiGraph* VG, int edge, int end);
EdgeStatus    VG_EdgeStatus(VoronoiGraph* VG, int edge);
int           VG_EdgeNext(VoronoiGraph* VG, int vertex, int edge);
UnitVector*   VG_EdgeDir(VoronoiGraph* VG, int edge);
int           VG_EdgeCode(VoronoiGraph* VG, int edge, int end);

/* ----------------------------------------------------------------- */

#define VG_nVertex(VG)  ((VG)->nvertex)
#define VG_nEdge(VG)    ((VG)->nedge)

#define VG_Edge(VG, vx)             ((VG)->vertex[vx].edge)
#define VG_VertexCoordinate(VG, vx) (&(VG)->vertex[vx].coord)

/* ----------------------------------------------------------------- */


#if 0

/**************************   Voronoi piste   *****************************/

typedef double   VorPiste[2];
typedef struct {
               VorPiste LinkkiKohde;
               char     LinkkiEhto;     /* suunta vai kohdepiste */
               int      CBnIndeksit[2]; /* kertoo mink„ kahden pisteen
                                           v„list„ linkki kulkee. */
               } Linkki;
typedef struct {
               VorPiste Piste;
               Linkki   Linkit[3];
               } VoronoiPiste;

typedef VoronoiPiste*  VoronoiPisteet;
typedef struct {
               VoronoiPisteet VPisteet;
               int            IdTunnus;
               int            Koko;      /* montako voronoipistett„ */
               int            Dimensio;  /* pisteen ulottuvuus */
               } VoronoiVektori;

/**************************************************************************/

typedef struct {
               int x;
               int y;
               } CBnIndeksit;

/**************************    SoluVektori    ****************************/

struct Alkio {
             int Vind;
             struct Alkio* seuraava;
             };

typedef struct {
               struct Alkio* Alkiot;
               } SoluVektori;

/*************************************************************************/

void AlustaVoronoiVektori(VoronoiVektori* VV, CODEBOOK* CB);
void VapautaMuistiaVoronoiVektorilta(VoronoiVektori* VV);
void AlustaSoluVektori(CODEBOOK*CB, SoluVektori* SV);
void VapautaSoluVektori(SoluVektori* SV, CODEBOOK* CB);

/* T„ll„ funktiolla luodaan voronoidiagrammi */
void VoronoiDiagrammi(CODEBOOK*       CB,
                      VoronoiVektori* VV,
                      CBnIndeksit     CBi);

/* T„ll„ funktiolla luodaan soluvektori */
void LuoSoluVektori(CODEBOOK*       CB,
                    VoronoiVektori* VV,
                    SoluVektori*    SV,
                    CBnIndeksit     CBi);

#endif

#endif /* __VORO_H */

