#if ! defined(__EDM_H)
#define __EDM_H

#include "owntypes.h"
#include "contree.h"

typedef BYTE* PBYTE;

typedef struct { int  ImageX;
                 int  ImageY;
                 size_t ImageTotal;
                 int  ClustersX;
                 int  ClustersY;
                 int  ClustersT;
                 int  BlocksX;
                 int  BlocksY;
                 int  Action;
                 int  ModellingType;
                 int  BlockModel;
                 int  BlockContext;
                 int  BlockSize;
                 int  ClusterSize;
                 int  Clustering;
                 int  ContextSize;
                 int  ContextTreeMode;
                 int  LayerContextSize;
                 size_t  NumberOfStates;
                 size_t  FirstBlockContext;

                 int  Position;       /* indicates the current position */
                                      /* in the EDM-file (EDMPOSITION)  */

                 int  Preview;
                 int  OverWrite;
                 int  QuietLevel;
                 int  FastMode;
                 int  ForceJBIG;
                 int  Stage;

                 /* technical purpose data */

                 int  ContextTreeDepth0;
                 long ContextsInTree;
                 CONTREE  ContextTree;  /* Context Tree */

                 size_t   IndexesPos;   /* Cluster Indexes begining position */
                 size_t   DataPos;      /* Pixel Data begining position      */
                 size_t   EndOfData;    /* End of compressed data            */

                 BYTE*   StateIndex;
                 size_t*   ClusterIndex;
                 unsigned short*  ClusterDelta;
                 PBYTE*  BlockCode;

                 FILE*  EDMFile;    /* used for EDMView */
                 FILE*  ModelFile;
                 FILE*  ContreeFile;
                 IMAGE* Image;
                 IMAGE* LayerImage;

                 } EDMPARAMETERS;

typedef struct {
                 long*   ContextUsed;
                 long*   ContextWhites;
                 double* ContextEntropy;
                 double  DynamicEntropy;

                 int     Cw_;
                 long    Bt_;
                 long    Bw_;
                 long    Bb_;
                 long    Bm_;

                 long    BytesTotal;
                 long    PixelDataTotal;
                 long    BlockDataTotal;

                 } EDMSTATISTIC;

extern EDMSTATISTIC EDMStatistic;

/* ----------------------- Global symbols ---------------------------*/

#define  FormatNameEDM    "edm"
#define  FormatNameQMT    "qmt"
#define  DEFAULTJBIGCONTEXT 10
#define  DEFAULTEDMCONTEXT  10


typedef  enum { ALLWHITE=0, ALLBLACK=1, MIXED=2 } BLOCKTYPE;
typedef  enum { NOACTION=0, COMPRESSION=1, DECOMPRESSION=2, PREVIEW=3 } ACTIONTYPE;
typedef  enum { LEFT=0, ABOVE=1 } DIRECTIONTYPE;
typedef  enum { NOWHERE=0, EDMBEGIN=1, EDMTEXTHEADER=2, EDMPREVIEW=3,
                EDMCLUSTERINDEX=4, EDMCONTEXTTREE=5, EDMMODEL=6, EDMPIXELDATA=7, EDMEND=8 }
                EDMPOSITION;
/* it means: position in the EDM-file after xxxx, status EDMPIXELDATA */
/* means position inside pixelwise data (between clusters)            */

/*-----------------------  ERROR CODES  -------------------------*/

typedef enum {
  EDM_OTHER=0, EDM_NOMEMORY, UMETHOD, UMODEL, UBLOCKTYPE,
  CLUSTERDELTA, WRMODEL, RMODEL, QMEXPECTED, NOMODEL, QMTCONTEXT,
  WRONGVERSION, STOREBLOCKT_1, STOREBLOCKT_2, SAMENAME, NESIZES,
  UNSUPMODEL, UNSUPLAYER, NOPREVIEW, NOCLUSTERS, WRONGSIZE, NOCONTREE,
  UNSUPCT, WRCT, RCT, CTCONTEXTS, CTWRONGCT
  } EDM_ERROR;

EDM_ERROR     EDMErrorStatus();

/*----------------------- INFORNATIVE FUNCTIONS -----------------*/

#define PreviewIsPossible(PEDMP) ( (PEDMP)->BlockModel != 0 )
/* PIXELWISE is 0 */

#define ClusteringIsPossible(PEDMP) ( (PEDMP)->Clustering > 1)

/* ======================== Interface *batch mode* ===================*/

YESNO EncodeImage(IMAGE* Image, FILE* OutputFile, IMAGE* LayerImage,
                 IMAGE* PreImage, FILE* ModelFile, FILE* ContreeFile);

YESNO DecodeImage(IMAGE* Image, FILE* InputFile, IMAGE* LayerImage,
                 IMAGE* PreImage, FILE* ModelFile, FILE* ContreeFile);

YESNO PreviewImage(IMAGE* Image, FILE* InputFile);

/* function return pointers to IMAGE anf FILE.                              */
/* you must pass to function file names and pointer to filled EDMPARAMETERS */

YESNO InitializeEDM  (EDMPARAMETERS* EDMpointer,
                     IMAGE* Image, IMAGE* LayerImage, IMAGE* PreImage,
                     char* InputName, char* OutputName,  char* LayerName,
                     char* PreviewName, char* ModelName, char* ContreeName,
                     FILE** InputFile,  FILE** OutputFile, FILE** ModelFile,
                     FILE** ContreeFile);

YESNO EDMDone(IMAGE* Image, IMAGE* LayerImage, IMAGE* PreImage,
             FILE* InputFile,  FILE* OutputFile, FILE* ModelFile,
             FILE* ContreeFile);

/* ===================== Interface for EDM View ======================= */
/*                                                                      */
/* EDM Image Structure:                                                 */
/*                                                                      */
/*    1) Text Header                                                    */
/*    2) Header                                                         */
/*       - Preview Data (QM-stream)                                     */
/*       - Image Data (Model, Cluster indecies)                         */
/*    3) Clusters (cluster by cluster, pixel by pixel, QM-coded)        */
/*                                                                      */
/* ===================== Interface for EDM View ======================= */

YESNO OpenEDMImage(EDMPARAMETERS* EDMpointer, char* EDMName,
                   /* char* LayerName, char* ModelName, char* ContreeName */
                   int QuietLevel);
/* Open EDM image, read EDM Text Header */

YESNO CloseEDMImage(EDMPARAMETERS* EDMpointer);

YESNO ReadEDMPreviewData(EDMPARAMETERS* EDMpointer);
/* Read Preview Data to the structure.                                 */
/* Now, one can use this preview data directly from the EDM-structure, */
/* or use WritePreview function directly after call this function to   */
/* store Preview Image into the PBM image (IMAGE-structure) into the   */
/* file or memory.                                                     */

YESNO WriteEDMPreviewImage (EDMPARAMETERS* EDMpointer, IMAGE* Image);
/* Store Preview Image into the Image.                                 */

YESNO ReadEDMImageData(EDMPARAMETERS* EDMpointer);
/* Read the rest of Header, or read all Header, if it isn't yet read.  */

YESNO ReadEDMCluster(EDMPARAMETERS* EDMpointer,
                    int ClusterX, int ClusterY,       /* human notion */
                    IMAGE* Image, int PosX, int PosY);
/* Read cluster and decompress it to the give IMAGE-structure.          */
/* DIRECT ACCESS to the EDM file must be alloved here.                  */

#endif /* __EDM_H */
