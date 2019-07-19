/*-------------------------------------------------------------------*/
/* EDM.C           Eugene Ageenko                                    */
/*                 Pasi Fr„nti                                       */
/*                                                                   */
/*  Binary images compressor for EDM (main library code)             */
/*                                                                   */
/* - Three-stage modelling scheme for binary image compression.      */
/* - Modelling: (1/2) cluster/block model (3) pixelwise model (JBIG) */
/* - Cluster modelling: fixed size (NxN)                             */
/* - Block modelling: none (pixelwise), all-white/non-white or       */
/*   all-white/all-block/mixed fixed size (MxM) blocks               */
/* - Block contexts: none, previous, or two neighboring blocks.      */
/* - Context modelling: local templates + external layers.           */
/* - Modelling type: FORWARD-adaptive, or adaptive (dynamic)         */
/*   (on first case the model is semi-adaptive or static)            */
/* - Fast mode supported (floating point operations minimized)       */
/* - Coding by QM-coder (including modified w/o ESC-codes one)       */
/* - JBIG mode supported (w/o clusters, pixelwise dynamic modelling) */
/* - CONTEXT TREE supported                                          */
/*                                                                   */
/*  SUPPORT ONE IMAGE COMPRESSED AT A TIME ONLY !!!!!!!              */
/*                                                                   */
/*  USE IT WITH QuietLevel <= 2 to prevent debug info appearing      */
/*-------------------------------------------------------------------*/

#define  ProgName       "EDM"
#define  VersionNumber  "Version 0.07"
#define  LastUpdated    "9.11.97"

#define IDString        "EDM-CT-2"

/* ----------------------------------------------------------------- */

#include <math.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <assert.h>

#include "file.h"
#include "binimage.h"
#include "memctrl.h"
#include "pgm.h"
#include "qm.h"
#include "error.h"
#include "owntypes.h"
#include "edm.h"
#include "contree.h"

extern void DecodeBlockType(FILE* f, int blockindex, int blockline, int* blocktype);
extern void EncodeBlockType(FILE* f, int blockindex, int blockline, int blocktype);

#define  MAXBLOCKS       65532
#define  MAXTEMPLATE        22
#define  MAXCLUSTERS     16383
#define  MAXCLUSTERSIZE    724
#define  MAXCLUSTERLENGTH (MAXCLUSTERSIZE*((long)MAXCLUSTERSIZE))/8
#define  CONTEXTHEIGHT      10
#define  DEFCLUSTERSIZE    256
#define  DEFBUFSIZE        128

/* typedef enum {CHESSTYPE=1, NOCOMPRESS=0xFFFE, FCLUSTERCODE=0xFFFF} SPECIALCLUSTERCODES;
*/

#define  FCLUSTERCODE   0xFFFF
#define  CHESSTYPE           1
#define  NOCOMPRESS     0xFFFE
#define  NUL            0L

/* !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! */
/* constants defined in your *.fac file must be sinchronized with these */
typedef enum {PIXELWISE=0, BLOCK=1, EXTENDEDBLOCK=2} BLOCKMODEL;
typedef enum {MODEL_0=0, MODEL_1=1, MODEL_2=2} BLOCKCONTEXT;
typedef enum {STATIC=0, SEMIADAPTIVE=1, ADAPTIVE=2} QMMODELLINGTYPE;
typedef enum {NO_CT=0, STATIC_CT=1, SEMIADAPTIVE_CT=2, ADAPTIVE_CT=3} CONTREETYPE;

/* !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! */
/* this pointer must be assigned in order to use functions from EDM.c  */
EDMPARAMETERS* PEDMParameters;

EDMSTATISTIC EDMStatistic;

#define Para_(x)   (PEDMParameters->x)
#define Result(x)  (EDMStatistic.x)

/* When you call function from EDM.c you must first assign EDMParameters, */
/* then pass pointer to this structure to PEDMParameters.                 */


/* ==================== ERROR HANDLING ==============================*/

/* In the case of serious internal error special error function is   */
/* called from this module to hadle the error situations, default it */
/* print the message (through the ErrorMsg), breaks the program and  */
/* exit to OS, otherwise (in case of non-serious) error,             */
/* function returns NO and error status can be obtained from the     */
/* EDMErrorStatus().                                                 */

EDM_ERROR    EDMErrorStatus_ = 0;

EDM_ERROR    EDMErrorStatus() {return EDMErrorStatus_;}

#if defined(ERROR_)
#undef ERROR_
#endif
#if defined(ERROR_V)
#undef ERROR_V
#endif

#define  ERROR_(x) edm_error((x),0,0,ProgName)
#define  ERROR_V(x,value) edm_error((x),((long)(value)),0,ProgName)
#define  ERROR_W(x,value1,value2) edm_error((x),((long)(value1)),((long)(value2)),ProgName)

/* const unsigned  MAXSHORT = ((long) 2 << (8 * sizeof(unsigned short) - 1)) - 1; */

void edm_error(int code, long value, long value2, char* module)

{
  EDMErrorStatus_ = (EDM_ERROR) code;
  switch (code)
  {
    case EDM_NOMEMORY:  ErrorMsg ("(%s) Not enough memory to complete operation, aborted.",module); break;
    /* (!!!) for me: exit to DOS, processing delayed */
    case CLUSTERDELTA:  ErrorMsg ("(%s) Cluster size exceed maximal value %u",module,MAXCLUSTERLENGTH); break;
    case WRMODEL:       ErrorMsg ("(%s) Error writing the model. Aborted.",module); break;
    case WRCT:          ErrorMsg ("(%s) Error writing the context tree. Aborted.",module); break;
    case RMODEL:        ErrorMsg ("(%s) Error reading the model. Aborted.",module); break;
    case RCT:           ErrorMsg ("(%s) Error reading the context tree. Aborted.",module); break;
    case QMEXPECTED:    ErrorMsg ("(%s) QMT identification flag not found. Given file isn't model file.",module); break;
    case NOMODEL:       ErrorMsg ("(%s) Model file (.qmt) must be given.",module); break;
    case NOCONTREE:     ErrorMsg ("(%s) Context tree file (.ct) must be given.",module); break;
    case QMTCONTEXT:    ErrorMsg ("(%s) Model file doesn't fit. Differense in CONTEXT SIZE.",module); break;
    case CTCONTEXTS:    ErrorMsg ("(%s) Context tree file doesn't fit. Differense in CONTEXTS.",module); break;
    case CTWRONGCT:     ErrorMsg ("(%s) Error building the C-Tree. %li contexts read, %li expected.",module,value,value2); break;
    case WRONGVERSION:  ErrorMsg ("(%s) Incorrect versions string in EDM file",module); break;
    case STOREBLOCKT_1: ErrorMsg ("(%s) StoreBlockType[%li] out of range.",module,value); break;
    case STOREBLOCKT_2: ErrorMsg ("(%s) StoreBlockType[][%li] out of range.",module,value); break;
    case UNSUPMODEL:    ErrorMsg ("(%s) Model stored in different file. FEATURE UNSUPPORTED!",module); break;
    case UNSUPCT:       ErrorMsg ("(%s) Context-tree is used for compression. FEATURE UNSUPPORTED!",module); break;
    case UNSUPLAYER:    ErrorMsg ("(%s) Decompressing demands extra layer information. FEATURE UNSUPPORTED!",module); break;
    case NOPREVIEW:     ErrorMsg ("(%s) Preview isn not possible. No preview information in file.",module); break;
    case NOCLUSTERS:    ErrorMsg ("(%s) Image isn't devided to clusters.",module); break;
    case WRONGSIZE:     ErrorMsg ("(%s) Output image has wrong size to procede with output.",module); break;

    default:            ErrorMsg ("(%s) Internal error - Call program vendor",module);
  }

#ifndef _WINDOWS
  printf("\n");
#endif

}


/* ================ INFORMATION OUTPUT HANDLING ===================== */
/*                                                                    */
/*                           INFO.C                                   */
/*                                                                    */
/* These functionions serve to output information while such module   */
/* as EDM is running, i.e. these provides output of information about */
/* stage and state of process and/or its progress.                    */
/*                                                                    */
/* Usualy, when process (e.g. compression) starts, special window is  */
/* created, where information output is provided. Information output  */
/* consist of Text Message, which indicates the stage of process, and */
/* Progress Bar, which indicates the state of process.                */
/*                                                                    */
/* Following functions are used to provide information output.        */
/*    (they can be separated into independent module INFO.C)          */
/*                                                                    */
/* . PrintMessage      - print Text Message into Information Window   */
/* . PrintBlankMessage - clean messages space                         */
/* . InitProgressBar   - initialize Progress Bar                      */
/* . UpdateProgressBar - update Progress Bar                          */
/* . DoneProgressBar   - clear Progress Bar                           */
/*                                                                    */
/* Following functions are not yet made (not necessary for DOS/UNIX)  */
/*                                                                    */
/* . InitInfoWindow - init. and draw Information Window               */
/* . DoneInfoWindow - close information window                        */
/*                                                                    */
/* CURRENTLY RELEASED ONLY FOR DOS/UNIX ENVIRONMENT                   */
/* ================================================================== */


#include <stdarg.h>
/* Windows definitions here */
/* #ifdef _WINDOWS          */
/* #endif                   */

int PrintMessage(char *format, ...)
{
        va_list ap;             /* points to each unnamed arg in turn */
        int retval;
        char msg[10000];

        va_start(ap, format);   /* make ap point to 1st unnamed arg */
        vsprintf(msg, format, ap);
        va_end(ap);             /* clean up when done */

  #ifdef _WINDOWS
        /* Windows code here */ return(0);
  #else
        retval = fprintf (stderr, "%s", msg);
        /* fprintf(stderr,"\r"); */
        fflush(stderr);
        return(retval);
  #endif  /* _WINDOWS */
}


int PrintBlankMessage(void)
{
  #ifdef _WINDOWS
        /* Windows code here */ return (0);
  #else
        fprintf(stderr,"\r                                                                      \r");
        fflush(stderr);
        return (0);
  #endif  /* _WINDOWS */
}


int InitProgressBar(void)
{
  #ifdef _WINDOWS
        /* Windows code here */ return (0);
  #else
        fprintf(stderr,"            ");   /* go to the next line to allove InfMsg stay above */
        fflush(stderr);
        return (0);
  #endif  /* _WINDOWS */
}


int UpdateProgressBar(int progress)
{
  static int x;
  char c;

  #ifdef _WINDOWS
        /* Windows code here */ return (0);
  #else

        switch (x)
        {
           case 0: c='-'; x++; break;
           case 1: c='\\'; x++; break;
           case 2: c='|'; x++; break;
           default: c='/'; x=0; break;
        }

        fprintf(stderr,"\b\b\b\b\b\b\b\b\b\b\b\b%c)%3i %% done",c,progress);
        fflush(stderr);
        return (0);
  #endif  /* _WINDOWS */
}


int DoneProgressBar(void)
{
  #ifdef _WINDOWS
        /* Windows code here */ return (0);
  #else
        PrintBlankMessage();
        return (0);
  #endif  /* _WINDOWS */
}

/*                                                                    */
/* ======================== END OF INFO.C =========================== */


/* ================ INFORMATION OUTPUT in EDM ======================= */
/*                                                                    */
/* internal functions which call functions from INFO.C              . */
/* ================================================================== */


typedef enum {
  BLANKMSG=0, ANALYZE=1, ENCODE=2, DECODE=3,
  COMPRESS, DECOMPRESS, PROCESS,
  CODETABLE, HEADERPR, PREVIEWWR,
  STAGE1, STAGE2, STAGE3, STAGE4, STAGEX
  } EDMMESSAGECODE;

/* print information in the window, which remains on the screen      */

void EDMOutput(int code, int value)
{
  if( Para_(QuietLevel) >= 2 )
  {
     switch (code)
     {
        case -1:           break;
        case 0:            PrintBlankMessage(); break;
        case ANALYZE:      PrintMessage ("Analysing:"); break;
        case PROCESS:      PrintMessage ("Processing:"); break;
        case ENCODE:       PrintMessage ("Encoding:"); break;
        case DECODE:       PrintMessage ("Restore image:"); break;
        case CODETABLE:    PrintMessage ("Code table is calculating ..."); break;
        case HEADERPR:     PrintMessage ("Header is processing ..."); break;
        case PREVIEWWR:    PrintMessage ("Preview is writing ..."); break;
        case STAGE1:       PrintMessage ("Stage one:"); break;
        case STAGE2:       PrintMessage ("Stage two:"); break;
        case STAGE3:       PrintMessage ("Stage three:"); break;
        case STAGE4:       PrintMessage ("Stage four:"); break;
        case STAGEX:       PrintMessage ("Stage %i:", value); break;
        default:           PrintMessage ("Processing... Be patient.");
     }
  }
}


/*-------------------------------------------------------------------*/
/* Init Progress Output                                              */

void InitProgress(int code)
{
  if( Para_(QuietLevel) >= 2 )
  {
     EDMOutput(code,0);          /* print the Information Message */
     InitProgressBar();          /* initialize Progress Bar       */
  }
}

/*-------------------------------------------------------------------*/
/* Done Progress Indicator                                           */

void DoneProgress(int code)
{
  if( Para_(QuietLevel) >= 2 )
  {
     EDMOutput(0,0);
     DoneProgressBar();
  }
}

/*-------------------------------------------------------------------*/
/* Print Progress                                                    */

void PrintProgress(int y, int total)
{
  if( Para_(QuietLevel) >= 2 )
  {
     UpdateProgressBar((int) (((long) y * 100L) / total));
  }
}


/* ======================= Global variables ======================== */

/* long    EDMTimer;              */
/* Internal Timer is now desabled */

int  TemplateX[31] = {0,-1,0,-1,+1,-2,0,-2,+2,-1,+1,-2,+2,-3,0,-3,+3,-1,+1,-3,+3,-2,+2,-4, 0,-4,+4,-1,+1,-3,+3 };
int  TemplateY[31] = {0,0,-1,-1,-1,0,-2,-1,-1,-2,-2,-2,-2,0,-3,-1,-1,-3,-3,-2,-2,-3,-3, 0,-4,-1,-1,-4,-4,-3,-3 };

/*    JBIG extended template   */
/*           27 24 28          */
/*     29 21 17 14 18 22 30    */
/*     19 11  9  6 10 12 20    */
/*  25 15  7  3  2  4  8 16 26 */
/*  23 13  5  1  ?             */
/*                             */

int     LayerX[10]     = { 0,0, 0,+1,+1,-1, 0,-1,+1,-1 };
int     LayerY[10]     = { 0,0,+1, 0,+1,+1,-1, 0,-1,-1 };

/*      Layer context setup    */
/*            9  6  8          */
/*            7  1  3          */
/*            5  2  4          */

/* ------------ Global variables / used for Compression ------------ */

/* In the case of serious internal error special error function is   */
/* called from this module to hadle the error situations, default it */
/* output message (through the Errormsg), breaks the program and     */
/* exit to OS, otherwise (non-serious) error,                        */
/* function returns NO and error status can be obtained from the     */
/* EDMErrorStatus()                                                  */

/*=======================  Misc routines  ============================*/

#if defined(max) || defined(min)
#undef max
#undef min
#endif

#define  max(a,b) ((a) > (b) ? (a) : (b))
#define  min(a,b) ((a) < (b) ? (a) : (b))

/*-------------------------------------------------------------------*/

int RoundClusterSize(ClusterSize,BlockSize)
{
  return ( ( (int) ((ClusterSize - 1) / BlockSize) + 1) * BlockSize );
}

/*-------------------------------------------------------------------*/

long EDMSetTimer(long* start)
{
  time_t tmp;
  return( *start = (long)time(&tmp) );
}

/*-------------------------------------------------------------------*/

double EDMCalculateEntropy(double p)
{
  return( p==0 ? 0 : -log(p)/log(2) );
}


/*======================= Q M T - Interface ===========================*/


int ReadQmtHeader(FILE* f)
/* return ContextSize or 0 */
{
  int Id1, Id2, contextsize;

  Id1 = getc(f);
  Id2 = getc(f);
  if( Id1 != 'Q' || Id2 != 'M' )
  {
     ERROR_(QMEXPECTED);
     return(0);
     /* not .QM file or broken */
  }
  fscanf(f, "%i", &contextsize);
  /* Skip white space; only one character can appear! */
  fgetc(f);
  return(contextsize);
}


/*-------------------------------------------------------------------*/


void WriteQmtHeader(FILE* f, int contextsize)
{
  putc('Q', f);
  putc('M', f);
  putc(10, f);
  fprintf(f, "%i", contextsize);  putc(10, f);
}


/*======================= Interface with QM ==========================*/


void GetModel(void)
{
  size_t c;

  switch (Para_(ModellingType))
  {
     case STATIC:       /* never happened */
                        break;

     case SEMIADAPTIVE:
                  EDMOutput(CODETABLE,0);
                  for (c=0;c<Para_(NumberOfStates); c++)
                  {
                     if (Result(ContextUsed[c]) == 0)
                        /* StateIndex[c] = GetStateIndex( 0.5 ); */
                        Para_(StateIndex[c]) = 0; /* the same */
                     else
                        Para_(StateIndex[c]) = GetFirstAttackStateIndex(
                        ( (float) Result(ContextWhites[c]) ) /
                                  Result(ContextUsed[c]) );
                  }
                  EDMOutput(0,0);
                  break;

     default:           break;
  }

}


/*-------------------------------------------------------------------*/


void RestoreModel(void)
{
  size_t c;

  switch (Para_(ModellingType))
  {
     case ADAPTIVE:   NewModel(Para_(NumberOfStates));
                      break;

     case STATIC:
     case SEMIADAPTIVE:
          for (c=0;c<Para_(NumberOfStates); c++)
          {
             RestoreFirstAttackStateIndex(c, Para_(StateIndex[c]));
          }
          break;
  }
}


/*-------------------------------------------------------------------*/


YESNO WriteModel(FILE* f)
{
  BITSTREAM bs;
  size_t c;

  InitializeBitStream(&bs, f);

  for (c=0;c<Para_(NumberOfStates);c++)
  {
     OutputValue(&bs, Para_(StateIndex[c]), 5);
  }

  FlushOutput(&bs);

/*  if (fwrite(Para_(StateIndex),sizeof(BYTE),Para_(NumberOfStates),f)!=Para_(NumberOfStates))
     {ERROR_(WRMODEL); return(NO);}
*/
  return(YES);
}


/*-------------------------------------------------------------------*/


YESNO ReadModel(FILE* f)
{
  BITSTREAM bs;
  size_t c;

  InitializeBitStream(&bs, f);
  FlushInput(&bs);

  for (c=0;c<Para_(NumberOfStates);c++)
  {
     Para_(StateIndex[c]) = InputValue(&bs, 5);
  }

  /* ??? check for ERROR */
  /*if (fread(Para_(StateIndex),sizeof(BYTE),Para_(NumberOfStates),f)!=Para_(NumberOfStates))
     {ERROR_(RMODEL); return(NO);}
  */

  return(YES);
}


/*======================= EDM file routines =========================*/


YESNO ReadEdmTextHeader(FILE* f)
{
  char  VersionString[80];

  fgets(VersionString, 80, f);
  if( strncmp(VersionString,IDString,strlen(IDString)) )
  {
     ERROR_(WRONGVERSION);
     return(NO);
  }
  fscanf(f, "%i, %i\n", &Para_(ImageX), &Para_(ImageY));
  fscanf(f, "%i\n",     &Para_(ModellingType));
  fscanf(f, "%i\n",     &Para_(ContextTreeMode));
  fscanf(f, "%li\n",    &Para_(ContextsInTree));
  fscanf(f, "%i\n",     &Para_(ContextSize));
  fscanf(f, "%i\n",     &Para_(LayerContextSize));
  fscanf(f, "%i\n",     &Para_(Clustering));
  fscanf(f, "%i\n",     &Para_(ClusterSize));
  fscanf(f, "%i\n",     &Para_(BlockModel));
  fscanf(f, "%i\n",     &Para_(BlockContext));
  fscanf(f, "%i\n",     &Para_(BlockSize));

  while( getc(f) == '-' );       /* Read dashes and '\n' */
  Para_(Position) = EDMTEXTHEADER;
  return(YES);

}


/*-------------------------------------------------------------------*/


YESNO WriteEdmTextHeader(FILE* f)
{
  fprintf(f, "%s - %s\n",  IDString, VersionNumber);
  fprintf(f, "%i, %i\n", Para_(ImageX), Para_(ImageY));
  fprintf(f, "%i\n",     Para_(ModellingType));
  fprintf(f, "%i\n",     Para_(ContextTreeMode));
  fprintf(f, "%li\n",    Para_(ContextsInTree));
  fprintf(f, "%i\n",     Para_(ContextSize));
  fprintf(f, "%i\n",     Para_(LayerContextSize));
  fprintf(f, "%i\n",     Para_(Clustering));
  fprintf(f, "%i\n",     Para_(ClusterSize));
  fprintf(f, "%i\n",     Para_(BlockModel));
  fprintf(f, "%i\n",     Para_(BlockContext));
  fprintf(f, "%i\n",     Para_(BlockSize));

  fprintf(f, "--");
  putc(10, f);
  Para_(Position) = EDMTEXTHEADER;
  return(YES);
}


/*-------------------------------------------------------------------*/


YESNO WriteClusterIndexes(FILE* f)
{
  int k;
  long delta, databegin;
  int k0 = 0;
  int n  = Para_(Clustering);

  Para_(IndexesPos) = ftell (f);

  /* let's find first cluster with pixel level codes */
  k=0;
  while ( (Para_(ClusterIndex[k])==NUL || Para_(ClusterIndex[k])==CHESSTYPE)
          && k<n-1 )  k++;
  databegin = Para_(ClusterIndex[k]);
  if( databegin==CHESSTYPE) databegin = NUL;
  fwrite(&databegin,sizeof(long),1,f); /* begin of data */
  if(databegin==NUL) {;} /* all clusters are coded w/o pixel data */

  /* only delta is coded: it defined as
     delta = ClusterIndex[current] - ClusterIndex [previous non-white];
     if cluster is all-white then delta coded as 0
     if cluster is chesstype then as FCLUSTERCODE
     for first non-white cluster delta coded as 1                     */

  if (Para_(Stage) == 2)
  {
     k0 = -1;
     Para_(ClusterIndex[n])=Para_(EndOfData);
     for (k=0;k<n;k++)
     {
        switch (Para_(ClusterIndex[k])) {
        case NUL:             Para_(ClusterDelta[k]) = NUL; break;
        case CHESSTYPE:       Para_(ClusterDelta[k]) = FCLUSTERCODE; break;
        default:
                 delta = k0 < 0 ? 1 : Para_(ClusterIndex[k]) - Para_(ClusterIndex[k0]);
                 if (delta > MAXCLUSTERLENGTH)
                 {
                    ERROR_(CLUSTERDELTA);
                    return(NO); /* cluster to long */
                 }
                 else Para_(ClusterDelta[k]) = (unsigned short) delta;
                 k0 = k;
        }
     }
  }

  /* DEBUG INFO !!! */
  /*
  if( Para_(QuietLevel) >= 3 )
  {
     printf ("Cluster:     Index - Delta\n");
     for (k=0;k<n;k++)
     {  printf (" %6i: %9li - %5hu\n",k,Para_(ClusterIndex[k]),Para_(ClusterDelta[k]));  }
     printf ("    end: %9li\n",Para_(ClusterIndex[n]));
  }
  */

  fwrite(Para_(ClusterDelta),sizeof(unsigned short),n,f);
  fwrite(&(Para_(EndOfData)),sizeof(long),1,f); /* end of data */

  if (Para_(Stage)==1 && Para_(QuietLevel)>=2)
  {
     Result(Cw_) = 0;
     for (k=0;k<n;k++)
     {
        if (Para_(ClusterIndex[k])==NUL) Result(Cw_)++;
     }
  }

  Para_(Position) = EDMCLUSTERINDEX;
  return (YES);

}


/*-------------------------------------------------------------------*/


YESNO ReadClusterIndexes(FILE* f)
{
  int k;
  int k0 = 0;
  long databegin;
  int n  = Para_(Clustering);

  fread(&databegin,sizeof(long),1,f); /* begin of data */

  if(databegin==NUL) ; /* all clusters are w/o pixel codes */

  /* because only delta is coded: ClusterIndex determined as
     ClusterIndex[current] = ClusterIndex [previous non-white] + delta;
     if delta is 0 then cluster is all-white (index is NULL)
     if delta is FCLUSTERCODE then cluster is chess-type (index=1(MIXED))
     for first non-white (==1) delta cluster index yields databegin  */

  fread(Para_(ClusterDelta),sizeof(unsigned short),n,f);
  fread(&(Para_(EndOfData)),sizeof(long),1,f); /* end of data */

  Result(Cw_) = 0; k0 = -1;

  for (k=0;k<n;k++)
  {
     switch (Para_(ClusterDelta[k])) {
     case NUL:            Para_(ClusterIndex[k]) = NUL;  Result(Cw_)++; break;
     case FCLUSTERCODE:   Para_(ClusterIndex[k]) = CHESSTYPE; break;
     default:
                  Para_(ClusterIndex[k]) = k0 < 0 ? databegin :
                  (long) Para_(ClusterDelta[k]) + Para_(ClusterIndex[k0]);
                  k0 = k;
     }
  }
  Para_(ClusterIndex[n])=Para_(EndOfData);

  /* DEBUG INFO!!! */
  /*
  if( Para_(QuietLevel) >= 3 )
  {
     printf ("Cluster:     Index - Delta\n");
     for (k=0;k<n;k++)
     {  printf (" %6i: %9li - %5hu\n",k,Para_(ClusterIndex[k]),Para_(ClusterDelta[k]));  }
     printf ("    end: %9li\n",Para_(ClusterIndex[n]));
  }
  */

  Para_(Position) = EDMCLUSTERINDEX;
  return (YES);

}


/*-------------------------------------------------------------------*/


YESNO WriteBlockCodes(FILE* f)
{
  int x,y,bc;
  long ind;

  ind=ftell(f);
  InitModelQM(18);
  InitEncodeQM();
  SetQMEscMode(1); /* it is default, unnecessary */

  Result(Bw_) = Result(Bm_) = Result(Bb_) = 0;


  for (y=1; y<=Para_(BlocksY); y++)
  {
     for (x=1; x<=Para_(BlocksX); x++)
     {
        bc = Para_(BlockCode[y][x]);
        EncodeBlockType(f,x,y,bc);
        if (Para_(QuietLevel) >= 2 && (Para_(FastMode) == 0) )
        {
           switch (bc) {
           case ALLWHITE: Result(Bw_)++; break;
           case MIXED:    Result(Bm_)++; break;
           case ALLBLACK: Result(Bb_)++; break;
           }
        }
     }
  }

  FlushEncodeQM(f);
  DoneQM();
  Result(BlockDataTotal)+=ftell(f)-ind;
  Para_(Position) = EDMPREVIEW;
  return (YES);

}


/*-------------------------------------------------------------------*/


YESNO ReadBlockCodes(FILE* f)
{
  int x,y,bt;

  InitModelQM(18);
  SetQMEscMode(1);
  InitDecodeQM(f);

  for (y=1; y<=Para_(BlocksY); y++)
  {
     for (x=1; x<=Para_(BlocksX); x++)
     {
        DecodeBlockType(f,x,y,&bt);
        Para_(BlockCode[y][x]) = bt;
     }
  }

  DoneQM();
  Para_(Position) = EDMPREVIEW;
  return (YES);

}


/*-------------------------------------------------------------------*/


YESNO WriteHeader(FILE* f,FILE* ModelFile)
{
   BYTE depth0;

   if (Para_(Position)<EDMTEXTHEADER) return (NO);
   EDMOutput(HEADERPR,0);

   if (Para_(BlockModel)!=PIXELWISE) WriteBlockCodes(f);
   if Para_(Clustering)
   {
      if (WriteClusterIndexes(f)==NO) return(NO);
   }

   if (Para_(ContextTreeMode)==SEMIADAPTIVE_CT)
   {
      /* ??? store the ContextTreeDepth0 (startimg depth) first */
      depth0 = (BYTE) Para_(ContextTreeDepth0);
      if (fwrite(&depth0,sizeof(BYTE),1,f)!=1) {ERROR_(WRCT); return(NO);}
      if (StoreContextTree(Para_(ContextTree),f,Para_(ContextTreeDepth0))==NO) return (NO);
   }

   if (Para_(ModellingType)==SEMIADAPTIVE) /* read from the compressed file */
   {
       if (WriteModel(f)==NO) return(NO);
       if (ModelFile!=NULL)
          if (WriteModel(ModelFile)==NO) return(NO);
          /* write model in the QMT file */
   }

   EDMOutput(0,0);
   Para_(Position) = EDMMODEL;
   return (YES);

}


/*-------------------------------------------------------------------*/


YESNO ReadHeader(FILE* f,FILE* ModelFile, FILE* ContreeFile)
{
   BYTE depth0;
   long contexts;

   if (Para_(Position)<EDMTEXTHEADER) return (NO);

   EDMOutput(HEADERPR,0);

   if (Para_(Position)<EDMPREVIEW)
   {
      if (Para_(BlockModel)!=PIXELWISE)  ReadBlockCodes(f);
      else                               Para_(Position) = EDMPREVIEW;
   }

   if ((Para_(Position)<EDMCLUSTERINDEX))
   {
      if Para_(Clustering)    { ReadClusterIndexes(f);   }
      else                    { Para_(ClusterIndex[0]) = MIXED;
                                Para_(Position) = EDMCLUSTERINDEX; }
   }

   /* position isn't cheked anymore:                                 */
   /* cluster table, CT and model are read as one block together     */

   if (Para_(ContextTreeMode)==STATIC_CT)   /* read from the CT file */
   {
      if ((Para_(ContextTree)=ReadContextTree(ContreeFile,Para_(ContextTreeDepth0),&contexts))==NULL) return (NO);
      if (contexts != Para_(ContextsInTree)) { ERROR_W(CTWRONGCT,contexts,Para_(ContextsInTree)); return(NO); }
   }
   else if (Para_(ContextTreeMode)==SEMIADAPTIVE_CT)   /* read from the QMT file */
   {
      if (Para_(Action)==COMPRESSION)
      {
         if ((Para_(ContextTree)=ReadContextTree(ContreeFile,Para_(ContextTreeDepth0),&contexts))==NULL) return (NO);
         if (contexts != Para_(ContextsInTree)) { ERROR_W(CTWRONGCT,contexts,Para_(ContextsInTree)); return(NO); }
         /* this code is never used */
      }
      else
      {
         if (fread(&depth0,sizeof(BYTE),1,f)!=1) {ERROR_(RCT); return(NO);}
         Para_(ContextTreeDepth0) = (int) depth0;
         /* read CT from EDM-file */
         if ((Para_(ContextTree)=ReadContextTree(f,Para_(ContextTreeDepth0),&contexts))==NULL) return (NO);
         if (contexts != Para_(ContextsInTree)) { ERROR_W(CTWRONGCT,contexts,Para_(ContextsInTree)); return(NO); }
      }
   }
   Para_(Position) = EDMCONTEXTTREE;

   if (Para_(ModellingType)==STATIC)   /* read from the QMT file */
   {
      if (ReadModel(ModelFile)==NO) return (NO);
   }
   else if (Para_(ModellingType)==SEMIADAPTIVE)
   {
      if (ReadModel(f)==NO) return (NO); /* read from the compressed file */
      if (ModelFile!=NULL)
         if (WriteModel(ModelFile)==NO) return(NO); /* write model in the QMT file */
   }
   Para_(Position) = EDMMODEL;

   EDMOutput(0,0);
   return(YES);
}


/*==================== Initializing variables =======================*/


void  ReinitContextInformation(void)
{
  size_t i;

  for(i=0; i<Para_(NumberOfStates); i++)
  {
     Result(ContextUsed[i])    = 0;
     Result(ContextWhites[i])  = 0;
     Result(ContextEntropy[i]) = 0;
  }
  Result(DynamicEntropy) = 0;
}


/*-------------------------------------------------------------------*/


void InitializeEDMStatistics(void)
{

  Result(ContextUsed)    = (long*) allocate(Para_(NumberOfStates)*sizeof(long));
  Result(ContextWhites)  = (long*) allocate(Para_(NumberOfStates)*sizeof(long));
  Result(ContextEntropy) = (double*) allocate(Para_(NumberOfStates)*sizeof(double));

  ReinitContextInformation();

}

/*-------------------------------------------------------------------*/

void InitEDMModel(void)
{
  size_t k,n;

  n = Para_(NumberOfStates);
  if ((Para_(StateIndex) = (BYTE*) allocate( n * sizeof(BYTE) ) ) == NULL)
     { ERROR_(EDM_NOMEMORY); exit(-1); }

  for (k=0;k<n;k++) Para_(StateIndex[k])=0;

}

/*-------------------------------------------------------------------*/

void InitClusterIndex(void)
{
  int k,n;

  n=Para_(Clustering)+1;
  if ((Para_(ClusterIndex) = (size_t*) allocate( n * sizeof(long) ) ) ==NULL) {ERROR_(EDM_NOMEMORY); exit(-1); }
  for (k=0;k<n;k++) Para_(ClusterIndex[k])=MIXED;

  n=Para_(Clustering)+1;
  if ((Para_(ClusterDelta) = (unsigned short*) allocate( n * sizeof(unsigned short) ) ) ==NULL) {ERROR_(EDM_NOMEMORY); exit(-1); }
  for (k=0;k<n;k++) Para_(ClusterDelta[k])=0;
}

/*-------------------------------------------------------------------*/

void InitBlockCodes(void)
{
  int k,i,n;

  if (Para_(BlockModel)!=PIXELWISE)
  {
     if ((Para_(BlockCode) = (PBYTE*) allocate( (Para_(BlocksY)+1) * sizeof(PBYTE) ) )==NULL) {ERROR_(EDM_NOMEMORY); exit(-1); }
     n = (Para_(BlocksX) + 1 );
     for(k=0; k <= (Para_(BlocksY)); k++)
     {
        if ((Para_(BlockCode[k]) = (BYTE*) allocate( n * sizeof(BYTE) ) ) == NULL ) {ERROR_(EDM_NOMEMORY); exit(-1); }
        for (i=0; i<n; i++)
        {
           Para_(BlockCode[k][i]) = ALLWHITE;
        }
     }
  }

}

/*-------------------------------------------------------------------*/

typedef enum {vBASIC=0, vALL=1} EDMVARINITCHOISE;

void InitializeEDMVariables(int flag)

/* intialize 0: basic variables only
             1: all variables (including basic, block codes,
                               model and cluster indecies)
*/

{
  int layerpixels;
  int contextpixels;

  Para_(IndexesPos) = 0L;
  Para_(DataPos)    = 0L;
  Para_(EndOfData)  = 0L;

  layerpixels = Para_(LayerContextSize);

  if (Para_(ContextTreeMode))
  {
     Para_(FirstBlockContext) = Para_(NumberOfStates) =
     ((size_t) Para_(ContextsInTree)) << layerpixels;
  }
  else
  {
     contextpixels = Para_(ContextSize);

     if(contextpixels /* +layerpixels */ > MAXTEMPLATE)
     {
         ErrorMsg("Template size (%i) exceeds recommended maximum (%2i).",
                  layerpixels+contextpixels, MAXTEMPLATE);
     }

     Para_(FirstBlockContext) = Para_(NumberOfStates) =
     1L << (contextpixels+layerpixels);
   }

  Result(ContextUsed)    = NULL;
  Result(ContextWhites)  = NULL;
  Result(ContextEntropy) = NULL;
  Result(DynamicEntropy) = 0;

  Result(Cw_) = 0;             /* number of white clusters */
  Result(Bw_) = 0;             /* number of white blocks   */
  Result(Bb_) = 0;             /* number of black blocks   */
  Result(Bm_) = 0;             /* number if mixed blocks   */

  Result(BytesTotal)     = 0;  /* bytes sent to file in total             */
  Result(PixelDataTotal) = 0;  /* pel coding bytes sended to file         */
  Result(BlockDataTotal) = 0;

/*
  printf("ImageX               = %i \n", Para_(ImageX));
  printf("ImageY               = %i \n", Para_(ImageY));
  printf("ModellingType        = %i \n", Para_(ModellingType));
  printf("BlockModel           = %i \n", Para_(BlockModel));
  printf("BlockContext         = %i \n", Para_(BlockContext));
  printf("BlockSize            = %i \n", Para_(BlockSize));
  printf("Clusters             = %i \n", Para_(Clustering));
  printf("Cluster Size         = %i \n", Para_(ClusterSize));
  printf("ContextSize          = %i \n", Para_(ContextSize));
  printf("LayerContextSize     = %i \n", Para_(LayerContextSize));
  printf("NumberOfStates       = %li \n", Para_(NumberOfStates));
  printf("FirstBlockContext    = %li \n", Para_(FirstBlockContext));
  printf("WholeImage           = %i \n", Para_(WholeImage));
  printf("ForceJBIG            = %i \n", Para_(ForceJBIG));
*/

  if (flag==vALL)
  {
     InitEDMModel();
     InitClusterIndex();
     InitBlockCodes();
  }

}


/*-------------------------------------------------------------------*/


void EDMVariablesDone(void)
{
int k;

  if (Para_(StateIndex)!= NULL) deallocate(Para_(StateIndex));
  if (Para_(ClusterIndex)!= NULL) deallocate(Para_(ClusterIndex));
  if (Para_(ClusterDelta)!= NULL) deallocate(Para_(ClusterDelta));

  if ( Para_(BlockModel)!=PIXELWISE && Para_(BlockCode)!=NULL )
  {
     for(k=0; k <= (Para_(BlocksY)); k++)
     {
        if (Para_(BlockCode[k])!=NULL) deallocate(Para_(BlockCode[k]));
     }
     deallocate(Para_(BlockCode));
  }

  if (Result(ContextUsed)!= NULL)    deallocate(Result(ContextUsed));
  if (Result(ContextWhites)!= NULL)  deallocate(Result(ContextWhites));
  if (Result(ContextEntropy)!= NULL) deallocate(Result(ContextEntropy));

}


/*-------------------------------------------------------------------*/


void InitEncodeParameters(IMAGE* Image)
{

  Para_(ImageX)           = Image->ImageSizeX;
  Para_(ImageY)           = Image->ImageSizeY;

  if (Para_(ForceJBIG))
  {
     /* special JBIG like mode: only pixel level code stream */
     Para_(ModellingType)   = ADAPTIVE;         /* no s/a modelling */
     if (Para_(ContextTreeMode) != STATIC_CT)
         Para_(ContextTreeMode) = NO_CT;          /* no CT modelling  */
     Para_(BlockModel)      = PIXELWISE;        /* no blockcodes    */
     Para_(ClusterSize)     = max(Para_(ImageX),Para_(ImageY));
     Para_(ClusterSize)     = RoundClusterSize(Para_(ClusterSize),Para_(BlockSize));
     Para_(ClustersX)       = Para_(ClustersY) = 1;
     Para_(Clustering)      = 0;                /* no clusters */
  }

  else if ( (Para_(ClusterSize) < Para_(BlockSize)) || (Para_(Clustering)==1) )  /* or WholeImage mode */
  {
     /* only one cluster - whole image mode */
     Para_(ClusterSize) = max(Image->ImageSizeX,Image->ImageSizeY);
     Para_(ClusterSize) = RoundClusterSize(Para_(ClusterSize),Para_(BlockSize));
     Para_(Clustering)  = Para_(ClustersX) = Para_(ClustersY) = 1;
  }

  else /* normal way - with multiple clusters */
  {
     Para_(ClusterSize) = RoundClusterSize(Para_(ClusterSize),Para_(BlockSize));
     /* cluster size must devided to block size without residue */
     /* if not then we round ClusterSize up */

     if (Para_(ClusterSize) > MAXCLUSTERSIZE)
     /* Oops. We exceed the maxima. Let us round down */
     {
        Para_(ClusterSize) -= Para_(BlockSize);
     }

     Para_(ClustersX) = (Para_(ImageX) -1) / Para_(ClusterSize) +1;
     Para_(ClustersY) = (Para_(ImageY) -1) / Para_(ClusterSize) +1;
     Para_(Clustering) = Para_(ClustersX) * Para_(ClustersY);
  }

  /* sophisticated clustersize adjusting scheme */

  if (Para_(Clustering) > MAXCLUSTERS)
  {

    /* To much clusters! Let us enlarge ClusterSize and decrese their quantity */
    ErrorMsg("Clusters to much. Cluster Size will be enlarged.");

    do {

       Para_(ClusterSize) += Para_(BlockSize);

       if (Para_(ClusterSize) > MAXCLUSTERSIZE )
       {
          /* Hm! Now we exceed ClusterSize => clustering doesn'tpossible */
          /* Let us switch to the Whole Image As One Cluster mode        */
          ErrorMsg ("Clustering is not possible. Switch to *whole image* mode.");
          Para_(ClusterSize) = max(Para_(ImageX),Para_(ImageY));
          Para_(ClusterSize) = RoundClusterSize(Para_(ClusterSize),Para_(BlockSize));
          Para_(Clustering)   = Para_(ClustersX) = Para_(ClustersY) = 1;
          break;
       }
       else
       {
          Para_(ClustersX) = (Para_(ImageX) -1) / Para_(ClusterSize) +1;
          Para_(ClustersY) = (Para_(ImageY) -1) / Para_(ClusterSize) +1;
          Para_(Clustering) = Para_(ClustersX) * Para_(ClustersY);
       }

     } while (Para_(Clustering) > MAXCLUSTERS);

  }

  /* we continue here */

  Para_(ImageTotal)  = (long) Para_(ImageX) * (long) Para_(ImageY);

  Para_(BlocksX)     = (Para_(ImageX) -1) / Para_(BlockSize) +1;
  Para_(BlocksY)     = (Para_(ImageY) -1) / Para_(BlockSize) +1;
  Result(Bt_)        = (long) Para_(BlocksX) * Para_(BlocksY);

  Para_(FastMode)    = Para_(QuietLevel)<=1 ? 1 : Para_(FastMode);
  Para_(Stage)       = 1;

}


/*-------------------------------------------------------------------*/


void InitDecodeParameters()
{

  Para_(ClustersX) = (Para_(ImageX) -1) / Para_(ClusterSize) +1;
  Para_(ClustersY) = (Para_(ImageY) -1) / Para_(ClusterSize) +1;

/* !!! it must be that */

  if ( (Para_(Clustering) == 1) &&
       (Para_(Clustering) != Para_(ClustersX) * Para_(ClustersY)) )

  {
    /* old EDM file version */
    Para_(Clustering) = Para_(ClustersX) * Para_(ClustersY);
  }

/*
  assert (Para_(Clustering) == Para_(Clustering) ?
                               Para_(ClustersX) * Para_(ClustersY) : 0 );
*/

  /* Now it must be Cx*Cy or 0 ! */

  Para_(ImageTotal)       = (long) Para_(ImageX) * (long) Para_(ImageY);

  Para_(BlocksX) = (Para_(ImageX) -1) / Para_(BlockSize) +1;
  Para_(BlocksY) = (Para_(ImageY) -1) / Para_(BlockSize) +1;
  Result(Bt_)   = (long) Para_(BlocksX) * Para_(BlocksY);

  Para_(FastMode)         = Para_(QuietLevel)<=1 ? 1 : Para_(FastMode);

  Para_(Stage)            = 1;

}


/*-------------------------------------------------------------------*/


YESNO InitializeContreeFile ( char* ContreeName, FILE** ContreeFile)

{
  long contexts;
  *ContreeFile = NULL;

  if( (Para_(ContextTreeMode)==STATIC_CT) ||
      ( (Para_(ContextTreeMode)==SEMIADAPTIVE_CT) &&
        (Para_(Action)==COMPRESSION) ) )
  {
     if (ContreeName[0]==0)
     {
        ERROR_(NOCONTREE);
        return (NO);
     }
     if( ((*ContreeFile)=FileOpen(ContreeName, INPUT, NO))==NULL )
     {
        return (NO);
     }
     Para_(ContreeFile)=(*ContreeFile);
     if (ReadCTHeader(*ContreeFile,&Para_(ContextTreeDepth0),&contexts)==NO) return(NO);
     /* here we check how much contxets is in EDM-header and in CT-file */
     if (Para_(Action)==COMPRESSION) Para_(ContextsInTree)=contexts;
     else if (Para_(ContextsInTree)!=contexts) {ERROR_(CTCONTEXTS); return(NO);}
  }

  return (YES);

}


/*-------------------------------------------------------------------*/


YESNO InitializeModelFile (char* ModelName, FILE** ModelFile)

{
  int cs;
  *ModelFile = NULL;

  if(Para_(ModellingType)==STATIC)
  {
     if (ModelName[0]==0)
     {
        ERROR_(NOMODEL);
        return (NO);
     }
     if( ((*ModelFile)=FileOpen(ModelName, INPUT, NO))==NULL )
     {
        return (NO);
     }
     Para_(ModelFile)=(*ModelFile);
     cs = ReadQmtHeader(*ModelFile);
     if (cs!=Para_(ContextSize))
     {
        ERROR_(QMTCONTEXT);
        return (NO);
     }
  }
  else if (Para_(ModellingType)==SEMIADAPTIVE && ModelName[0]!=0)
  {
     if ( ((*ModelFile)=FileOpen(ModelName, OUTPUT, Para_(OverWrite)))==NULL)
     {
        return (NO);
     }
     Para_(ModelFile)=(*ModelFile);
     WriteQmtHeader(*ModelFile,Para_(ContextSize));
  }

  return (YES);

}


/*-------------------------------------------------------------------*/


YESNO InitializeEDM  (EDMPARAMETERS* EDMpointer,
                     IMAGE* Image, IMAGE* LayerImage, IMAGE* PreImage,
                     char* InputName, char* OutputName,  char* LayerName,
                     char* PreviewName, char* ModelName, char* ContreeName,
                     FILE** InputFile,  FILE** OutputFile, FILE** ModelFile,
                     FILE** ContreeFile)
/* in the case of error it does Deinitialization function automatically */

{

  int BSize=DEFCLUSTERSIZE+2*CONTEXTHEIGHT;

  /* pass pointer to the EDMPARAMETERS structure */
  PEDMParameters=EDMpointer;

  /* Safety trick */
  Para_(Position)     = NOWHERE;
  Para_(ContextsInTree) = 0L;

  Para_(Image)        = NULL;
  Para_(LayerImage)   = NULL;
  Para_(EDMFile)      = NULL;
  Para_(ModelFile)    = NULL;
  Para_(ContreeFile)  = NULL;
  Para_(StateIndex)   = NULL;
  Para_(ClusterIndex) = NULL;
  Para_(ClusterDelta) = NULL;
  Para_(BlockCode)    = NULL;
  (*InputFile)        = NULL;
  (*OutputFile)       = NULL;

  /* In the case of DJGPP it is by default */
  Image->TargetType = MEMORY;
  Image->FilePointer = NULL;
  Image->data = NULL;
  LayerImage->TargetType = MEMORY;
  LayerImage->FilePointer = NULL;
  LayerImage->data = NULL;
  PreImage->TargetType = MEMORY;
  PreImage->FilePointer = NULL;
  PreImage->data = NULL;

  if( Para_(Action)==COMPRESSION )
  {
     if (ImageInitScrolled(InputName, Image, INPUT, NO, AUTOMATIC, EBB, BSize, CONTEXTHEIGHT, DEFAULTCOLOR)==NO)
        {
            EDMDone(Image, LayerImage, PreImage, *InputFile, *OutputFile, *ModelFile, *ContreeFile);
            return(NO);
        }
     InitEncodeParameters(Image);

     if (Para_(ClusterSize) > DEFCLUSTERSIZE && Para_(Clustering) > 1)
     {
        /* :) You must specife buffers size slightly bigger than cluster size,    */
        /* but until image is opened one can't determine Cluster Size parameter value */
        /* therefore we close image and reopen it !!!                             */
        BSize = Para_(ClusterSize)+2*CONTEXTHEIGHT;
        ImageDone(Image);
        if (ImageInitScrolled(InputName, Image, INPUT, NO, AUTOMATIC, EBB, BSize, CONTEXTHEIGHT, DEFAULTCOLOR)==NO)
        {
            EDMDone(Image, LayerImage, PreImage, *InputFile, *OutputFile, *ModelFile, *ContreeFile);
            return(NO);
        }
     }

     (*InputFile)=Image->FilePointer;
     if (((*OutputFile)=FileOpen(OutputName, OUTPUT, Para_(OverWrite)))==NULL)
         {
            EDMDone(Image, LayerImage, PreImage, *InputFile, *OutputFile, *ModelFile, *ContreeFile);
            return (NO);
         }
     Para_(Image)=Image;
     Para_(EDMFile)=(*OutputFile);
     Para_(Position)=EDMBEGIN;
     if (InitializeContreeFile(ContreeName, ContreeFile) == NO)
        {
            EDMDone(Image, LayerImage, PreImage, *InputFile, *OutputFile, *ModelFile, *ContreeFile);
            return(NO);
        }
     if (InitializeModelFile(ModelName, ModelFile) == NO)
        {
            EDMDone(Image, LayerImage, PreImage, *InputFile, *OutputFile, *ModelFile, *ContreeFile);
            return(NO);
        }
     InitializeEDMVariables(vALL);
     InitializeEDMStatistics(); /* for COMPRESSION only */

     WriteEdmTextHeader(*OutputFile);
  }
  else /* DECOMPRESSION: clear, preview or clusters */
  {
     if (((*InputFile) = FileOpen(InputName, INPUT, NO))==NULL)
         {
            EDMDone(Image, LayerImage, PreImage, *InputFile, *OutputFile, *ModelFile, *ContreeFile);
            return (NO);
         }
     Para_(EDMFile)=(*InputFile);
     Para_(Position)=EDMBEGIN;
     if (ReadEdmTextHeader(*InputFile)==NO) return (NO);
     /* initialization failed due to error in EDM file */
     InitDecodeParameters();

     switch (Para_(Action)) {
     case PREVIEW:
         SetImageSize(Image, Para_(BlocksX), Para_(BlocksY));
         BSize=DEFBUFSIZE;
         break;
     default:
         SetImageSize(Image, Para_(ImageX), Para_(ImageY));
         if (Para_(Clustering) > 1)
         {
            BSize=Para_(ClusterSize) + 2*CONTEXTHEIGHT;
         }
         break;
     }
     Image->FileType = PBM;
     if (InitializeContreeFile(ContreeName, ContreeFile) == NO)
        {
            EDMDone(Image, LayerImage, PreImage, *InputFile, *OutputFile, *ModelFile, *ContreeFile);
            return(NO);
        }
     if (InitializeModelFile(ModelName, ModelFile) == NO)
        {
            EDMDone(Image, LayerImage, PreImage, *InputFile, *OutputFile, *ModelFile, *ContreeFile);
            return(NO);
        }

     InitializeEDMVariables(vALL);
     if (ImageInitScrolled(OutputName, Image, OUTPUT, Para_(OverWrite), AUTOMATIC, EBB, BSize, CONTEXTHEIGHT, DEFAULTCOLOR)==NO)
        {
            EDMDone(Image, LayerImage, PreImage, *InputFile, *OutputFile, *ModelFile, *ContreeFile);
            return(NO);
        }
     (*OutputFile)=Image->FilePointer;
     Para_(Image)=Image;

  }

  if (Para_(Preview) && Para_(Action)!=PREVIEW)
  {
     SetImageSize(PreImage, Para_(BlocksX), Para_(BlocksY));
     PreImage->FileType = PBM;
     if (ImageInit(PreviewName, PreImage, OUTPUT, Para_(OverWrite), AUTOMATIC, EBB, DEFBUFSIZE, DEFAULTCOLOR)==NO)
        {
            EDMDone(Image, LayerImage, PreImage, *InputFile, *OutputFile, *ModelFile, *ContreeFile);
            return(NO);
        }
  }

  if((Para_(LayerContextSize)>0) && Para_(Action)!=PREVIEW)
  {
     if (ImageInitScrolled(LayerName, LayerImage, INPUT, NO, AUTOMATIC, EBB, BSize, CONTEXTHEIGHT, DEFAULTCOLOR)==NO)
        {
            EDMDone(Image, LayerImage, PreImage, *InputFile, *OutputFile, *ModelFile, *ContreeFile);
            return(NO);
        }
     Para_(LayerImage)=LayerImage;
  }

return (YES);

}


/*-------------------------------------------------------------------*/


YESNO EDMDone(IMAGE* Image, IMAGE* LayerImage, IMAGE* PreImage,
              FILE* InputFile,  FILE* OutputFile, FILE* ModelFile,
              FILE* ContreeFile)
{
  Para_(Position)=NOWHERE;
  ImageDone(Image);
  if ( (Para_(Action)==COMPRESSION )&&( (OutputFile!=NULL)) ) fclose(OutputFile);
  if ( ((Para_(Action)==DECOMPRESSION ) || (Para_(Action)==PREVIEW ))
       &&( (InputFile!=NULL)) ) fclose(InputFile);
  if(Para_(LayerContextSize)>0)
  {
     ImageDone(LayerImage);
  }
  if (Para_(Preview) && Para_(Action)!=PREVIEW)
  {
     ImageDone(PreImage);
  }

  if (ContreeFile!=NULL) fclose(ContreeFile);
  if (ModelFile!=NULL) fclose(ModelFile);

  EDMVariablesDone();

  Para_(Image)      = NULL;
  Para_(LayerImage) = NULL;
  Para_(EDMFile)    = NULL;
  Para_(ModelFile)  = NULL;
  Para_(ContreeFile)= NULL;

  return (YES); /* no problem at all */

}


/*==================== Write Preview Image ==========================*/


YESNO WritePreview(IMAGE* Image)
{
int x,y,pel,bt;

  if (Para_(BlockModel)==PIXELWISE) return(NO);

  for (y=1; y<=Para_(BlocksY); y++)
  {
     for (x=1; x<=Para_(BlocksX); x++)
     {
        bt = Para_(BlockCode[y][x]);
        pel = (bt == ALLWHITE) ? WHITE : BLACK;
        PutImageBitPixel(Image,x,y,pel);
     }
  }

  return (YES);
}


/*-------------------------------------------------------------------*/


YESNO PreviewImage(IMAGE* Image, FILE* InputFile)
/* return YES if O.K., NO - if preview isn't supported */

{

  if (Para_(BlockModel)!=PIXELWISE)
  {
     if (ReadHeader(InputFile,NULL,NULL)==NO) return(NO);
     EDMOutput(PREVIEWWR,0);
     WritePreview(Image);
     EDMOutput(0,0);
     return(YES);
  }

  ErrorMsg ("Pixelwise modelling. Preview does not supported.");
  return(NO);

}



/*======================= Context modelling =========================*/


void Update(int bit, size_t c)
/* !!! call it for compression only, otherwise STATISTICS are not initialized */

{
  double e;


/*   if( Para_(QuietLevel)>=5 )
       {
          printf("Update: bit=%i, context=%li. \n",bit,c);
       }
*/

  if( Para_(Stage)==1 )
  ;  /* no encoding here -> no dynamic entropy */
  else
  {

 /*  if(c<0 || c>=Para_(NumberOfStates))
       {
          printf("ERROR: Update(%i,%li) where max = %li.\n",bit,c,Para_(NumberOfStates));
       }
 */
    if ( Para_(FastMode) == 0 )
    {
       /* only at the stage 2 and not in the fast mode */
       e = EDMCalculateEntropy(GetProbabilityQM(bit,(int) c));
       Result(DynamicEntropy) += e;
       Result(ContextEntropy[c]) += e;
    }
  }

  Result(ContextUsed[c])++;
  if(bit==WHITE) Result(ContextWhites[c])++;
}


/*-------------------------------------------------------------------*/

size_t PixelContext(IMAGE* Image, IMAGE* LayerImage, int x, int y)
{
  int  i,bit;
  size_t c;
  int  x_,y_;

  /* For macro */
  int by,bx,tmp;
  BYTE value,mask;

  if (Para_(ContextTreeMode)==NO_CT)
  {
     bit = GetImageBitPixel(Image, x+TemplateX[1], y+TemplateY[1] );
     c = (size_t) bit;

     for(i=2; i<=Para_(ContextSize); i++)
     {
        x_  = x+TemplateX[i];
        y_  = y+TemplateY[i];
        GetImageBitPixel_(bit,Image,x_,y_)
        c = (c<<1) + bit;
     }
  }
  else
  {
     c = (size_t) GetVariableContext(Para_(ContextTree), Image, x, y);
  }

  for(i=1; i<=Para_(LayerContextSize); i++)
  {
     bit = GetImageBitPixel( LayerImage, x+LayerX[i], y+LayerY[i] );
     c = (c<<1) + bit;
  }

  /*
  if( Para_(QuietLevel) >= 5 )
     printf("Context for (%i,%i) = %i \n",x,y,c);
  */

  return(c);
}



/* --- Context pixels are taken in the order: ---- */
/*            .  . 21 17 14 18 22  .  .            */
/*            . 19 11  9  6 10 12 20  .            */
/*            . 15  7  3  2  4  8 16  .            */
/*            . 13  5  1  X  .  .  .  .            */
/* ---- Layer pixels are taken in the order: ----- */
/*            .  .  .  .  5  .  .  .  .            */
/*            .  .  .  4  1  2  .  .  .            */
/*            .  .  .  .  3  .  .  .  .            */
/* ----------------------------------------------- */


/*================= Blockwise Context Modelling =====================*/


int BlockContext(int blockindex, int blockline)
{
  int left  = Para_(BlockCode[blockline][blockindex-1]);
  int above = Para_(BlockCode[blockline-1][blockindex]);

  switch( Para_(BlockContext) )
     {
     case MODEL_0: return( 0 );
     case MODEL_1: return( 2*left );
     case MODEL_2: return( 2*left + 6*above);
     default:      return( 0 );
     }
}

/* -----------------MODEL_2: --------------- */
/* above:    W  W  W    B  B  B    M  M  M   */
/* left:     W  B  M    W  B  M    W  B  M   */
/* RESULT:   0  2  4    6  8 10   12 14 16   */
/* ----------------------------------------- */
/* Note that two decisions might be needed   */
/* for each context: ALLWHITE/ALLBLACK/MIXED */


/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
/*                                                                   */
/*               ššš   šš   ššš   š  šš  š   ššš                     */
/*              š     š  š  š  š  š  š š š  š                        */
/*              š     š  š  š  š  š  š š š  š šš                     */
/*              š     š  š  š  š  š  š  šš  š  š                     */
/*               ššš   šš   ššš   š  š  šš   ššš                     */
/*                                                                   */
/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */


/*======================== Block modelling ==========================*/

int AllWhite(IMAGE* Image, int FirstX, int FirstY, int blocksize)
{
  int  x,y;
/*
  int  LastX = min( FirstX-1+blocksize, Para_(ImageX) );
  int  LastY = min( FirstY-1+blocksize, Para_(ImageY) );
*/
  int  LastX = min( FirstX-1+blocksize, Image->Bound.Right );   /* ??? */
  int  LastY = min( FirstY-1+blocksize, Image->Bound.Lower );

  for(y=FirstY ; y<=LastY ; y++)
    for(x=FirstX ; x<=LastX ; x++)
       {
       if( GetImageBitPixel(Image,x,y)==BLACK) return(MIXED);
       }

  return(ALLWHITE);
}


/*-------------------------------------------------------------------*/


int GetBlockType(IMAGE* Image, int FirstX, int FirstY, int blocksize)
{
  int  x,y;

/*
  int  LastX = min( FirstX-1+blocksize, Para_(ImageX) );
  int  LastY = min( FirstY-1+blocksize, Para_(ImageY) );
*/
  int  LastX = min( FirstX-1+blocksize, Image->Bound.Right );   /* ??? */
  int  LastY = min( FirstY-1+blocksize, Image->Bound.Lower );

  int  whites = 0;
  int  blacks = 0;

  for(y=FirstY ; y<=LastY ; y++)
    for(x=FirstX ; x<=LastX ; x++)
       {
       if( GetImageBitPixel(Image,x,y)==BLACK) blacks++; else whites++;
       if(blacks && whites) return(MIXED);
       }

  if(blacks) return(ALLBLACK);
  else       return(ALLWHITE);
}


/*-------------------------------------------------------------------*/


void StoreBlockType(int blockindex, int blockline, int blocktype)
{
  /* for debugging only - cut it */
  /*
  if(blockline<0 || blockline>Para_(BlocksY))
  {
     ERROR_V(STOREBLOCKT_1,blockline);
     exit(-1);
  }
  if(blockindex<0 || blockindex>Para_(BlocksX))
  {
     ERROR_V(STOREBLOCKT_2,blockindex);
     exit(-1);
  }
  */
  Para_(BlockCode[blockline][blockindex]) = blocktype;
}


/*-------------------------------------------------------------------*/


void EncodeBlockType(FILE* f, int blockindex, int blockline, int blocktype)
{
  int  c;
  int  bit;

  c = BlockContext(blockindex, blockline);
  if( Para_(BlockModel)==BLOCK )
  {
     bit = (blocktype==ALLWHITE ? 0 : 1);
     EncodeBitByQM(f, c, bit);
  }
  else switch(blocktype)
  {
     case ALLWHITE:  EncodeBitByQM(f, c,   0); break;
     case ALLBLACK:  EncodeBitByQM(f, c,   1);
                     EncodeBitByQM(f, c+1, 0); break;
     case MIXED:     EncodeBitByQM(f, c,   1);
                     EncodeBitByQM(f, c+1, 1); break;
     default:        break;
  }
}


/*-------------------------------------------------------------------*/


int AnalyseBlocks(IMAGE* Image, int FirstY, int blockline)
{
  int x,bt;
  int b = (Image->Bound.Left - 1) / Para_(BlockSize) + 1;
  int nwc = 0;

  for(x=Image->Bound.Left; x<=Image->Bound.Right ; x+=Para_(BlockSize))
  {
    if(Para_(BlockModel)==BLOCK || Para_(BlockModel)==PIXELWISE)
    /* block analysis is needed for cluster analysis */
    {
       bt = AllWhite(Image,x,FirstY,Para_(BlockSize));
    }
    else
    {
       bt = GetBlockType(Image,x,FirstY,Para_(BlockSize));
    }
    if (bt != ALLWHITE) nwc++;
    if (Para_(BlockModel)!=PIXELWISE) StoreBlockType(b++,blockline,bt);
  }
  return (nwc);
  /* return number of the non-white blocks in the srip */
}


/*-------------------------------------------------------------------*/


void EncodeSlize(FILE* f, IMAGE* Image, IMAGE* LayerImage,
                 int FirstX, int y, int blocksize)
{
  int  x,bit=0;
  size_t c;
  int  LastX = min(FirstX+blocksize-1,Image->Bound.Right);

  for(x=FirstX; x<=LastX; x++)
  {
     bit = GetImageBitPixel(Image,x,y);
     c = PixelContext(Image,LayerImage,x,y);

     Update(bit,c);

     if( (Para_(Stage)==2) )
     {
        EncodeBitByQM(f, (int) c, bit);   /* encoding only at stage 2 */
     }

  }
}


/*-------------------------------------------------------------------*/


int EncodePixels(FILE* f, IMAGE* Image, IMAGE* LayerImage,
                 int FirstY, int blockline, int blocksize)
/* return TRUE (1) if at least one pel is coded, FALSE otherwise */
{

  int  y,b,FirstX;
  int  LastY = min( FirstY-1+blocksize, Image->Bound.Lower );
  int  wascoded=0;

  if (Para_(BlockModel)==PIXELWISE) wascoded = 1;

  /* Process the blocks row by row */
  for(y=FirstY ; y<=LastY ; y++)
  {
     /* Process the blocks slize by slize */
     b = (Image->Bound.Left - 1) / Para_(BlockSize) + 1;
     for(FirstX=Image->Bound.Left; FirstX<=Image->Bound.Right ; FirstX+=blocksize, b++)
     {
        /* Only slizes of mixed blocks are encoded */
        if (Para_(BlockModel)==PIXELWISE)
        {
           EncodeSlize(f,Image,LayerImage,FirstX,y,blocksize);
        }
        else if (Para_(BlockCode[blockline][b])==MIXED)
        {
           EncodeSlize(f,Image,LayerImage,FirstX,y,blocksize);
           wascoded = 1;
        }
     }
  }
  return (wascoded);
}


/*-------------------------------------------------------------------*/


void AnalyseImage(IMAGE* Image, IMAGE* LayerImage, FILE* ModelFile)
{
  int y, blockline, ci, cx, cy;
  int nwc_, nwc = 0;

  /* initializing QM */

  Para_(Stage) = 1;
  InitProgress(ANALYZE);

  ci=0;
  for (cy=1; cy<=Para_(ImageY); cy+=Para_(ClusterSize))
  {
     for (cx=1; cx<=Para_(ImageX); cx+=Para_(ClusterSize))
     {
        SetImageBounds(Image,cx,cy,Para_(ClusterSize),Para_(ClusterSize));

        /* Cluster processing */
        blockline = (Image->Bound.Upper - 1) / Para_(BlockSize) + 1;
        for(y=Image->Bound.Upper; y<=Image->Bound.Lower; y+=Para_(BlockSize), blockline++)
        {
           nwc_ = AnalyseBlocks(Image, y, blockline);
           nwc += nwc_;
           if (Para_(Clustering) <= 1) PrintProgress(y,Para_(ImageY));

           if (Para_(ModellingType) == SEMIADAPTIVE) /* update() needed */
           {
              if ( (Para_(BlockModel)==PIXELWISE) || (nwc_ != 0) )
              {
                 EncodePixels(NULL, Image, LayerImage, y, blockline, Para_(BlockSize));
              }
           }
        }
        if (Para_(Clustering) > 1) PrintProgress(ci+1,Para_(Clustering));
        if (nwc == 0) Para_(ClusterIndex[ci]) = NUL;
        else          Para_(ClusterIndex[ci]) = MIXED;

        ci++; nwc=0;
     }
  }

  if (Para_(ModellingType)==SEMIADAPTIVE) /* get model */
  {
     GetModel();
  }
  else if (Para_(ModellingType)==STATIC) /* read from the QMT file */
  {
     ReadModel(ModelFile);
     /* ??? if error happened ??? */
  }

  DoneProgress(ANALYZE);

}


/*-------------------------------------------------------------------*/


void CompressImage(IMAGE* Image, FILE* OutputFile, IMAGE* LayerImage)
{
  int y, blockline, ci, cx, cy, ci_, wascoded;
  int EscMode;

  Para_(Stage) = 2;
  InitProgress(ENCODE);

  /* initializing QM */
  InitModelQM(Para_(NumberOfStates));

  ci=0; ci_=0;
  for (cy=1; cy<=Para_(ImageY); cy+=Para_(ClusterSize))
  {
     for (cx=1; cx<=Para_(ImageX); cx+=Para_(ClusterSize))
     {
        SetImageBounds(Image,cx,cy,Para_(ClusterSize),Para_(ClusterSize));
        wascoded = 0;
        if (Para_(ClusterIndex[ci]) != NUL)
        {
           Para_(ClusterIndex[ci]) = ftell(OutputFile);

           InitEncodeQM();
           RestoreModel();
           EscMode = !(Para_(Clustering)); /* original QM for original JBIG */
           SetQMEscMode(EscMode);

           blockline = (Image->Bound.Upper - 1) / Para_(BlockSize) + 1;
           for(y=Image->Bound.Upper; y<=Image->Bound.Lower; y+=Para_(BlockSize), blockline++)
           {
              if (Para_(Clustering) <= 1) PrintProgress(y, Para_(ImageY));
              wascoded |= EncodePixels(OutputFile, Image, LayerImage,
                                       y, blockline, Para_(BlockSize));
           }
           FlushEncodeQM(OutputFile);
           if (!wascoded)
           {
              /* used assumption that if no bits were QM-coded then   */
              /* no bytes will be written in FlushEncodeQM() function */

              assert (Para_(ClusterIndex[ci]) == ftell(OutputFile));

              Para_(ClusterIndex[ci]) = CHESSTYPE;
           }

        }

        Para_(Position) = EDMPIXELDATA;
        if (Para_(ClusterIndex[ci++]) != NUL)
           if (Para_(Clustering) > 1)
               PrintProgress(++ci_, (Para_(Clustering)-Result(Cw_)));
     }
  }

  DoneQM();

  DoneProgress(ENCODE);
}


/*-------------------------------------------------------------------*/


YESNO EncodeImage(IMAGE* Image, FILE* OutputFile, IMAGE* LayerImage,
                 IMAGE* PreImage, FILE* ModelFile, FILE* ContreeFile)
{
  long contexts;
  /* EDMSetTimer(&EDMTimer); */

  if ( Para_(ContextTreeMode)==STATIC_CT ||
       Para_(ContextTreeMode)==SEMIADAPTIVE_CT )
  {
     if ((Para_(ContextTree)=ReadContextTree(ContreeFile,Para_(ContextTreeDepth0),&contexts))==NULL) return (NO);
     if (contexts != Para_(ContextsInTree)) { ERROR_W(CTWRONGCT,contexts,Para_(ContextsInTree)); return(NO); }
  }

  if ( Para_(ModellingType)==ADAPTIVE &&         /* no s/a modelling */
       Para_(BlockModel)==PIXELWISE &&           /* no blockcodes    */
       Para_(Clustering)<=1 )            /* 1 cluster or no clusters */

  {
       if (Para_(Clustering)!=0)
       {
          Para_(ClusterIndex[0]) = MIXED;
          if (WriteHeader(OutputFile,NULL)==NO) return (NO);

          /* one cluster index (non-white) is written */
          /* for compatibility with EDM file format   */
       }

       /* so, if JBIG mode then no cluster codes */
       goto Stage2;
  }

  /* Stage I */

  AnalyseImage(Image, LayerImage, ModelFile);
  if (WriteHeader(OutputFile,ModelFile)==NO) return(NO);
  if (Para_(Preview))
  {
     WritePreview(PreImage);
  }

  /* if (Para_(QuietLevel) > 1) EDMPrintTime("Analysing time",EDMTimer); */

  /* Reinitializing */
  ReinitContextInformation();
  ImageRewind(Image);
  if (Para_(LayerContextSize) >0 && Para_(ModellingType)==SEMIADAPTIVE)
  {
     ImageRewind(LayerImage);
  }

  /*  Stage II */

; Stage2: ;
  /* EDMSetTimer(&EDMTimer); */
  Para_(DataPos)=ftell(OutputFile);

  CompressImage(Image, OutputFile, LayerImage);
  /* if (Para_(QuietLevel) > 1) EDMPrintTime("Encoding time",EDMTimer); */

  /* epilog */

  fflush(OutputFile);
  Para_(EndOfData)=ftell(OutputFile);
  Result(BytesTotal)=Para_(EndOfData);
  Result(PixelDataTotal)=Result(BytesTotal)-Para_(DataPos);

  /* write actutal cluster indexes */

  if (Para_(Clustering))
  {
     fseek(OutputFile,Para_(IndexesPos),SEEK_SET);
     WriteClusterIndexes(OutputFile);
  }

  return (YES);
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
/*                                                                   */
/*          ššš   šššš   ššš   šš   ššš   š  šš  š   ššš             */
/*          š  š  š     š     š  š  š  š  š  š š š  š                */
/*          š  š  ššš   š     š  š  š  š  š  š š š  š šš             */
/*          š  š  š     š     š  š  š  š  š  š  šš  š  š             */
/*          ššš   šššš   ššš   šš   ššš   š  š  šš   ššš             */
/*                                                                   */
/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */


/*======================== Block modelling ==========================*/


void DecodeBlockType(FILE* f, int blockindex, int blockline, int* blocktype)
{
  int  bit;
  int  c;

  c = BlockContext(blockindex,blockline);
  if( Para_(BlockModel)==BLOCK )
  {
     bit = DecodeBitByQM(f, (int) c);
     assert(bit==0 || bit==1);
     (*blocktype) = bit ? MIXED : ALLWHITE;
  }
  else
  {
     bit = DecodeBitByQM(f, c);
     assert(bit==0 || bit==1);
     if(bit) (*blocktype)  = DecodeBitByQM(f, c+1) ? MIXED : ALLBLACK;
     else    (*blocktype)  = ALLWHITE;
  }
}


/*======================== Pixel modelling ==========================*/


void MakeUniformSlize(IMAGE* Image, int FirstX, int y, int blocksize, int color)
{
  int  LastX = min( FirstX-1+blocksize, Image->Bound.Right );
  int  x;

  for(x=FirstX; x<=LastX; x++)
     {
     PutImageBitPixel(Image, x, y, color);
     }
}


/*-------------------------------------------------------------------*/


void DecodeSlize(FILE* f, IMAGE* Image,  IMAGE* LayerImage,
                 int FirstX, int y, int blocksize)
{
  size_t  c;
  int  x,bit;
  int  LastX = min(FirstX+blocksize-1,Image->Bound.Right);

  for(x=FirstX; x<=LastX; x++)
     {
     c = PixelContext(Image,LayerImage,x,y);
     bit = DecodeBitByQM(f, (int) c);
     assert(bit==0 || bit==1);
     /* !!! Update(bit,c); */
     PutImageBitPixel(Image, x, y, bit);
     }
}


/*-------------------------------------------------------------------*/


void DecodePixels(FILE* f, IMAGE* Image, IMAGE* LayerImage,
                  int FirstY, int blockline, int blockpos, int blocksize)
{
  int  y,b,FirstX;
  int  LastY = min( FirstY-1+blocksize, Image->Bound.Lower );

  /* Process the blocks row by row */
  for(y=FirstY ; y<=LastY ; y++)
  {
     /* Decode slize by slize */

     /* b = (Image->Bound.Left - 1) / Para_(BlockSize) + 1; */
     b = blockpos;
     for(FirstX=Image->Bound.Left; FirstX<=Image->Bound.Right ; FirstX+=blocksize, b++)
     {
         /* Only slizes of mixed blocks are decoded */

         if ((Para_(BlockModel)==PIXELWISE))
         {
             DecodeSlize(f,Image,LayerImage,FirstX,y,blocksize);
         }
         else
         {
             switch (Para_(BlockCode[blockline][b]))
             {
                case ALLWHITE: /* MakeUniformSlize(Image,FirstX,y,blocksize,WHITE); */ break;
                     /* !!! it is by default */
                case MIXED:    DecodeSlize(f,Image,LayerImage,FirstX,y,blocksize); break;
                default: /* case ALLBLACK: */ MakeUniformSlize(Image,FirstX,y,blocksize,BLACK); break;
             }
         }
     }
  }
}


/*-------------------------------------------------------------------*/


YESNO DecodeImage(IMAGE* Image, FILE* InputFile, IMAGE* LayerImage,
                 IMAGE* PreImage, FILE* ModelFile, FILE* ContreeFile)
{
  int y,blockline,blockpos,ci,cx,cy,ci_,ci__;
  int EscMode;
  long EndCodePos;

  if (ReadHeader(InputFile,ModelFile,ContreeFile)==NO) return(NO);
  if (Para_(Preview))
  {
     WritePreview(PreImage);
  }

  InitProgress(DECODE);
  InitModelQM(Para_(NumberOfStates));

  /* Reinitializing */

  ci=0; ci_=0;
  for (cy=1; cy<=Para_(ImageY); cy+=Para_(ClusterSize))
  {
     for (cx=1; cx<=Para_(ImageX); cx+=Para_(ClusterSize))
     {
        SetImageBounds(Image,cx,cy,Para_(ClusterSize),Para_(ClusterSize));

        if (Para_(ClusterIndex[ci]) != NUL)
        {

        /* if (ClusterIndex [ci] != CHESSTYPE)          */
        /* fseek (InputFile,ClusterIndex[ci],SEEK_SET); */
        /* in sequential using it is not necessary      */
        /* if (Para_(Clustering)==0) you must avoid it  */

           if (Para_(ClusterIndex[ci]) != CHESSTYPE)
           /* pixels aren't decoded from code stream, only blocks */
           {
              RestoreModel();
              EscMode = !(Para_(Clustering)); /* original QM for original JBIG */
              SetQMEscMode(EscMode);
              if (!EscMode)
              {
                 ci__=ci+1;
                 EndCodePos=Para_(ClusterIndex[ci__]);
                 while (EndCodePos==NUL || EndCodePos==CHESSTYPE)
                 {
                    ci__++;
                    assert (ci__<=Para_(Clustering));
                    /* cluster ci__ out of range.                        */
                    /* other is IMPOSSIBLE (!!!),                        */
                    /* because ClusterIndex[Para_(ClustersT)] != NULL .  */
                    EndCodePos = Para_(ClusterIndex[ci__]);        /* next one */
                 }
                 SetQMEndCodePos(EndCodePos);
              }
              InitDecodeQM(InputFile);
           }

           blockline = (Image->Bound.Upper - 1) / Para_(BlockSize) + 1;
           blockpos = (Image->Bound.Left - 1) / Para_(BlockSize) + 1;
           for(y=Image->Bound.Upper; y<=Image->Bound.Lower; y+=Para_(BlockSize), blockline++)
           {
              if (Para_(Clustering) <= 1) PrintProgress(y, Para_(ImageY));
              DecodePixels(InputFile, Image, LayerImage, y, blockline, blockpos, Para_(BlockSize));
           }

        }
        /* else cluster is white (white by default - no problem) */

        Para_(Position) = EDMPIXELDATA;
        if (Para_(ClusterIndex[ci++]) != NUL)
           if (Para_(Clustering) > 1)
               PrintProgress(++ci_, (Para_(Clustering)-Result(Cw_)));

     }
  }

  DoneQM();

  if ( Para_(ContextTreeMode)==STATIC_CT ||
       Para_(ContextTreeMode)==SEMIADAPTIVE_CT )
  {
     DeallocateContextTree(Para_(ContextTree));
     Para_(ContextTree) = NULL;
  }

  DoneProgress(DECODE);
  Para_(Position) = EDMEND;
  return (YES);

}


/*============================ INTERFACE =============================*/
/*                                                                    */
/*   INTERFACE FUNCTIONS FOR EDM-VIEW SOFTWARE                        */
/*                                                                    */

YESNO OpenEDMImage(EDMPARAMETERS* EDMpointer,
                   char* EDMName,
                /* char* LayerName, */
                /* char* ModelName, */
                   int QuietLevel)

{
  YESNO retstatus=YES;

  PEDMParameters = EDMpointer;

  Para_(Position)     = NOWHERE;
  Para_(Image)        = NULL;
  Para_(LayerImage)   = NULL;
  Para_(EDMFile)      = NULL;
  Para_(ModelFile)    = NULL;
  Para_(ContreeFile)  = NULL;

  Para_(FastMode)         = 1;
  Para_(Action)           = DECOMPRESSION;
  Para_(Preview)          = 0;
  Para_(OverWrite)        = 0;
  Para_(QuietLevel)       = QuietLevel;

  Para_(StateIndex)       = NULL;
  Para_(ClusterIndex)     = NULL;
  Para_(ClusterDelta)     = NULL;
  Para_(BlockCode)        = NULL;

  if( (Para_(EDMFile) = FileOpen(EDMName, INPUT, NO)) == NULL) return(NO);
  Para_(Position)   = EDMBEGIN;
  /* initialization failed due to error in EDM file */

  if (ReadEdmTextHeader(Para_(EDMFile))==NO) return (NO);
  /* initialization failed due to error in EDM file */

  InitDecodeParameters();
  /* initialize the rest of parameters on the base of Text Header */

  InitializeEDMVariables(vBASIC);

  return(retstatus);

}


/*-------------------------------------------------------------------*/


YESNO CloseEDMImage(EDMPARAMETERS* EDMpointer)
{
  YESNO retstatus=YES;

  PEDMParameters = EDMpointer;
  Para_(Position)   = NOWHERE;

  fclose(Para_(EDMFile));

  if( (Para_(LayerContextSize)>0) && (Para_(LayerImage)!=NULL) )
  {
     ImageDone(Para_(LayerImage));
     deallocate(Para_(LayerImage));
  }

  if (Para_(ModelFile)!=NULL) fclose(Para_(ModelFile));
  if (Para_(ContreeFile)!=NULL) fclose(Para_(ContreeFile));

  EDMVariablesDone();

  Para_(Image)        = NULL;
  Para_(LayerImage)   = NULL;
  Para_(EDMFile)      = NULL;
  Para_(ModelFile)    = NULL;
  Para_(ContreeFile)  = NULL;

  return (retstatus);
}


/*-------------------------------------------------------------------*/


YESNO ReadEDMCluster(EDMPARAMETERS* EDMpointer, int ClusterX, int ClusterY, IMAGE* Image, int PosX, int PosY)
/* Decompress one cluster and store it to the given IMAGE.                */
/* HUOM! Only clustered EDM images are supported!                         */
/*       Otherwise use DecodeImage() to decode whole image.               */
/* Image must be initialised and have aproriate sizes.                    */
/* Image can be "in memory only" - see BINIMAGE.h for details.            */
/* Return NO if clustering doesn't supported or output image is to small, */
/* check EDMErrorStatus() to identify the error.                          */
/* HUOM! Image must be cleaned (to white) before this operation.          */
{
  int FirstX, FirstY, LastX, LastY, DeltaX, DeltaY;
  int blockline, blockpos, y, ci, ci__;
  long EndCodePos;
  int EscMode;
  YESNO retstatus=YES;

  if (Para_(Position)<EDMMODEL) ReadEDMImageData(EDMpointer);
  /* initialize and read necessary data first */

  PEDMParameters = EDMpointer;

  if (Para_(Clustering)==0) {EDMErrorStatus_=NOCLUSTERS; return (NO); }

  FirstX = (ClusterX-1) * Para_(ClusterSize) + 1;
  FirstY = (ClusterY-1) * Para_(ClusterSize) + 1;
  LastX  = min((ClusterX * Para_(ClusterSize)),Para_(ImageX));
  LastY  = min((ClusterY * Para_(ClusterSize)),Para_(ImageY));
  DeltaX = LastX - FirstX + 1;
  DeltaY = LastY - FirstY + 1;

  ci = Para_(ClustersX) * (ClusterY - 1) + ClusterX - 1;

  Image->Bound.Left  = PosX;
  Image->Bound.Upper = PosY;
  if ( (Image->Bound.Right = PosX + DeltaX - 1) > Image->ImageSizeX) { EDMErrorStatus_=WRONGSIZE; return (NO);}
  if ( (Image->Bound.Lower = PosY + DeltaY - 1) > Image->ImageSizeY) { EDMErrorStatus_=WRONGSIZE; return (NO);}

  /* Image is prepared - now we begin decoding process */

  if (Para_(ClusterIndex[ci]) != NUL)
  {

     InitModelQM(Para_(NumberOfStates));

     if (Para_(ClusterIndex[ci]) != CHESSTYPE)
     /* pixels aren't decoded from code stream, only blocks */
     {
        fseek (Para_(EDMFile),Para_(ClusterIndex[ci]),SEEK_SET);
        RestoreModel();
        EscMode = !(Para_(Clustering)); /* original QM for original JBIG */
        SetQMEscMode(EscMode);
        if (!EscMode)
        {
           ci__=ci+1;
           EndCodePos = Para_(ClusterIndex[ci__]);
           while (EndCodePos==NUL || EndCodePos==CHESSTYPE)
           {
              ci__++;
              assert (ci__<=Para_(Clustering));
              EndCodePos = Para_(ClusterIndex[ci__]);        /* next one */
           }
           SetQMEndCodePos(EndCodePos);

        }
        InitDecodeQM(Para_(EDMFile));
     }

     blockline = (FirstY - 1) / Para_(BlockSize) + 1;
     blockpos = (FirstX - 1) / Para_(BlockSize) + 1;
     for(y=Image->Bound.Upper; y<=Image->Bound.Lower; y+=Para_(BlockSize), blockline++)
     {
        DecodePixels(Para_(EDMFile), Image, Para_(LayerImage), y, blockline, blockpos, Para_(BlockSize));
     }

     DoneQM();
     Para_(Position) = EDMPIXELDATA;

  }
  /* else cluster is white (white by default - no problem) */

  return(retstatus);

}


/*-------------------------------------------------------------------*/


YESNO ReadEDMPreviewData(EDMPARAMETERS* EDMpointer)
/* return NO if Preview doesn't possible */

{
  PEDMParameters = EDMpointer;

  /* Inittialize EDM variables */
  if (Para_(BlockCode) == NULL) InitBlockCodes();

  /* Read the preview data */
  if (PreviewIsPossible(EDMpointer))
  {
     if (Para_(Position)==EDMTEXTHEADER)
        ReadBlockCodes(Para_(EDMFile));
     else return(YES); /* suppose that data is already READ succesfully */

  } else
  {
     /* ERROR_(NOPREVIEW); */
     EDMErrorStatus_=NOPREVIEW;
     return(NO);
  }

  return(YES);

}


/*-------------------------------------------------------------------*/


YESNO WriteEDMPreviewImage (EDMPARAMETERS* EDMpointer, IMAGE* Image)
/* store Preview Image into the Image */
{
  PEDMParameters = EDMpointer;
  return(WritePreview(Image));
}


/*-------------------------------------------------------------------*/


YESNO ReadEDMImageData(EDMPARAMETERS* EDMpointer)
/* Return NO if EDM image isn't clustered or extra model file or layer */
/* file is needed, or any other errors.                                */
{
  PEDMParameters = EDMpointer;

  if (!(ClusteringIsPossible(EDMpointer))) return(NO);

  /* Inittialize EDM variables */
  if (Para_(BlockCode) == NULL) InitBlockCodes();
  if (Para_(StateIndex) == NULL) InitEDMModel();
  if (Para_(ClusterIndex) == NULL) InitClusterIndex();

  /* Initialize Model File */
  if (Para_(ModellingType)==STATIC) {ERROR_(UNSUPMODEL); return (NO);}
  /* YET NOTHING - CURRENTLY UNSUPPORTED */

  /* Initialize Model File */
  if (Para_(ContextTreeMode)!=NO_CT) {ERROR_(UNSUPCT); return (NO);}
  /* YET NOTHING - CURRENTLY UNSUPPORTED */

  /* !!! */

  /* Read the rest of header */
  if (ReadHeader(Para_(EDMFile),NULL,NULL)==NO) return(NO);

  /* Initialize Layer Image */
  if (Para_(LayerContextSize)>0) {ERROR_(UNSUPLAYER); return (NO); }
  /* YET NOTHING - CURRENTLY UNSUPPORTED */

   return(YES);
}


/*-------------------------------------------------------------------*/
/*-------------------------------------------------------------------*/
/*-------------------------------------------------------------------*/

  /* CODE which is RESERVED for FUTURE */
  /*

      if(Para_(LayerContextSize)>0)
      {
         if (Para_(Clustering) > 1)
         {
            BSize=Para_(ClusterSize) + 2*CONTEXTHEIGHT;
         }
         if ( (Para_(LayerImage)=(IMAGE*) allocate( sizeof(IMAGE) ) ) == NULL ) {ERROR_(EDM_NOMEMORY); exit(-1); }
         ImageInit(LayerName, Para_(LayerImage), INPUT, NO, AUTOMATIC, EBB, BSize, DEFAULTCOLOR);
      }

      if(Para_(ModellingType)==STATIC)
      {
         if (ModelName[0]==0)
         {
            ERROR_(NOMODEL);
            return(NO);
         }
         if (((Para_(ModelFile)=FileOpen(ModelName, INPUT, NO))==NULL)
            return (NO);
         cs=ReadQmtHeader(Para_(ModelFile));
         if (cs!=Para_(ContextSize))
         {
            ERROR_(QMTCONTEXT);
            return(NO);
         }
      }

  */


