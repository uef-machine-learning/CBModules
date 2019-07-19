/*-------------------------------------------------------------------*/
/* QM.C            Pasi Fr„nti                                       */
/*                 Eugene Ageenko                                    */
/*                                                                   */
/* QM-coder module; modified from the source code from AT&T.         */
/*    - Semi-adaptive model supported.                               */
/*    - W/o escape sequences mode supported                          */
/*      (EscMode=0; EndCodePos - position in the code file,          */
/*       right after QM-code ends)                                   */
/*                                                                   */
/* 32-bit version                                                    */
/*                                                                   */
/*-------------------------------------------------------------------*/

#define  ModuleName     "QM"
#define  VersionNumber  "Version 0.13a"
#define  LastUpdated    "28.7.99"

/* ----------------------------------------------------------------- */

#include <stdio.h>
#include <stdlib.h>
#include "qm.h"
#include "interfc.h"

/* ----------------------------------------------------------------- */

long   bytesinqm=0;      /* Number of bytes read  at decoding stage */
long   bytesoutqm=0;     /* Number of bytes write at encoding stage */
int    EscMode=1;
long   EndCodePos=0;

/* ----------------------------------------------------------------- */

typedef  unsigned short    U16;
typedef  unsigned long     U32;
typedef  unsigned char     BYTE;

static float prob[128] = {
  .49690, .20691, .09417, .04435, .02120, .01021, .00493, .00239, .00116,
  .00056, .00028, .00013, .00006, .00002, .49901, .34819, .24784, .17912,
  .13081, .09654, .07132, .05310, .03961, .02955, .02219, .01661, .01241,
  .00933, .00698, .00528, .00394, .00297, .00224, .00168, .00127, .00095,
  .50112, .39866, .32010, .25884, .21021, .17204, .14147, .11631, .09630,
  .07970, .06606, .05497, .04620, .03873, .03199, .02684, .02238, .01867,
  .01559, .01301, .01086, .00905, .00758, .00631, .00530, .00437, .00368,
  .00308, .50218, .42468, .35937, .30793, .26416, .22737, .19560, .17023,
  .14701, .12851, .11106, .09710, .08502, .07343, .06458, .05652, .48632,
  .42519, .37251, .33010, .29186, .25740, .22940, .20450, .47112, .42272,
  .37964, .34261, .30957, .27958, .25415, .47784, .43713, .39644, .36288,
  .33216, .30530, .45322, .41940, .38722, .36044, .47506, .44611, .41643,
  .47196, .44283, .49662, .46944, .49582, .00000, .00000, .00000, .00000,
  .00000, .00000, .00000, .00000, .00000, .00000, .00000, .00000, .00000,
  .00000, .00000};

static U16 lsz[128] = {
  0x5a1d,0x2586,0x1114,0x080b,0x03d8,0x01da,0x00e5,0x006f,0x0036,
  0x001a,0x000d,0x0006,0x0003,0x0001,0x5a7f,0x3f25,0x2cf2,0x207c,
  0x17b9,0x1182,0x0cef,0x09a1,0x072f,0x055c,0x0406,0x0303,0x0240,
  0x01b1,0x0144,0x00f5,0x00b7,0x008a,0x0068,0x004e,0x003b,0x002c,
  0x5ae1,0x484c,0x3a0d,0x2ef1,0x261f,0x1f33,0x19a8,0x1518,0x1177,
  0x0e74,0x0bfb,0x09f8,0x0861,0x0706,0x05cd,0x04de,0x040f,0x0363,
  0x02d4,0x025c,0x01f8,0x01a4,0x0160,0x0125,0x00f6,0x00cb,0x00ab,
  0x008f,0x5b12,0x4d04,0x412c,0x37d8,0x2fe8,0x293c,0x2379,0x1edf,
  0x1aa9,0x174e,0x1424,0x119c,0x0f6b,0x0d51,0x0bb6,0x0a40,0x5832,
  0x4d1c,0x438e,0x3bdd,0x34ee,0x2eae,0x299a,0x2516,0x5570,0x4ca9,
  0x44d9,0x3e22,0x3824,0x32b4,0x2e17,0x56a8,0x4f46,0x47e5,0x41cf,
  0x3c3d,0x375e,0x5231,0x4c0f,0x4639,0x415e,0x5627,0x50e7,0x4b85,
  0x5597,0x504f,0x5a10,0x5522,0x59eb,0x0000,0x0000,0x0000,0x0000,
  0x0000,0x0000,0x0000,0x0000,0x0000,0x0000,0x0000,0x0000,0x0000,
  0x0000,0x0000};

static int swtch[128] = {
     1,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   1,
     0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
     0,   0,   0,   0,   0,   0,   1,   0,   0,   0,   0,   0,   0,   0,   0,
     0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
     0,   0,   0,   0,   1,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
     0,   0,   0,   0,   0,   1,   0,   0,   0,   0,   0,   0,   0,   1,   0,
     0,   0,   0,   0,   0,   1,   0,   0,   0,   0,   0,   0,   0,   0,   0,
     1,   0,   0,   0,   0,   1,   0,   1,   0,   0,   0,   0,   0,   0,   0,
     0,   0,   0,   0,   0,   0,   0,   0};

static int nmps[128] = {
     1,   2,   3,   4,   5,   6,   7,   8,   9,  10,  11,  12,  13,  13,  15,
    16,  17,  18,  19,  20,  21,  22,  23,  24,  25,  26,  27,  28,  29,  30,
    31,  32,  33,  34,  35,   9,  37,  38,  39,  40,  41,  42,  43,  44,  45,
    46,  47,  48,  49,  50,  51,  52,  53,  54,  55,  56,  57,  58,  59,  60,
    61,  62,  63,  32,  65,  66,  67,  68,  69,  70,  71,  72,  73,  74,  75,
    76,  77,  78,  79,  48,  81,  82,  83,  84,  85,  86,  87,  71,  89,  90,
    91,  92,  93,  94,  86,  96,  97,  98,  99, 100,  93, 102, 103, 104,  99,
   106, 107, 103, 109, 107, 111, 109, 111,   0,   0,   0,   0,   0,   0,   0,
     0,   0,   0,   0,   0,   0,   0,   0};

static int nlps[128] = {
     1,  14,  16,  18,  20,  23,  25,  28,  30,  33,  35,   9,  10,  12,  15,
    36,  38,  39,  40,  42,  43,  45,  46,  48,  49,  51,  52,  54,  56,  57,
    59,  60,  62,  63,  32,  33,  37,  64,  65,  67,  68,  69,  70,  72,  73,
    74,  75,  77,  78,  79,  48,  50,  50,  51,  52,  53,  54,  55,  56,  57,
    58,  59,  61,  61,  65,  80,  81,  82,  83,  84,  86,  87,  87,  72,  72,
    74,  74,  75,  77,  77,  80,  88,  89,  90,  91,  92,  93,  86,  88,  95,
    96,  97,  99,  99,  93,  95, 101, 102, 103, 104,  99, 105, 106, 107, 103,
   105, 108, 109, 110, 111, 110, 112, 112,   0,   0,   0,   0,   0,   0,   0,
     0,   0,   0,   0,   0,   0,   0,   0};

/* ------------------- Semi-adaptive modelling --------------------- */

static int IndexForOrderedProbability[113] = {
  13,  12,  11,  10,   9,  35,   8,  34,  33,  32,
   7,  31,  63,  62,  30,  61,   6,  29,  60,  59,
  28,  58,  57,  27,   5,  56,  26,  55,  54,  25,
  53,   4,  24,  52,  51,  23,  50,  49,  22,   3,
  48,  21,  47,  79,  78,  46,  20,  77,  45,  76,
   2,  44,  19,  75,  74,  43,  73,  18,  42,  72,
  71,  41,  17,  70,  87,   1,  40,  69,  86,  16,
  94,  85,  39,  68,  93,  84, 100,  67,  92,  38,
  83,  99,  91,  15,  66, 104,  98,  82,  90, 103,
  97,  37, 107, 102,  89,  65,  81,  96, 109, 106,
 101, 111,  88, 108, 105,  95,  80, 112, 110,   0,
  14,  36,  64};

static float OrderedProbability[113] = {
  .00002,.00006,.00013,.00028,.00056,.00095,.00116,.00127,.00168,.00224,
  .00239,.00297,.00308,.00368,.00394,.00437,.00493,.00528,.00530,.00631,
  .00698,.00758,.00905,.00933,.01021,.01086,.01241,.01301,.01559,.01661,
  .01867,.02120,.02219,.02238,.02684,.02955,.03199,.03873,.03961,.04435,
  .04620,.05310,.05497,.05652,.06458,.06606,.07132,.07343,.07970,.08502,
  .09417,.09630,.09654,.09710,.11106,.11631,.12851,.13081,.14147,.14701,
  .17023,.17204,.17912,.19560,.20450,.20691,.21021,.22737,.22940,.24784,
  .25415,.25740,.25884,.26416,.27958,.29186,.30530,.30793,.30957,.32010,
  .33010,.33216,.34261,.34819,.35937,.36044,.36288,.37251,.37964,.38722,
  .39644,.39866,.41643,.41940,.42272,.42468,.42519,.43713,.44283,.44611,
  .45322,.46944,.47112,.47196,.47506,.47784,.48632,.49582,.49662,.49690,
  .49901,.50112,.50218};

/* ------------------- Forward-adaptive modelling ------------------ */

static float FirstAttackStateBounds[13] = {
  .30891, .14590, .06891, .03255, .01537, .00726, .00343, .00162, .00076,
  .00036, .00017, .00008, .00004};

/* ----------------------------------------------------------------- */

static void ByteIn(FILE* f);


#define EUP   { cstate[state]  = nmps[cstate[state]]; }
#define SUP   { mps[state] = swtch[cstate[state]] ? 1-mps[state]: mps[state]; \
                cstate[state]  = nlps[cstate[state]]; }

#define ESC         0xff        /* escape        */
#define STUFF       0x00        /* escape escape */
#define SDNORM      0x02        /* normal stripe data end */

#define RENORME(f)  while(a<0x8000) { \
                          a <<= 1; c <<= 1; ct--; \
                          if(ct==0) (void) ByteOut(f); }
#define RENORMD(f)  do { \
                          if(ct==0) ByteIn(f); \
                          a <<= 1; c <<= 1; ct--; } while(a<0x8000);
#define OUT(f, A)   if(stflag) stflag = 0; else { \
                    if((A)==0) nzero++; \
                    else { \
                       while(nzero>0) {nzero--; OutputByteToFile(f,0x00);} \
                       OutputByteToFile(f,(BYTE)A); } }
#define OUT00           if(stflag) stflag = 0; else nzero++;
#define OUTFF(f)        if(stflag) stflag = 0; else { \
                        while(nzero>0) {nzero--; OutputByteToFile(f,0x00);} \
                        OutputByteToFile(f,0xFF); }


static int* cstate;          /* State tables */
static int* mps;             /* MPS tables */
static U32  c;               /* Coding register */
static U16  a;               /* Size of coding interval*/
static int  stflag;          /* Flag to inhibit first write */
static int  nzero;           /* Potential trailing zeros */
static int  sc;              /* Number of FF bytes */
static int  buffer;          /* Bits buffered for output or input */
static int  pacfeed;         /* End-of-Code-Flag ?? */
static int  ct;              /* Shift counter */


/*==================== Access to gloval variables ====================*/

long   BytesInQM() { return bytesinqm; }
long   BytesOutQM() { return bytesoutqm; }
void   SetQMEscMode(int escmode) {EscMode=escmode;}
void   SetQMEndCodePos(long endcodepos) {EndCodePos=endcodepos;}

/*======================== General routines ==========================*/


float GetProbabilityQM(int bit, int state)
{
/*PrintMessage("state=%3i  cstate=%3i  lsz=0x%04x  mps=%i \n",
                state, cstate[state], lsz[ cstate[state] ], mps[ state ] ); */
if(bit==mps[state])  return( 1 - prob[ cstate[state] ] );
else                 return(     prob[ cstate[state] ] );
}


/*=========================== I/O routines ===========================*/


static void OutputByteToFile(FILE* f, BYTE x)
{
  bytesoutqm++;
  if(f==NULL) return;      /* Used for Pseudo encoding */
  putc(x, f);
  if((EscMode)&&(x==ESC))  /* ESC must be coded: ESC => ESC STUFF */
     {
     bytesoutqm++;
     putc(STUFF, f);
     }
 }


static void ByteOut(FILE* f)
{
  U32 temp;

  temp = (U32) (c>>19)&0x1ff;
  if(temp>0xff) {
    OUT(f, buffer+1) while(sc>0) { sc--; OUT00 } buffer = (int) temp&0xff; }
  else {
    if(temp==0xff) sc++;
    else {
      OUT(f, buffer) while(sc>0) { sc--; OUTFF(f) } buffer = (int) temp; } }
  c &= 0x7ffffL; ct = 8;
}


static BYTE InputByteFromFile(FILE* f)
{
  int x;

  if(pacfeed) return(0);

  if (EscMode)
  {
     x = getc(f);
     if(x != ESC) { bytesinqm++; return(x); }
     switch(getc(f))
     {
        case STUFF:   bytesinqm++; return(ESC); /* ESC STUFF => ESC */
        case SDNORM:  pacfeed = 1; return(0);
        case EOF:     ErrorMessage("ERROR: EOF detected. \n"); ExitProcessing(-1);
        default:      ErrorMessage("ERROR: Unsupported ESC-code detected\n"); ExitProcessing(-1);
     }
     return(0);
  }
  else
  {
     if(EndCodePos==ftell(f))
       { pacfeed = 1;  return(0); /* end-of-code */ }
     else
       { x = getc(f); bytesinqm++; return(x); }
  }
}


static void ByteIn(FILE *f)
{
  U32 temp;

  temp = (U32) InputByteFromFile(f) << 8;
  c += temp;
  ct = 8;
}


/*========================= Initialization ===========================*/


void NewModel(int MaxStates)
{
  int p;

  for(p=0 ; p<MaxStates ; p++)
    {
    mps[p]          = 0;
    cstate[p]       = 0;
    }
}

/*--------------------------------------------------------------------*/

void InitModelQM(int MaxStates)
{
  /* PrintMessage("\nINIT modelling: states=%i \n", MaxStates); */

  cstate = (int*) malloc( MaxStates * sizeof(int) );
  mps    = (int*) malloc( MaxStates * sizeof(int) );
  if( cstate == NULL )
    {
    ErrorMessage("ERROR: cstate[%i] is out of memory", MaxStates);
    ExitProcessing(-1);
    }
  if( mps == NULL )
    {
    ErrorMessage("ERROR: mps[%i] is out of memory", MaxStates);
    ExitProcessing(-1);
    }

  bytesinqm = bytesoutqm = 0;
  NewModel(MaxStates);
}

/*--------------------------------------------------------------------*/

void DoneQM(void)
{
  free(cstate);
  free(mps);
}

/*--------------------------------------------------------------------*/


void InitDecodeQM(FILE* f)
{
   /* static reinitialization for sequential use*/
   c=0;               /* Coding register */
   a=0;               /* Size of coding interval*/
   stflag=0;          /* Flag to inhibit first write */
   nzero=0;           /* Potential trailing zeros */
   sc=0;              /* Number of FF bytes */
   buffer=0;          /* Bits buffered for output or input */
   pacfeed=0;         /* End-of-Code-Flag ?? */
   ct=0;              /* Shift counter */

  /* printf("\nINIT decoding \n"); */

  c = 0; ByteIn(f); c <<= 8; ByteIn(f); c <<= 8; ByteIn(f);
  a = 0;

  /* Note: it is assumed Short = 16 bits. In case of 32 bit use this: */
  /* a = 0x10000L; */
}


/*--------------------------------------------------------------------*/


void InitEncodeQM(void)
{
   /* static reinitialization for sequential use*/
   c=0;               /* Coding register */
   a=0;               /* Size of coding interval*/
   stflag=0;          /* Flag to inhibit first write */
   nzero=0;           /* Potential trailing zeros */
   sc=0;              /* Number of FF bytes */
   buffer=0;          /* Bits buffered for output or input */
   pacfeed=0;         /* End-of-Code-Flag ?? */
   ct=0;              /* Shift counter */

  /* printf("\nINIT encoding \n"); */

  c = 0; a = 0;
  buffer = 0x00; nzero = 0; stflag = 1; sc = 0; ct = 11;

  /* Note: it is assumed Short = 16 bits. In case of 32 bit use this: */
  /* a = 0x10000L; */
}


/*=========================== QM-decoder =============================*/


int DecodeBitByQM(FILE* f, int state)
{
  int Bit;

/* printf("1: state=%3i,  cstate=%3i,  mps=%i  Bit=?  a=%08x  c=%08lx\n",
             state,cstate[state],mps[state],a,c); */

  a -= lsz[cstate[state]];
  if((c>>16) < a) {
    if(a<0x8000) {
      if(a < lsz[cstate[state]]) { Bit = 1-mps[state]; SUP; } else { Bit = mps[state]; EUP; }
      RENORMD(f) }
    else Bit = mps[state]; }
  else {
    if(a<lsz[cstate[state]]) { Bit= mps[state];
      c -= (U32)a<<16; a=lsz[cstate[state]]; EUP; RENORMD(f) }
    else { Bit = 1-mps[state];
      c -= (U32)a<<16; a=lsz[cstate[state]]; SUP; RENORMD(f) } }

/* printf("2: state=%3i,  cstate=%3i,  mps=%i  Bit=%i  a=%08x  c=%08lx\n",
             state,cstate[state],mps[state],Bit,a,c); */

  return Bit;
}


/*=========================== QM-encoder =============================*/


void EncodeBitByQM(FILE*  f,
                   int    state,
                   int    Bit)
{
/* printf("1: state=%3i,  cstate=%3i,  mps=%i  Bit=%i  a=%08x  c=%08lx\n",
              state,cstate[state],mps[state],Bit,a,c); */

  a -= lsz[cstate[state]];
  if(Bit==mps[state])
    {
    if(a<0x8000)
       { if(a<lsz[cstate[state]]) { c += a; a = lsz[cstate[state]]; }
         EUP; RENORME(f)
       }
    }
  else
    {
    if(a>=lsz[cstate[state]]) { c += a; a = lsz[cstate[state]]; }
    SUP;
    RENORME(f)
    }

/*printf("2: state=%3i,  cstate=%3i,  mps=%i  Bit=%i  a=%04x  c=%08lx \n\n",
               state,cstate[state],mps[state],Bit,a,c); */
}


/*--------------------------------------------------------------------*/


void FlushEncodeQM(FILE* f)
{
  U32 temp;

  /* printf("\nFLUSH encoding \n\n"); */

  temp = (c+a-1)&0xffff0000L;
  if(temp<c) c = temp+0x8000L; else c = temp;
  c <<= ct;
  if(c>0x7ffffffL) { OUT(f, buffer+1) while(sc>0) { sc--; OUT00 } }
  else { OUT(f, buffer) while(sc>0) { sc--; OUTFF(f) } }
  OUT(f, (int)((c>>19)&0xff))
  OUT(f, (int)((c>>11)&0xff))

  /* End-of-Code-Symbol */
  if (EscMode)
  {
    putc(ESC, f);
    putc(SDNORM, f);
  }
}


/*==================== Semi-adaptive QM-coder ========================*/


int FindNearestState(float Prob)
{
  float  Halfpoint;
  int    i;

  for(i=1; i<112; i++)
    {
    if( Prob<OrderedProbability[i] )
       {
       Halfpoint = (OrderedProbability[i-1] + OrderedProbability[i]) / 2;
       return( IndexForOrderedProbability[ Prob<Halfpoint ? i-1 : i ] );
       }
    }
  return(0);
}


/*-------------------------------------------------------------------*/


int GetStateIndex(float WhiteProb)
{
  float  LpsProb;
  int    index;

  LpsProb = WhiteProb < 0.5 ? WhiteProb : 1-WhiteProb;
  index = WhiteProb < 0.5 ? 0x80 : 0x00;
  index = index | FindNearestState(LpsProb);
  return(index);
}


/*-------------------------------------------------------------------*/


void RestoreStateIndex(int state, int value)
{
  mps[state]    = (value & 0x80) ? 1 : 0;
  cstate[state] = (value & 0x7f);
}


/*==================== Forward-adaptive QM-coder =====================*/


int FindNearestFirstAttackState(float Prob)
{
  int i;

  for(i=0; i<14; i++)
    {
	if( Prob>FirstAttackStateBounds[i] )  return(i);
    }
  return(13);
}

/*-------------------------------------------------------------------*/

int GetFirstAttackStateIndex(float WhiteProb)
{
  float  LpsProb;
  int    index;

  LpsProb   = WhiteProb < 0.5 ? WhiteProb : 1-WhiteProb;
  index     = WhiteProb < 0.5 ? 0x10 : 0x00;
  index     = index | FindNearestFirstAttackState(LpsProb);
  return (index);
}


/*-------------------------------------------------------------------*/


void RestoreFirstAttackStateIndex(int context, int index)
{
  mps[context]    = (index & 0x10) ? 1 : 0;
  cstate[context] = (index & 0x0f);
}

