/*-------------------------------------------------------------------*/
/* QMCODER.C       Pavel Kopylov                                     */
/*                                                                   */
/* QM-coder module;                                                  */
/* modified from the source code from E.Ageenko & P.Franti           */
/* modified from the source code from AT&T.                          */
/*                                                                   */
/*                                                                   */
/* Supports multiple instances                                       */
/*-------------------------------------------------------------------*/

#define  ModuleName     "QMCODER"
#define  VersionNumber  "Version 0.02"
#define  LastUpdated    "27.03.2000"

// -----------------------------------------------------------------

#include "qmcoder.h"

// ----------------------------------------------------------------- 

static const double prob[128] = {
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

static const unsigned int lsz[128] = {
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

static const int swtch[128] = {
     1,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   1,
     0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
     0,   0,   0,   0,   0,   0,   1,   0,   0,   0,   0,   0,   0,   0,   0,
     0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
     0,   0,   0,   0,   1,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
     0,   0,   0,   0,   0,   1,   0,   0,   0,   0,   0,   0,   0,   1,   0,
     0,   0,   0,   0,   0,   1,   0,   0,   0,   0,   0,   0,   0,   0,   0,
     1,   0,   0,   0,   0,   1,   0,   1,   0,   0,   0,   0,   0,   0,   0,
     0,   0,   0,   0,   0,   0,   0,   0};

static const int nmps[128] = {
     1,   2,   3,   4,   5,   6,   7,   8,   9,  10,  11,  12,  13,  13,  15,
    16,  17,  18,  19,  20,  21,  22,  23,  24,  25,  26,  27,  28,  29,  30,
    31,  32,  33,  34,  35,   9,  37,  38,  39,  40,  41,  42,  43,  44,  45,
    46,  47,  48,  49,  50,  51,  52,  53,  54,  55,  56,  57,  58,  59,  60,
    61,  62,  63,  32,  65,  66,  67,  68,  69,  70,  71,  72,  73,  74,  75,
    76,  77,  78,  79,  48,  81,  82,  83,  84,  85,  86,  87,  71,  89,  90,
    91,  92,  93,  94,  86,  96,  97,  98,  99, 100,  93, 102, 103, 104,  99,
   106, 107, 103, 109, 107, 111, 109, 111,   0,   0,   0,   0,   0,   0,   0,
     0,   0,   0,   0,   0,   0,   0,   0};

static const int nlps[128] = {
     1,  14,  16,  18,  20,  23,  25,  28,  30,  33,  35,   9,  10,  12,  15,
    36,  38,  39,  40,  42,  43,  45,  46,  48,  49,  51,  52,  54,  56,  57,
    59,  60,  62,  63,  32,  33,  37,  64,  65,  67,  68,  69,  70,  72,  73,
    74,  75,  77,  78,  79,  48,  50,  50,  51,  52,  53,  54,  55,  56,  57,
    58,  59,  61,  61,  65,  80,  81,  82,  83,  84,  86,  87,  87,  72,  72,
    74,  74,  75,  77,  77,  80,  88,  89,  90,  91,  92,  93,  86,  88,  95,
    96,  97,  99,  99,  93,  95, 101, 102, 103, 104,  99, 105, 106, 107, 103,
   105, 108, 109, 110, 111, 110, 112, 112,   0,   0,   0,   0,   0,   0,   0,
     0,   0,   0,   0,   0,   0,   0,   0};

// ------------------- Semi-adaptive modelling --------------------- 

static const int IndexForOrderedProbability[113] = {
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

static const double OrderedProbability[113] = {
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

// ------------------- Forward-adaptive modelling ------------------

static const double FirstAttackStateBounds[13] = {
  .30891, .14590, .06891, .03255, .01537, .00726, .00343, .00162, .00076,
  .00036, .00017, .00008, .00004};

// -----------------------------------------------------------------

static void ByteIn(QMDEC *);
void SendByte(QMENC *, BYTE, int);

#define EUP  {*pI = nmps[*pI];}
#define SUP  {*pM = swtch[*pI] ? 1-*pM : *pM; *pI = nlps[*pI];}

#define ESC         0xff        // escape       
#define STUFF       0x00        // escape escape
#define SDNORM      0x02        // Normal end


#define RENORME(QM)  while(QM->A<0x8000) {QM->A <<= 1; QM->C <<= 1; \
                     QM->CT--; if(QM->CT==0) ByteOut(QM); }

#define RENORMD(QM)  do { if(QM->CT==0) ByteIn(QM); \
                       QM->A <<= 1; QM->C<<=1; QM->CT--;} \
                     while(QM->A<0x8000);

#define OUT(QM, A)   if(QM->stflag) QM->stflag = 0; \
                     else { \
                       if((A)==0) QM->nzero++; \
                       else { \
                         while(QM->nzero>0) { \
                           QM->nzero--; \
                           SendByte(QM, (BYTE)0x00, 0); } \
                         SendByte(QM, (BYTE)A, 0); } }

#define OUT00(QM)   if(QM->stflag) QM->stflag = 0; else QM->nzero++;

#define OUTFF(QM)   if(QM->stflag) QM->stflag = 0; \
                    else { \
                      while(QM->nzero>0) { QM->nzero--; \
                      SendByte(QM, (BYTE)0x00, 0);} \
                      SendByte(QM, (BYTE)0xFF, 0); }



/*--------------------------------------------------------------------*/

double GetProbabilityQM(int bit, BYTE *pI, BYTE *pM)
{
  //printf("Index=%3i  lsz=0x%04x  mps=%i \n", *pI, lsz[*pI], *pM);
  if(bit== *pM) return(1 - prob[*pI]);
  else          return(prob[*pI]);
}




/*========================= Initialization ===========================*/

void QMSetInFunction(QMDEC *pQMDec, BYTE (*f)(void *),void *pDec )
{
  pQMDec->fByteIn = f;
  pQMDec->pDecoder = pDec;
}

void QMSetOutFunction(QMENC *pQMEnc, void (*f)(void *, BYTE),void *pEnc)
{
  pQMEnc->fByteOut = f;
  pQMEnc->pEncoder = pEnc;
}


/*=========================== I/O routines ===========================*/

static void ByteOut(QMENC *pQMEnc)
{
  ULONG temp;

  temp = (ULONG)(pQMEnc->C>>19)&0x1ff;
  if(temp>0xff)
  {
    OUT(pQMEnc, pQMEnc->buffer+1)
    while(pQMEnc->sc>0)
    { 
      pQMEnc->sc --;
      OUT00(pQMEnc)
    } 
    pQMEnc->buffer = (int)temp&0xff;
  }
  else
  {
    if(temp==0xff) pQMEnc->sc ++;
    else
    {
      OUT(pQMEnc, pQMEnc->buffer)
      while(pQMEnc->sc > 0)
      { 
        pQMEnc->sc--;
        OUTFF(pQMEnc)
      } 
      pQMEnc->buffer = (int)temp; 
    }
  }
  pQMEnc->C &= 0x7ffffL;
  pQMEnc->CT = 8;
}

/*--------------------------------------------------------------------*/

void SendByte(QMENC *pQMEnc, BYTE x, int mode)
{
  pQMEnc->lNumBytesOut ++;
  pQMEnc->fByteOut(pQMEnc->pEncoder, x);
  
  
  if(x==ESC)
  {
    if (mode != 1)
    {
      pQMEnc->lNumBytesOut ++;
      pQMEnc->fByteOut(pQMEnc->pEncoder, STUFF);
    }
    else
    {
      pQMEnc->lNumBytesOut ++;
      pQMEnc->fByteOut(pQMEnc->pEncoder, SDNORM);
    }
  }
}

/*-------------------------------------------------------------------*/
     
void QMInitEnc(QMENC *pQMEnc)
{
  pQMEnc->C=0;               /* Coding register */
  pQMEnc->buffer=0x00;       /* Bits buffered for output or input */
  pQMEnc->stflag=1;          /* Flag to inhibit first write */
  pQMEnc->nzero=0;           /* Potential trailing zeros */
  pQMEnc->sc=0;              /* Number of FF bytes */
  pQMEnc->CT=11;             /* Shift counter */
  pQMEnc->lNumBytesOut=0;

  pQMEnc->A = 0;
  //pQMEnc->A = 0x10000l;

  /* Note: it is assumed Short = 16 bits. In case of 32 bit use this: */
  /* A = 0x10000L; */
}

/*=========================== QM-encoder =============================*/


void QMEncode(QMENC *pQMEnc, BOOLN Bit, BYTE *pI, BYTE *pM)
{
  //printf("1: Index=%3i,  mps=%i  Bit=%i  a=%04x  c=%08lx \n",*pI,*pM,Bit,pQMEnc->A,pQMEnc->C);

  pQMEnc->A -= lsz[*pI];
  if(Bit==*pM)
  {
    if(pQMEnc->A < 0x8000U)
    { 
      if(pQMEnc->A < lsz[*pI])
      { 
        pQMEnc->C += pQMEnc->A; 
        pQMEnc->A = lsz[*pI];
      }
      EUP; 
      RENORME(pQMEnc)
    }
  }
  else
  {
    if(pQMEnc->A >= lsz[*pI])
    { 
      pQMEnc->C += pQMEnc->A;
      pQMEnc->A = lsz[*pI];
    }
    SUP;
    RENORME(pQMEnc)
  }

  //printf("2: Index=%3i,  mps=%i  Bit=%i  a=%04x  c=%08lx \n\n",*pI,*pM,Bit,pQMEnc->A,pQMEnc->C);
}

/*--------------------------------------------------------------------*/

ULONG QMFlush(QMENC *pQMEnc)
{
  ULONG temp;

  /* printf("\nFLUSH encoding \n\n"); */

  temp = (pQMEnc->C + pQMEnc->A - 1) & 0xffff0000L;
  if(temp < pQMEnc->C) pQMEnc->C = temp + 0x8000L;
  else pQMEnc->C = temp;
  pQMEnc->C <<= pQMEnc->CT;
  if(pQMEnc->C > 0x7ffffffL)
  { 
    OUT(pQMEnc, pQMEnc->buffer+1)
    while(pQMEnc->sc>0)
    {
      pQMEnc->sc--;
      OUT00(pQMEnc)
    }
  }
  else 
  { 
    OUT(pQMEnc, pQMEnc->buffer) 
    while(pQMEnc->sc>0) 
    { 
      pQMEnc->sc--; 
      OUTFF(pQMEnc)
    } 
  }
  OUT(pQMEnc, (int)((pQMEnc->C>>19)&0xff))
  OUT(pQMEnc, (int)((pQMEnc->C>>11)&0xff))

  SendByte(pQMEnc, ESC, 1);        // End-of-Code-Sequence

  temp = pQMEnc->lNumBytesOut;
  pQMEnc->lNumBytesOut = 0;
  return(temp);
}
       

/*--------------------------------------------------------------------*/

void QMInitDec(QMDEC *pQMDec)
{
  pQMDec->buffer=0x00;       /* Bits buffered for output or input */
  pQMDec->stflag=0;          /* Flag to inhibit first write */
  pQMDec->nzero=0;           /* Potential trailing zeros */
  pQMDec->sc=0;              /* Number of FF bytes */
  pQMDec->CT=0;              /* Shift counter */
  pQMDec->lNumBytesIn=0;
  pQMDec->pacfeed = 0;

  pQMDec->C = 0;
  ByteIn(pQMDec);
  pQMDec->C <<= 8;
  ByteIn(pQMDec);
  pQMDec->C <<= 8;
  ByteIn(pQMDec);
  pQMDec->A = 0;

  /* !!! This commennt is incorrect !!! */
  /* Note: it is assumed Short = 16 bits. In case of 32 bit use this: */
  /* A = 0x10000L; */
}

/*=========================== QM-decoder =============================*/

BOOLN QMDecode(QMDEC *pQMDec, BYTE *pI, BYTE *pM)
{
  int Bit;

  pQMDec->A -= lsz[*pI];
  if((pQMDec->C >> 16) < pQMDec->A)
  {
    if(pQMDec->A < 0x8000) 
    {
      if(pQMDec->A < lsz[*pI]) { Bit = 1-*pM; SUP; } else { Bit = *pM; EUP; }
      RENORMD(pQMDec) 
    }
    else Bit = *pM;
  }
  else
  {
    if(pQMDec->A < lsz[*pI]) 
    { 
      Bit = *pM;
      pQMDec->C -= (ULONG) pQMDec->A << 16;
      pQMDec->A = lsz[*pI];
      EUP;
      RENORMD(pQMDec) 
    }
    else
    { 
      Bit = 1-*pM;
      pQMDec->C -= (ULONG) pQMDec->A << 16;
      pQMDec->A = lsz[*pI];
      SUP;
      RENORMD(pQMDec) 
    }
  }

  return Bit;
}

/*-------------------------------------------------------------------*/

BYTE InputByte(QMDEC *pQMDec)
{
  int x;


  if (pQMDec->pacfeed == 1 ) return(0);

  x = pQMDec->fByteIn(pQMDec->pDecoder);
  if(x != ESC) 
  {
    pQMDec->lNumBytesIn++;
    return (x);
  }
 
  x = pQMDec->fByteIn(pQMDec->pDecoder);
  switch(x)
  {
    case STUFF:  pQMDec->lNumBytesIn++;
                 return(ESC);
    case SDNORM: pQMDec->pacfeed = 1;
                 pQMDec->lNumBytesIn++;
                 return(0);
    case EOF:    ErrorMessage("QM: EOF detected. \n");
                 exit(-1);
    default:     ErrorMessage("QM: Unsupported ESC-code detected at %d\n", ftell(pQMDec->pDecoder));
                 exit(-1);
  }
  return(0);
}


void ByteIn(QMDEC *pQMDec)
{
  ULONG temp;

  temp = (ULONG) InputByte(pQMDec) << 8;
  pQMDec->C += temp;
  pQMDec->CT = 8;
}


/*==================== Semi-adaptive QM-coder ========================*/

int FindNearestState(double Prob)
{
  double  Halfpoint;
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

int GetStateIndex(double WhiteProb)
{
  double  LpsProb;
  int     index;

  LpsProb = WhiteProb < 0.5 ? WhiteProb : 1-WhiteProb;
  index = WhiteProb < 0.5 ? 0x80 : 0x00;
  index = index | FindNearestState(LpsProb);
  return(index);
}


/*-------------------------------------------------------------------*/

void RestoreStateIndex(BYTE *pI, BYTE *pM, int value)
{
  *pM = (value & 0x80) ? 1 : 0;
  *pI = (value & 0x7f);
}



/*==================== Forward-adaptive QM-coder =====================*/

int FindNearestFirstAttackState(double Prob)
{
  int    i;

  for(i=0; i<14; i++)
  {
    if( Prob>FirstAttackStateBounds[i] )
    {
      return(i);
    }
  }
  return(13);
}

/*-------------------------------------------------------------------*/

int GetFirstAttackStateIndex(double WhiteProb)
{
  double  LpsProb;
  int     index;

  LpsProb   = WhiteProb < 0.5 ? WhiteProb : 1-WhiteProb;
  index     = WhiteProb < 0.5 ? 0x10 : 0x00;
  index     = index | FindNearestFirstAttackState(LpsProb);
  return (index);
}


/*-------------------------------------------------------------------*/

void RestoreFirstAttackStateIndex(BYTE *pI, BYTE *pM, int index)
{
  *pM = (index & 0x10) ? 1 : 0;
  *pI = (index & 0x0f);
}



/*======================= Modelling part ============================*/

void InitModel(QMMODEL* model, int MaxStates)
{
  int i;

  model->Index = (BYTE*)allocate(sizeof(BYTE)*MaxStates);
  model->MPS   = (BYTE*)allocate(sizeof(BYTE)*MaxStates);

  for(i=MaxStates-1; i>=0; i--)
  {
    model->MPS[i] = 0;
    model->Index[i] = 0;
  }
}


void ReInitModel(QMMODEL* model, int MaxStates)
{
  int i;

  for(i=MaxStates-1; i>=0; i--)
  {
    model->MPS[i] = 0;
    model->Index[i] = 0;
  }
}


void DoneModel(QMMODEL* model)
{
  free(model->MPS);
  free(model->Index);
}



/*=========================== User part =============================*/

void OutputByteToFile(FILE* f, BYTE x)
{   
  putc(x,f);
}

BYTE InputByteFromFile(FILE* f)
{
  BYTE x;

  x = (BYTE)getc(f);
  return(x);
}

void QMInitEncoder(QMENC *pQMEnc, FILE *f)
{
  QMSetOutFunction(pQMEnc, OutputByteToFile, f);
  QMInitEnc(pQMEnc);
}

void QMInitDecoder(QMDEC *pQMDec, FILE *f)
{
  QMSetInFunction(pQMDec, InputByteFromFile, f);
  QMInitDec(pQMDec);
}
