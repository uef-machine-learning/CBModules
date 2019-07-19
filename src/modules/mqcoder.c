/*-------------------------------------------------------------------*/
/* MQCODER.C       Pavel Kopylov                                     */
/*                 Eugene Ageenko                                    */
/*                                                                   */
/* MQ-coder module                                                   */
/*                                                                   */
/*                                                                   */
/* Supports multiple instances                                       */
/*-------------------------------------------------------------------*/

#define  ModuleName     "MQCODER"
#define  VersionNumber  "Version 0.05"
#define  LastUpdated    "06.02.2001"

#include "mqcoder.h"

static void MQByteOut(MQENC *);
static void MQSendByte(MQENC *);
static void MQSetBits(MQENC *);
static void MQByteIn(MQDEC *);


/* The MQ Coder routines use a pointer to Decode() & Encode() functions
to send & receive bits.  These "IN" & "OUT" functions are defined externally,
and also include a pointer to allow for decoder & encoder structures

For example, the simplist support functions could be:

void MyOut (void *Ignore, BYTE a)
	{*encpointer++ = a}

BYTE MyIn (void *Ignore)
	{return *decpointer++}

and would be initialized by

MQDEC MyMQDecoder;
MQENC MyMQEncoder;

MQSetOutFunction(&MyMQEncoder,MyOut,NULL);
MQSetInFunction(&MyMQDecoder,MyIn,NULL);

*/

void MQSetInFunction(MQDEC *pMQDec, BYTE (*f)(void *),void *pDec )
{
	pMQDec->fByteIn = f;
	pMQDec->pDecoder = pDec;
}

void MQSetOutFunction(MQENC *pMQEnc, void (*f)(void *, BYTE),void *pEnc)
{
	pMQEnc->fByteOut = f;
	pMQEnc->pEncoder = pEnc;
}

/* Lookup Data tables for the MQ Coder */

static const ULONG Qe[47] = {
  0x5601,0x3401,0x1801,0x0ac1,0x0521,0x0221,0x5601,0x5401,
  0x4801,0x3801,0x3001,0x2401,0x1c01,0x1601,0x5601,0x5401,
  0x5101,0x4801,0x3801,0x3401,0x3001,0x2801,0x2401,0x2201,
  0x1c01,0x1801,0x1601,0x1401,0x1201,0x1101,0x0ac1,0x09c1,
  0x08a1,0x0521,0x0441,0x02a1,0x0221,0x0141,0x0111,0x0085,
  0x0049,0x0025,0x0015,0x0009,0x0005,0x0001,0x5601 };

static const BOOLN SWTCH[47] = {
  1,0,0,0,0,0,1,0,
  0,0,0,0,0,0,1,0,
  0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,  
  0,0,0,0,0,0,0 };

static const BYTE NMPS[47] = {
   1, 2, 3, 4, 5,38, 7, 8,
   9,10,11,12,13,29,15,16,
  17,18,19,20,21,22,23,24,
  25,26,27,28,29,30,31,32,
  33,34,35,36,37,38,39,40,
  41,42,43,44,45,45,46 };

static const BYTE NLPS[47] = {
   1, 6, 9,12,29,33, 6,14,
  14,14,17,18,20,21,14,14,
  15,16,17,18,19,19,20,21,
  22,23,24,25,26,27,28,29,
  30,31,32,33,34,35,36,37,
  38,39,40,41,42,43,46 };

/* Approximate bit costs MQ Coder */

static const double MPSCost[47] = {
1.01140474, 0.52432363, 0.21868226, 0.09389752, 0.04402218, 0.01810936, 1.01140474, 0.97771884,
0.79061900, 0.57379537, 0.47649224, 0.34183697, 0.25857702, 0.19914109, 1.01140474, 0.97771884,
0.92862173, 0.79061900, 0.57379537, 0.52432363, 0.47649224, 0.38533915, 0.34183697, 0.32056719,
0.25857702, 0.21868226, 0.19914109, 0.17986108, 0.16083532, 0.15141652, 0.09389752, 0.08490442,
0.07485278, 0.04402218, 0.03641638, 0.02239622, 0.01810936, 0.01063861, 0.00904369, 0.00439826,
0.00241276, 0.00122248, 0.00069410, 0.00029723, 0.00016448, 0.00003318, 1.01140474 };

static const double LPSCost[47] = {
0.98868471, 1.71446757, 2.82981854, 3.98822959, 5.05634718, 6.32493203, 0.98868471, 1.02263068,
1.24501333, 1.60755947, 1.82993651, 2.24493126, 2.60745836, 2.95532891, 0.98868471, 1.02263068,
1.07509441, 1.24501333, 1.60755947, 1.71446757, 1.82993651, 2.09294937, 2.24493126, 2.32738457,
2.60745836, 2.82981854, 2.95532891, 3.09280782, 3.24478082, 3.32721805, 3.98822959, 4.12902696,
4.30583119, 5.05634718, 5.32620480, 6.02055116, 6.32493203, 7.08862901, 7.32215894, 8.35981593,
9.22507255, 10.20535041, 11.02167549, 12.24506804, 13.09867856, 15.40800661, 0.98868471 };



/* ------------------- Forward-adaptive modelling ------------------ */

static double FirstAttackStateBounds[13] = {
  .391863705, .207022136, .094141584, .043516659, .019361847, .008828931,
  .004361417, .002255332, .001189679, .000638284, .000314779, .000153245,
  5.12055E-05};


/*======================== MQ - Encoder =============================*/

void MQInitEnc(MQENC *pMQEnc)
{
  /* NOTE: I[] & MPS[] must be initialized elsewhere */
  pMQEnc->A = 0x8000;
  pMQEnc->nFF = 0;
  pMQEnc->n7F = 0;
  pMQEnc->bLastByte = 0;
  pMQEnc->lNumBytesOut = 0;
  pMQEnc->outBuffer = 0;
  pMQEnc->bFirstByte = 1;
  pMQEnc->C = 0;
  pMQEnc->CT = 12;
}

/*-------------------------------------------------------------------*/

/* This routine can be used to determine how much it costs to encode data */

void MQEncodeStats(MQENC *pMQEnc, BOOLN b, BYTE *pI ,BYTE *pM, double *Cost)
{
  if (b == *pM) *Cost += MPSCost[*pI];
  else *Cost += LPSCost[*pI];
  
  MQEncode(pMQEnc,b,pI,pM);
}

/*-------------------------------------------------------------------*/

void MQEncode(MQENC *pMQEnc, BOOLN b, BYTE *pI ,BYTE *pM)
{
  BOOLN bRenormE = 0;
  BYTE MPSVal = *pM;
  BYTE IndexVal = *pI;
  ULONG QeVal = Qe[IndexVal];

  if (b == *pM)
  {
    /* MQCodeMPS */
    pMQEnc->A -= QeVal;
    if ((pMQEnc->A & 0x8000)==0)
    {
      if (pMQEnc->A < QeVal) pMQEnc->A = QeVal;
      else pMQEnc->C += QeVal;
      *pI = NMPS[IndexVal];
      bRenormE = 1;
    } else pMQEnc->C += QeVal;
  } 
  else
  {
    /* MQCodeLPS */
    pMQEnc->A -= QeVal;
    if (pMQEnc->A < QeVal) pMQEnc->C += QeVal;
    else pMQEnc->A = QeVal;
    if (SWTCH[IndexVal]) *pM = (BYTE)(1 - MPSVal);
    *pI = NLPS[IndexVal];
    bRenormE = 1;
  }

  if (bRenormE)
    do
    {
      pMQEnc->A <<= 1;
      pMQEnc->C <<= 1;
      pMQEnc->CT --;
      if (pMQEnc->CT ==0) MQByteOut(pMQEnc);
    } while ((pMQEnc->A & 0x8000)==0);
}

/*-------------------------------------------------------------------*/

static void MQByteOut(MQENC *pMQEnc)
{
  BOOLN bPad = 0;
  if (pMQEnc->outBuffer == 0xFF) bPad = 1;
  else 
  if (pMQEnc->C > 0x7FFFFFF)
  {
    pMQEnc->outBuffer ++;
    if (pMQEnc->outBuffer == 0xFF)
    {
      pMQEnc->C &= 0x7FFFFFF;
      bPad = 1;
    }
  }

  MQSendByte(pMQEnc);

  if (bPad==1)
  {
    pMQEnc->outBuffer = (BYTE) (pMQEnc->C >> 20);
    pMQEnc->C &= 0xFFFFF;
    pMQEnc->CT = 7;
  }
  else
  {
    pMQEnc->outBuffer = (BYTE) (pMQEnc->C >> 19);
    pMQEnc->C &= 0x7FFFF;
    pMQEnc->CT = 8;
  }
}

/*-------------------------------------------------------------------*/

ULONG MQFlushEncoder(MQENC *pMQEnc)
{
  MQSetBits(pMQEnc);
  pMQEnc->C <<= pMQEnc->CT;
  MQByteOut(pMQEnc);

  pMQEnc->C = pMQEnc->C | 0x7FFF;

  pMQEnc->C <<= pMQEnc->CT;
  MQByteOut(pMQEnc);

  pMQEnc->C = pMQEnc->C | 0x7FFF;

  if (pMQEnc->outBuffer != 0xFF)
  {
    pMQEnc->C <<= pMQEnc->CT;
    MQByteOut(pMQEnc);
  }

  pMQEnc->bLastByte = 1;
  MQSendByte(pMQEnc);
	
  pMQEnc->fByteOut(pMQEnc->pEncoder,0xAC);
  pMQEnc->lNumBytesOut++;

  return pMQEnc->lNumBytesOut;
}

/*-------------------------------------------------------------------*/

static void MQSendByte(MQENC *pMQEnc)
{
  if (pMQEnc->bFirstByte==1) pMQEnc->bFirstByte = 0;
  else 
  {  
    /* check for repeating FF7F in output */
    if ( ((pMQEnc->outBuffer==0x7F) && (pMQEnc->n7F == pMQEnc->nFF - 1)) || ((pMQEnc->outBuffer==0xFF) && (pMQEnc->nFF == pMQEnc->n7F)) )
    {
      if (pMQEnc->outBuffer==0x7F) pMQEnc->n7F ++;
      else pMQEnc->nFF ++;

      if (pMQEnc->bLastByte)
      {
        pMQEnc->fByteOut(pMQEnc->pEncoder,0xFF);
        pMQEnc->lNumBytesOut++;
      }
    }
    else
    {
      while (pMQEnc->nFF > 0)
      {
        pMQEnc->fByteOut(pMQEnc->pEncoder,0xFF);
        pMQEnc->lNumBytesOut++;
        pMQEnc->nFF --;
        if (pMQEnc->n7F > 0)
        {
          pMQEnc->fByteOut(pMQEnc->pEncoder,0x7F);
          pMQEnc->lNumBytesOut++;
          pMQEnc->n7F --;
        }
      }
      pMQEnc->fByteOut(pMQEnc->pEncoder,pMQEnc->outBuffer);
      pMQEnc->lNumBytesOut++;
    }
  }
}

/*-------------------------------------------------------------------*/

static void MQSetBits(MQENC *pMQEnc)
{
  ULONG T;

  T = pMQEnc->C + pMQEnc->A;
  pMQEnc->C = pMQEnc->C | 0xFFFF;
  if (pMQEnc->C >= T) pMQEnc->C -= 0x8000;
}


/*======================== MQ - Decoder =============================*/

void MQInitDec(MQDEC *pMQDec)
{
  /* NOTE: I[] & MPS[] must be initialized elsewhere */

  pMQDec->bEOF = 0;
  pMQDec->C = 0;
  pMQDec->inBuffer = pMQDec->fByteIn(pMQDec->pDecoder);
  pMQDec->lNumBytesOut = 1;
  
  pMQDec->C += pMQDec->inBuffer << 16;
  MQByteIn(pMQDec);

  pMQDec->C <<= 7;
  pMQDec->CT -= 7;

  pMQDec->A = 0x8000;
}


BOOLN MQDecode(MQDEC *pMQDec, BYTE *pI, BYTE *pM)
{
  BOOLN bBit = 0;
  BOOLN bRenormD = 0;
  BYTE MPSVal = *pM;
  BYTE IndexVal = *pI;
  ULONG QeVal = Qe[IndexVal];

  pMQDec->A -= QeVal;

  if ((pMQDec->C>>16) >= QeVal)
  {
    pMQDec->C -= (QeVal << 16);
    if ((pMQDec->A & 0x8000)==0)
    {
      /* ExchangeMPS() */
      if (pMQDec->A < QeVal)
      {
        bBit = 1 - MPSVal;
        if (SWTCH[IndexVal]==1) *pM = (BYTE)(1 - MPSVal);
        *pI = NLPS[IndexVal];
      }
      else
      {
        bBit = MPSVal;
        *pI = NMPS[IndexVal];
      }
      /**/
      bRenormD = 1;
    }
    else bBit = MPSVal;
  }
  else
  {
    /* MQExchangeLPS() */
    if (pMQDec->A < QeVal)
    { 
      pMQDec->A = QeVal;
      bBit = MPSVal;
      *pI = NMPS[IndexVal];
    }
    else
    {
      pMQDec->A = QeVal;
      bBit = 1 - MPSVal;
      if (SWTCH[IndexVal]==1) *pM = (BYTE)(1 - MPSVal);
      *pI = NLPS[IndexVal];
    }
    /**/
    bRenormD = 1;
  }
	
  if (bRenormD)
    do
    {
      if (pMQDec->CT ==0) MQByteIn(pMQDec);
      pMQDec->A <<= 1;
      pMQDec->C <<= 1;
      pMQDec->CT --;
    }
    while ((pMQDec->A & 0x8000)==0);

  return bBit;
}

/*-------------------------------------------------------------------*/

static void MQByteIn(MQDEC *pMQDec)
{
  BYTE prevBuffer = pMQDec->inBuffer;

  if (pMQDec->bEOF==0)
  {
    pMQDec->inBuffer = pMQDec->fByteIn(pMQDec->pDecoder);
    pMQDec->lNumBytesOut++;

    if (prevBuffer == 0xFF)
    {
      if (pMQDec->inBuffer > 0x8F) 
        pMQDec->bEOF = 1;
      else
      {
        pMQDec->C += pMQDec->inBuffer << 9;
        pMQDec->CT = 7;
      }
    }
    else
    {
      pMQDec->C += pMQDec->inBuffer << 8;
      pMQDec->CT = 8;
    }
  }

  if (pMQDec->bEOF)
  {
    pMQDec->C += 0xFF00;
    pMQDec->CT = 8;
  }
}

/*-------------------------------------------------------------------*/

ULONG MQFlushDecoder(MQDEC *pMQDec)
{
  BYTE tmpBuffer;

  if (pMQDec->inBuffer == 0xFF && !pMQDec->bEOF)
    do
    {
      tmpBuffer = pMQDec->fByteIn(pMQDec->pDecoder);
      pMQDec->lNumBytesOut++;
    }
    while (tmpBuffer!=0xAC);

  pMQDec->bEOF = 1;
  return pMQDec->lNumBytesOut;
}

/*======================= Modelling part ============================*/

void InitModel(MQMODEL* model, int MaxStates)
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

void ReinitModel(MQMODEL* model, int MaxStates)
{
  int i;

  for(i=MaxStates-1; i>=0; i--)
  {
    model->MPS[i] = 0;
    model->Index[i] = 0;
  }
}

void DoneModel(MQMODEL* model)
{
  free(model->MPS);
  free(model->Index);
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

int GetFirstAttackStateIndex(double WhiteProb)
{
  double  LpsProb;
  int     index;

  LpsProb   = WhiteProb < 0.5 ? WhiteProb : 1-WhiteProb;
  index     = WhiteProb < 0.5 ? 0x10 : 0x00;
  index     = index | FindNearestFirstAttackState(LpsProb);
  return (index);
}

void RestoreFirstAttackStateIndex(BYTE *pI, BYTE *pM, int index)
{
  *pM = (index & 0x10) ? 1 : 0;
  *pI = (index & 0x0f);
}



/*=========================== User part =============================*/

void OutputByteToFile(FILE* f, BYTE x)
{
  if(f==NULL) return;
  putc(x, f);
}

BYTE InputByteFromFile(FILE* f)
{
  BYTE x;

  x = (BYTE)getc(f);
  return(x);
}


void MQInitEncoder(MQENC *pMQEnc, FILE *f)
{
  MQSetOutFunction(pMQEnc, OutputByteToFile, f);
  MQInitEnc(pMQEnc);
}

void MQInitDecoder(MQDEC *pMQDec, FILE *f)
{
  MQSetInFunction(pMQDec, InputByteFromFile, f);
  MQInitDec(pMQDec);
}

