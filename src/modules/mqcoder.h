#ifndef MQCODER_H
#define MQCODER_H

#include "file.h"
#include "owntypes.h"
#include "memctrl.h"

typedef unsigned long ULONG;
typedef signed   long SLONG;
typedef unsigned long BOOLN;
typedef unsigned char BYTE;

typedef struct tagMQMODEL
{
  BYTE *Index;
  BYTE *MPS;
} MQMODEL;

typedef struct tagMQDEC
{
  ULONG C;
  ULONG A;
  ULONG CT;
  BYTE inBuffer;
  BOOLN bEOF;

  ULONG lNumBytesOut;

  void *pDecoder;
  BYTE (*fByteIn)(void *);
} MQDEC;

typedef struct tagMQENC
{
  ULONG C;
  ULONG A;
  ULONG CT;

  BYTE outBuffer;
  BOOLN bFirstByte;
  BOOLN bLastByte;
  ULONG nFF;
  ULONG n7F;

  ULONG lNumBytesOut;

  void *pEncoder;
  void (*fByteOut)(void *, BYTE);
} MQENC;


void MQInitEncoder(MQENC *, FILE *);
void MQInitEnc(MQENC *);
void MQEncode(MQENC *, BOOLN,BYTE *,BYTE *);
void MQSetOutFunction(MQENC *, void (*f)(void *, BYTE), void *);
void MQEncodeStats(MQENC *, BOOLN,BYTE *,BYTE *, double *);
ULONG MQFlushEncoder(MQENC *);


void MQInitDecoder(MQDEC *, FILE *);
void MQInitDec(MQDEC *);
void MQSetInFunction(MQDEC *, BYTE (*f)(void *), void *);
BOOLN MQDecode(MQDEC *, BYTE *, BYTE *);
ULONG MQFlushDecoder(MQDEC *);

/*----------------- Interface for forward-adative QM -----------------*/

int GetFirstAttackStateIndex(double);
void RestoreFirstAttackStateIndex(BYTE *, BYTE*, int);

/*------------------ Interface for Modeling -------------------------*/
void InitModel(MQMODEL *, int);
void ReinitModel(MQMODEL *, int);
void DoneModel(MQMODEL *);

#endif

