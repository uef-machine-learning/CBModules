#ifndef QMCODER_H
#define QMCODER_H

#include "file.h"
#include "owntypes.h"
#include <stdlib.h>
#include "interfc.h"
#include "memctrl.h"

typedef unsigned long ULONG;
typedef signed   long SLONG;
typedef unsigned long BOOLN;
typedef unsigned char BYTE;


typedef struct tagQMMODEL
{
  BYTE *Index;
  BYTE *MPS;
} QMMODEL;


typedef struct tagQMDEC
{
  ULONG C;
  UINT16 A;
  ULONG CT;             /* Shift counter */
  ULONG lNumBytesIn;

  int  buffer;          /* Bits buffered for output or input */
  int  stflag;          /* Flag to inhibit first write */
  int  nzero;           /* Potential trailing zeros */
  int  sc;              /* Number of FF bytes */
  int  pacfeed;         /* End Of Code flag */

  void *pDecoder;
  BYTE (*fByteIn)(void *);
} QMDEC;

typedef struct tagQMENC
{
  ULONG C;
  UINT16 A;
  ULONG CT;             /* Shift counter */
  ULONG lNumBytesOut;

  int  buffer;          /* Bits buffered for output or input */
  int  stflag;          /* Flag to inhibit first write */
  int  nzero;           /* Potential trailing zeros */
  int  sc;              /* Number of FF bytes */
  int  pacfeed;         /* End Of Code flag */

  void *pEncoder;
  void (*fByteOut)(void *, BYTE);
} QMENC;


void QMInitEncoder(QMENC *, FILE *);
void QMInitEnc(QMENC *);
void QMEncode(QMENC *, BOOLN, BYTE *, BYTE *);
void QMSetOutFunction(QMENC *, void (*f)(void *, BYTE), void *);
void QMEncodeStats(QMENC *, BOOLN, BYTE *, BYTE *, double *);
ULONG QMFlush(QMENC *);

void QMInitDecoder(QMDEC *, FILE *);
void QMInitDec(QMDEC *);
void QMSetInFunction(QMDEC *, BYTE (*f)(void *), void *);
BOOLN QMDecode(QMDEC *, BYTE *, BYTE *);

/*----------------- Interface for forward-adative QM -----------------*/

int GetFirstAttackStateIndex(double);
void RestoreFirstAttackStateIndex(BYTE *, BYTE*, int);

/*------------------- Interface for Modeling -------------------------*/

void InitModel(QMMODEL *, int);
void ReInitModel(QMMODEL *, int);
void DoneModel(QMMODEL *);

#endif
