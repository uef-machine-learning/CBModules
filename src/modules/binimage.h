#if ! defined(__BINIMAGE_H)
#define __BINIMAGE_H

#include "file.h"
#include "pgm.h"
#include "owntypes.h"

/*---------------------  Image data structure  -----------------------*/

#define  LinesMax         20000
#define  LinesAutoMax     128
#define  LinesMin         10
#define  NameMaxLength    127

#define  OBB EXPANDED
#define  EBB PACKED

#define  TMPNAME "TEMPFILE"

/* typedef  enum  {  WHITE=0, BLACK=1 } BWCOLORTYPE */
 #define  BLACK 1
 #define  WHITE 0

#define  DEFAULTCOLOR WHITE

typedef  enum  { DISK=0, TEMPORARY, MEMORY }              TARGETTYPE;

typedef  enum  { WHOLEIMAGE=0, BUFACCESS=1, AUTOMATIC=2 } BUFFERINGYESNO;
typedef  enum  { EXPANDED=0, PACKED=1 }                   MODETYPE;
typedef  BYTE*                                            BINLINE;

typedef  struct { int            Left;
                  int            Right;
                  int            Upper;
                  int            Lower;
                } BOUNDS;

typedef  struct { FILE*          FilePointer;
                  TARGETTYPE     TargetType;
                  int            ImageSizeX;
                  int            ImageSizeY;
                  int            ImageBytesX;
                  BOUNDS         Bound;
                  FILETYPE       FileType;
                  FILEMODE       FileMode;
                  BUFFERINGYESNO Buffering;
                  MODETYPE       ModType;
                  int            DefaultColor;
                  int            LinesInBuffer;
                  int            ScrollSize;
                  int            FirstBufLine;
                  int            LastBufLine;
                  int            BufCount;
                  int            LinesReadTotal;
                  BINLINE*       data;
                  } IMAGE;

/*-----------------------  ERROR CODES  -------------------------*/

typedef enum {
  BI_OTHER=0, BI_NOMEMORY, IMSIZE, IMBOUNDS, TRBLOCK, UBIMAGEMODEL,
  LINEPASSED, LINESTORED, UFMODE, UBUFFERING, LARGEBUFFER }
  BINIMAGE_ERROR;

extern BINIMAGE_ERROR BIErrorStatus();

/*-----------------------  General routines  -------------------------*/

extern YESNO ImageInit(char InputName[], IMAGE* image, FILEMODE mode,
                     int AllowOverWrite, BUFFERINGYESNO Bufferingyesno,
                     MODETYPE ModType, int BufferSize, int DefaultColor);

extern YESNO ImageInitScrolled(char InputName[], IMAGE* image,
                               FILEMODE mode, int AllowOverWrite,
                               BUFFERINGYESNO Bufferingyesno,
                               MODETYPE ModType, int BufferSize,
                               int ScrollSize, int DefaultColor);

/* File will be erased after closing - not yet realesed */
extern YESNO TmpImageInit(char InputName[], IMAGE* image,
                     int AllowOverWrite, BUFFERINGYESNO Bufferingyesno,
                     MODETYPE ModType, int BufferSize, int DefaultColor);

/* Image will be stored in memory only */
extern YESNO MemImageInit(IMAGE* image, int DefaultColor);

extern YESNO ImageDone(IMAGE* image);
extern YESNO ImageRewind(IMAGE* image);

extern YESNO SetImageSize(IMAGE* image, int x, int y);
extern YESNO SetImageBounds(IMAGE* image, int bx, int by, int cx, int cy);
extern YESNO EmptyWholeImage(IMAGE* image);

/*------------------  Pixel level; absolute  -------------------------*/

extern int  GetImagePixel(IMAGE* image, int x, int y);
extern void PutImagePixel(IMAGE* image, int x, int y, int value);
extern int  GetImageBitPixel(IMAGE* image, int x, int y);
extern void PutImageBitPixel(IMAGE* image, int x, int y, int value);
extern int GetUnboundedImageBitPixel(IMAGE* image, int x, int y);

/*------------------  Pixel level; relative  -------------------------*/

extern int GetRelativeImageBitPixel(IMAGE* image, int x, int y, int deltax, int deltay);

/*--------------------------------------------------------------------*/

extern int BytesPerLine(int PelsPerLine, int packedbits);

/*--------------------------------------------------------------------*/

  /* For macro: */
  /* int by,bx,tmp; */
  /* BYTE value,mask; */

#define GetImageBitPixel_(bit,image,x,y)               \
{                                                      \
  if ( x<image->Bound.Left || x>image->Bound.Right ||  \
       y<image->Bound.Upper || y>image->Bound.Lower )  \
  { bit=(image->DefaultColor); }                       \
  else                                                 \
  {                                                    \
  by= (y - image->FirstBufLine + image->BufCount);     \
  if ( (tmp=by-image->LinesInBuffer) >= 0) by=tmp;     \
                                                       \
  bx= (x-1)>>3;                                        \
  mask= 0x80;                                          \
  value= image->data[by][bx];                          \
  if ( (tmp= ((x-1) & 7)) > 0 )                        \
     mask>>=tmp;                                       \
                                                       \
  bit = (value & mask) != 0 ? BLACK : WHITE;           \
  }                                                    \
}                                                      \

#endif /* __BINIMAGE_H */
