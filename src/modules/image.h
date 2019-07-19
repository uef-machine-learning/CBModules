#if ! defined(__IMAGE_H)
#define __IMAGE_H

#include "file.h"
#include "pgm.h"
#include "owntypes.h"

/*---------------------  Image data structure  -----------------------*/

#define  LinesMax        2048
#define  NameMaxLength    127
#define  PlanesMax          3

typedef  enum  {  WHOLEIMAGE=0, WINDOWBUFFER=1, BLOCKBUFFER=2 }  BUFFERINGTYPE;
typedef  BYTE                             PIXELTYPE;
typedef  int                              EXTENDEDPIXEL;
typedef  struct { BYTE Red,Green,Blue; }  RGBPIXEL;

typedef  struct { FILE*          FilePointer;
                  int            ImageSizeX;
                  int            ImageSizeY;
                  int            MaxValue;
                  FILETYPE       FileType;
                  FILEMODE       FileMode;
                  BUFFERINGTYPE  Buffering;
                  int            LinesInBuffer;
                  int            PreviousLines;
                  int            FutureLines;
                  int            CurrentLine;
                  int            LinesReadTotal;
                  int            HeaderSize;                /*P. Kopylov */
                  PIXELTYPE*     data[PlanesMax][LinesMax]; /* Changed TKa */
                } IMAGE;

/*-------------------------  Abstractions  ---------------------------*/

#define   ImageSize(Image)    ((long)Image->ImageSizeX * Image->ImageSizeY)
#define   ColorBands(Image)   ((Image)->FileType==PPM ? 3 : 1)
int       SameColor(RGBPIXEL v1, RGBPIXEL v2);

/*-----------------------  General routines  -------------------------*/

void FreeImage(IMAGE* image);
void CopyImageHeader(IMAGE* source, IMAGE* destination);
void SetImageSize(IMAGE* image, int x, int y);
void InitializeImage(FILE* f, IMAGE* image, FILEMODE mode,
                     BUFFERINGTYPE BufferType, int BufferSize);
void CopyLine(IMAGE* source, IMAGE* destination, int line);
void InitializeInternalImage(IMAGE* image, int x, int y, FILETYPE FileType);

/*-------------------  I/O routines; whole image  --------------------*/

void ReadWholeImage(IMAGE* image);
void WriteWholeImage(IMAGE* image);
void SaveImageToFile(IMAGE* image, FILE* f);

/*-----------------  I/O routines; window buffering  -----------------*/

void ScrollBuffer(IMAGE* image);
void WriteCurrentLine(IMAGE* image);
void RewindImage(IMAGE* image);

/*-----------------  I/O routines; block buffering  ------------------*/

void ReadNewBuffer(IMAGE* image);
void WriteCurrentBuffer(IMAGE* image);

/*------------------  Pixel level; absolute  -------------------------*/

BYTE GetImagePixel(IMAGE* image, int x, int y, int z);
void PutImagePixel(IMAGE* image, int x, int y, int z, BYTE value);
void GetImagePnmPixel(IMAGE* image, int x, int y, BYTE* R, BYTE* G, BYTE* B);
void PutImagePnmPixel(IMAGE* image, int x, int y, BYTE  R, BYTE  G, BYTE  B);
RGBPIXEL GetImageRGBPixel(IMAGE* image, int x, int y);
void     PutImageRGBPixel(IMAGE* image, int x, int y, RGBPIXEL pixel);
EXTENDEDPIXEL GetImageExtendedPixel(IMAGE* image, int x, int y);
void          PutImageExtendedPixel(IMAGE* image, int x, int y, EXTENDEDPIXEL pixel);

#define  GetImageGray(image,x,y)        GetImagePixel(image,x,y,0)
#define  GetImageRed(image,x,y)         GetImagePixel(image,x,y,0)
#define  GetImageGreen(image,x,y)       GetImagePixel(image,x,y,1)
#define  GetImageBlue(image,x,y)        GetImagePixel(image,x,y,2)

#define  PutImageGray(image,x,y,v)      PutImagePixel(image,x,y,0,v)
#define  PutImageRed(image,x,y,v)       PutImagePixel(image,x,y,0,v)
#define  PutImageGreen(image,x,y,v)     PutImagePixel(image,x,y,1,v)
#define  PutImageBlue(image,x,y,v)      PutImagePixel(image,x,y,2,v)

/*-----------  Pixel level; relative to current line (0)  ------------*/

BYTE GetBufferPixel(IMAGE* image, int x, int y, int z);
void PutBufferPixel(IMAGE* image, int x, int y, int z, BYTE value);
void GetBufferPnmPixel(IMAGE* image, int x, int y, BYTE* R, BYTE* G, BYTE* B);
void PutBufferPnmPixel(IMAGE* image, int x, int y, BYTE  R, BYTE  G, BYTE  B);
RGBPIXEL GetBufferRGBPixel(IMAGE* image, int x, int y);
void     PutBufferRGBPixel(IMAGE* image, int x, int y, RGBPIXEL pixel);
EXTENDEDPIXEL GetBufferExtendedPixel(IMAGE* image, int x, int y);
void          PutBufferExtendedPixel(IMAGE* image, int x, int y, EXTENDEDPIXEL pixel);

#define  GetBufferGray(image,x,y)       GetBufferPixel(image,x,y,0)
#define  GetBufferRed(image,x,y)        GetBufferPixel(image,x,y,0)
#define  GetBufferGreen(image,x,y)      GetBufferPixel(image,x,y,1)
#define  GetBufferBlue(image,x,y)       GetBufferPixel(image,x,y,2)

#define  PutBufferGray(image,x,y,v)     PutBufferPixel(image,x,y,0,v)
#define  PutBufferRed(image,x,y,v)      PutBufferPixel(image,x,y,0,v)
#define  PutBufferGreen(image,x,y,v)    PutBufferPixel(image,x,y,1,v)
#define  PutBufferBlue(image,x,y,v)     PutBufferPixel(image,x,y,2,v)

/*--------------------------------------------------------------------*/

#endif /* __IMAGE_H */

