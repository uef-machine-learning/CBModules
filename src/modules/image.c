/*-------------------------------------------------------------------*/
/* IMAGE.C      Pasi Fr„nti                                          */
/*              Timo Kaukoranta                                      */
/*                                                                   */
/* Data structures and interface to image buffer.                    */
/* The image can be processed in the following ways:                 */
/*   - Reading whole image into memory                               */
/*   - Using scrolling buffer of a fixed number of lines             */
/*                                                                   */
/* Image module structure has been changed. See the changes to       */
/* ver. 0.07 in image.h.                                             */
/*-------------------------------------------------------------------*/

#define ProgName        "IMAGE"
#define VersionNumber   "Version 0.17"
#define LastUpdated     "02.07.2001"

/* ----------------------------------------------------------------- */
/*      Changelog:                                                   */
/*                                                                   */
/* PF  0.16  Rewind: Bug corrected (read first line two times).      */
/* PK  0.17  Rewind fixed                                            */
/* ------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "interfc.h"
#include "file.h"
#include "image.h"
#include "memctrl.h"
#include "pgm.h"


/* -------------------------- Misc -----------------------------------*/

#if defined(max) || defined(min)
#undef max
#undef min
#endif
#define  max(a,b) ((a) > (b) ? (a) : (b))
#define  min(a,b) ((a) < (b) ? (a) : (b))
#define  ScaleBetween(a,b,c)   max((a), min((b), (c)))


/*======================  Memory allocation  ========================*/


static void AllocateMemoryForImage(IMAGE* image)
{
  int  LineAmount  = (image->ImageSizeX+1) * sizeof(PIXELTYPE);
  int  PlaneAmount;
  int  y, z, n;

  switch( image->FileType )
    {
    case PPM:         PlaneAmount = 3; break;
    case PGMEXTENDED: PlaneAmount = 2; break;
    default:          PlaneAmount = 1; break;
    }

  n = image->LinesInBuffer;
  if( n >= LinesMax )
	{
    ErrorMessage("ERROR: Image module supports image sizes up to %i lines only.\n",LinesMax);
    ExitProcessing(-1);
	}

  /* Allocate memory for each line for each plane */
  for( z = 0; z < PlaneAmount; z++ ) /* Added: Loop for each plane (0-2) */
    {
    for( y = 1; y <= n; y++ )
      {
      image->data[z][y] = (PIXELTYPE*) allocate( LineAmount ); /* Changed */
      if( image->data[z][y] == NULL )
        {
        ErrorMessage("\nERROR: allocate image->data\n");
        ExitProcessing(-1);
        }
      }
    }
}


/*-------------------------------------------------------------------*/


void FreeImage(IMAGE* image)
{
  int  y, z, n;
  int  PlaneAmount;

  switch( image->FileType )
    {
    case PPM:         PlaneAmount = 3; break;
    case PGMEXTENDED: PlaneAmount = 2; break;
    default:          PlaneAmount = 1; break;
    }

  n = image->LinesInBuffer;

  for( z = 0; z < PlaneAmount; z++ )
  {
    for( y = 1; y <= n; y++ )
    {
      deallocate(image->data[z][y]);
    }
  }
}


/*===================  Buffer operation routines  ===================*/


BYTE GetImagePixel(IMAGE* image, int x, int y, int z)
{
  x = ScaleBetween(1, x, image->ImageSizeX);
  y = ScaleBetween(1, y, image->LinesInBuffer);
  z = ScaleBetween(0, z, 2);
  return( image->data[z][y][x] );	/* Changed */
}


/*-------------------------------------------------------------------*/


void PutImagePixel(IMAGE* image, int x, int y, int z, BYTE value)
{
  if( x<1 || x>image->ImageSizeX || y<1 || y>image->LinesInBuffer || z<0 || z>2 )
  {
    return;
  }
  else
  {
    image->data[z][y][x] = value; /* Changed */
  }
}


/* ------------------------------------------------------------------*/


void GetImagePnmPixel(IMAGE* image, int x, int y, BYTE* R, BYTE* G, BYTE* B)
{
  switch(image->FileType)
	{
    case PBM: (*R) = GetImagePixel(image,x,y,0); break;
    case PGMASCII:
    case PGM: (*R) = GetImagePixel(image,x,y,0); break;
    case PGMEXTENDED:
              (*R) = GetImagePixel(image,x,y,0);
              (*G) = GetImagePixel(image,x,y,1); break;
    case PPM: (*R) = GetImagePixel(image,x,y,0);
              (*G) = GetImagePixel(image,x,y,1);
              (*B) = GetImagePixel(image,x,y,2); break;
    default:  break;
	}
}


/* ------------------------------------------------------------------*/


void PutImagePnmPixel(IMAGE* image, int x, int y, BYTE R, BYTE G, BYTE B)
{
  switch(image->FileType)
	{
    case PBM: PutImagePixel(image,x,y,0,R); break;
    case PGMASCII:
    case PGM: PutImagePixel(image,x,y,0,R); break;
    case PGMEXTENDED:
              PutImagePixel(image,x,y,0,R);
              PutImagePixel(image,x,y,1,G); break;
    case PPM: PutImagePixel(image,x,y,0,R);
              PutImagePixel(image,x,y,1,G);
              PutImagePixel(image,x,y,2,B); break;
    default:  break;
	}
}


/*-------------------------------------------------------------------*/


RGBPIXEL GetImageRGBPixel(IMAGE* image, int x, int y)
{
  RGBPIXEL pixel;

  x = ScaleBetween(1, x, image->ImageSizeX);
  y = ScaleBetween(1, y, image->LinesInBuffer);
  pixel.Red   = image->data[0][y][x];	/* Changed */
  pixel.Green = image->data[1][y][x];	/* Changed */
  pixel.Blue  = image->data[2][y][x];	/* Changed */
  return(pixel);
}


/*-------------------------------------------------------------------*/


void PutImageRGBPixel(IMAGE* image, int x, int y, RGBPIXEL pixel)
{
  if( x<1 || x>image->ImageSizeX || y<1 || y>image->LinesInBuffer )
    {
    return;
    }
  else
    {
    image->data[0][y][x] = pixel.Red;	/* Changed */
    image->data[1][y][x] = pixel.Green;	/* Changed */
    image->data[2][y][x] = pixel.Blue;	/* Changed */
    }
}


/*-------------------------------------------------------------------*/


EXTENDEDPIXEL GetImageExtendedPixel(IMAGE* image, int x, int y)
{
  BYTE R, G ,B;

  GetImagePnmPixel(image, x, y, &R, &G, &B);

  if( image->FileType == PGMEXTENDED )
    {
    return( ((EXTENDEDPIXEL)R << 8) | (EXTENDEDPIXEL)G );
    }
  else
    {
    return( (EXTENDEDPIXEL)R );
    }
}


/*-------------------------------------------------------------------*/


void PutImageExtendedPixel(IMAGE* image, int x, int y, EXTENDEDPIXEL pixel)
{
  BYTE R, G ,B;

  if( x<1 || x>image->ImageSizeX ||
      y<1 || y>image->LinesInBuffer )
    {
    return;
    }
  else
    {
    if( image->FileType == PGMEXTENDED )
      {
      R = (pixel & 0x0000FF00) >> 8;
      G = pixel & 0x000000FF;
      B = 0;
      }
    else
      {
      R = pixel & 0x000000FF;
      G = 0;
      B = 0;
      }
    PutImagePnmPixel(image, x, y, R, G, B);
    }
}


/*===================================================================*/


BYTE GetBufferPixel(IMAGE* image, int x, int y, int z)
{
  int line;

  x = ScaleBetween(1, x, image->ImageSizeX);
  z = ScaleBetween(0, z, 2);
  line = image->CurrentLine - 1 + y + image->LinesInBuffer;
  line = 1 + (line % image->LinesInBuffer);
  return( image->data[z][line][x] );	/* Changed */
}


/*-------------------------------------------------------------------*/


void PutBufferPixel(IMAGE* image, int x, int y, int z, BYTE value)
{
  int line;

  line = image->CurrentLine - 1 + y + image->LinesInBuffer;
  line = 1 + (line % image->LinesInBuffer);
  if( x<1 || x>image->ImageSizeX || y<image->PreviousLines ||
	 y>image->FutureLines || z<0 || z>2 )
     {
     return;
     }
  else
     {
     image->data[z][line][x] = value;	/* Changed */
     }
}


/*-------------------------------------------------------------------*/


void GetBufferPnmPixel(IMAGE* image, int x, int y, BYTE* R, BYTE* G, BYTE* B)
{
  switch(image->FileType)
	{
    case PBM: (*R) = GetBufferPixel(image,x,y,0); break;
    case PGMASCII:
    case PGM: (*R) = GetBufferPixel(image,x,y,0); break;
    case PGMEXTENDED:
              (*R) = GetBufferPixel(image,x,y,0);
              (*G) = GetBufferPixel(image,x,y,1); break;
    case PPM: (*R) = GetBufferPixel(image,x,y,0);
              (*G) = GetBufferPixel(image,x,y,1);
              (*B) = GetBufferPixel(image,x,y,2); break;
    default:  break;
	}
}


/*-------------------------------------------------------------------*/


void PutBufferPnmPixel(IMAGE* image, int x, int y, BYTE  R, BYTE  G, BYTE  B)
{
  switch(image->FileType)
	{
    case PBM: PutBufferPixel(image,x,y,0,R); break;
    case PGMASCII:
    case PGM: PutBufferPixel(image,x,y,0,R); break;
    case PGMEXTENDED:
              PutBufferPixel(image,x,y,0,R);
              PutBufferPixel(image,x,y,1,G); break;
    case PPM: PutBufferPixel(image,x,y,0,R);
              PutBufferPixel(image,x,y,1,G);
              PutBufferPixel(image,x,y,2,B); break;
    default:  break;
	}
}


/*-------------------------------------------------------------------*/


RGBPIXEL GetBufferRGBPixel(IMAGE* image, int x, int y)
{
  RGBPIXEL pixel;
  int      line;

  x    = ScaleBetween(1, x, image->ImageSizeX);
  line = image->CurrentLine - 1 + y + image->LinesInBuffer;
  line = 1 + (line % image->LinesInBuffer);
  pixel.Red   = image->data[0][line][x];
  pixel.Green = image->data[1][line][x];
  pixel.Blue  = image->data[2][line][x];
  return(pixel);
}


/*-------------------------------------------------------------------*/


void PutBufferRGBPixel(IMAGE* image, int x, int y, RGBPIXEL pixel)
{
  int line;

  line = image->CurrentLine - 1 + y + image->LinesInBuffer;
  line = 1 + (line % image->LinesInBuffer);
  if( x<1 || x>image->ImageSizeX || y<image->PreviousLines || y>image->FutureLines )
     {
     return;
     }
  else
    {
	image->data[0][line][x] = pixel.Red;
	image->data[1][line][x] = pixel.Green;
	image->data[2][line][x] = pixel.Blue;
    }
}


/*-------------------------------------------------------------------*/


EXTENDEDPIXEL GetBufferExtendedPixel(IMAGE* image, int x, int y)
{
  BYTE R, G ,B;

  GetBufferPnmPixel(image, x, y, &R, &G, &B);

  if( image->FileType == PGMEXTENDED )
    {
    return( ((EXTENDEDPIXEL)R << 8) | (EXTENDEDPIXEL)G );
    }
  else
    {
    return( (EXTENDEDPIXEL)R );
    }
}


/*-------------------------------------------------------------------*/


void PutBufferExtendedPixel(IMAGE* image, int x, int y, EXTENDEDPIXEL pixel)
{
  BYTE R, G ,B;

  if( x<1 || x>image->ImageSizeX ||
      y<1 || y>image->LinesInBuffer )
    {
    return;
    }
  else
    {
    if( image->FileType == PGMEXTENDED )
      {
      R = (pixel & 0x0000FF00) >> 8;
      G = pixel & 0x000000FF;
      B = 0;
      }
    else
      {
      R = pixel & 0x000000FF;
      G = 0;
      B = 0;
      }
    PutBufferPnmPixel(image, x, y, R, G, B);
    }
}


/*-------------------------------------------------------------------*/


static void ReadBuffer(IMAGE* image, int line)
{
  int   x;
  BYTE  R=0, G=0, B=0;

  (image->LinesReadTotal)++;
  for(x=1; x<=image->ImageSizeX; x++)
    {
    ReadPnmPixel(image->FilePointer, image->FileType, &R, &G, &B);
    PutImagePnmPixel(image,x,line,R,G,B);
    }
  if(image->FileType==PBM && image->FileMode!=OUTPUT)
    {
    FlushInputBitBuffer(image->FilePointer);
    }
}


/*-------------------------------------------------------------------*/


static void WriteBuffer(IMAGE* image, int line)
{
  int   x;
  BYTE  R=0, G=0, B=0;

  (image->LinesReadTotal)++;
  for(x=1; x<=image->ImageSizeX; x++)
    {
    GetImagePnmPixel(image, x, line, &R, &G, &B);
    WritePnmPixel(image->FilePointer, image->FileType, R, G, B);
    }
  if(image->FileType==PBM && image->FileMode!=INPUT)
    {
    FlushOutputBitBuffer(image->FilePointer);
    }
}


/*-------------------------------------------------------------------*/


static void MakeEmptyLine(IMAGE* image, int line)
{
  int x;

  for(x=1; x<=image->ImageSizeX; x++)
    {
    PutImagePnmPixel(image,x,line,0,0,0);
    }
}


/*-------------------------------------------------------------------*/


static void CopyBufferLine(IMAGE* image, int source, int destination)
{
  int   x;
  BYTE  R=0, G=0, B=0;

  for(x=1; x<=image->ImageSizeX; x++)
    {
    GetImagePnmPixel(image, x, source, &R, &G, &B);
    PutImagePnmPixel(image, x, destination, R, G, B);
    }
}


/*-------------------------------------------------------------------*/


void CopyLine(IMAGE* source, IMAGE* destination, int line)
{
  int   x;
  BYTE  R=0, G=0, B=0;

  for(x=1; x<=source->ImageSizeX; x++)
    {
    GetImagePnmPixel(source, x, line, &R, &G, &B);
    PutImagePnmPixel(destination, x, line, R, G, B);
    }
}


/*=======================  General routines  ========================*/


void CopyImageHeader(IMAGE* source, IMAGE* destination)
{
  destination->FilePointer    = NULL;
  destination->ImageSizeX     = source->ImageSizeX;
  destination->ImageSizeY     = source->ImageSizeY;
  destination->MaxValue       = source->MaxValue;
  destination->FileType       = source->FileType;
  destination->FileMode       = source->FileMode;
  destination->Buffering      = source->Buffering;
  destination->LinesInBuffer  = 0;
  destination->PreviousLines  = 0;
  destination->FutureLines    = 0;
  destination->CurrentLine    = 1;
  destination->LinesReadTotal = 0;
}


/*-------------------------------------------------------------------*/


void SetImageSize(IMAGE* image, int x, int y)
{
  image->ImageSizeX = x;
  image->ImageSizeY = y;
}


/*-------------------------------------------------------------------*/


void InitializeImage(FILE* f, IMAGE* image, FILEMODE mode,
                     BUFFERINGTYPE BufferType, int BufferSize)
{
  if(mode==INPUT)
    {
    ReadPnmHeader(f, &(image->ImageSizeX), &(image->ImageSizeY),
                     &(image->MaxValue),   &(image->FileType));
    image->HeaderSize = ftell(f);
    }
  else if(mode==OUTPUT)
    {
    WritePnmHeader(f, image->ImageSizeX, image->ImageSizeY,
                      image->MaxValue,   image->FileType);
    image->HeaderSize = ftell(f);
    }
  image->FilePointer    = f;
  image->Buffering      = BufferType;
  switch(BufferType)
    {
    case WHOLEIMAGE:    image->LinesInBuffer  = image->ImageSizeY;
                        image->PreviousLines  = 0;
                        image->FutureLines    = 0;
                        image->CurrentLine    = 1;
                        break;
    case WINDOWBUFFER:  image->PreviousLines  = BufferSize / 2;
                        image->FutureLines    = (BufferSize-1) / 2;
                        image->LinesInBuffer  = BufferSize;
                        image->CurrentLine    = (BufferSize / 2) + 1;
                        break;
    case BLOCKBUFFER:   image->LinesInBuffer  = BufferSize;
                        image->PreviousLines  = 0;
                        image->FutureLines    = BufferSize - 1;
                        image->CurrentLine    = 1;
                        break;
    }
  image->FileMode = mode;
  image->LinesReadTotal = 0;
  AllocateMemoryForImage(image);
}

/*-------------------------------------------------------------------*/

void InitializeInternalImage(IMAGE* image, int x, int y, FILETYPE FileType) {
	SetImageSize(image, x, y);
	image->FileType = FileType;
	InitializeImage(0, image, INTERNAL, WHOLEIMAGE, 0);
}


/*===================  Image I/O - whole image  =====================*/


void ReadWholeImage(IMAGE* image)
{
  int  y;

  for(y=1; y<=image->ImageSizeY; y++)
    {
    ReadBuffer(image, y);
    }
}


/* ------------------------------------------------------------------*/


void WriteWholeImage(IMAGE* image)
{
  int  y;

  for(y=1; y<=image->ImageSizeY; y++)
    {
    WriteBuffer(image, y);
    }
}


/* ------------------------------------------------------------------*/


void SaveImageToFile(IMAGE* image, FILE* f)
{
  image->FileType = OUTPUT;
  image->FilePointer = f;
  WritePnmHeader(f, image->ImageSizeX, image->ImageSizeY,
                    image->MaxValue,   image->FileType);
  WriteWholeImage(image);
}


/*==================  Image I/O - window buffer  ====================*/


static void InitializeInputBuffer(IMAGE* image)
{
  int  i;
  int  c=image->CurrentLine;

  /* Read current line */
  ReadBuffer(image, c);

  /* Fill in FutureLines in buffer */
  for(i=1; i<=image->FutureLines; i++)
	 {
	 ReadBuffer(image, c+i);
	 }

  /* Copy current line to PreviousLines in buffer */
  for(i=1; i<=image->PreviousLines; i++)
	 {
	 CopyBufferLine(image, c, i);
	 }

  /* Just checking... */
  if(image->LinesReadTotal != image->LinesInBuffer - image->PreviousLines )
	{
    ErrorMessage("ERROR: Initialize buffer.\n");
    ExitProcessing(-1);
	}
}


/*-------------------------------------------------------------------*/


void ScrollBuffer(IMAGE* image)
{
  int  n,p;

  if( image->LinesReadTotal == 0 )
    {
    InitializeInputBuffer(image);
    return;
    }
  image->CurrentLine = 1 + image->CurrentLine  % image->LinesInBuffer;
  n = image->CurrentLine + image->FutureLines;
  n = (n > image->LinesInBuffer) ? n-image->LinesInBuffer : n;
  if( image->LinesReadTotal < image->ImageSizeY )
    {
    ReadBuffer(image, n);
    }
  else if(image->FutureLines > 0)
  /* All lines read => Make copy of the last line */
    {
    p = (n==1 ? image->LinesInBuffer : n-1);
    CopyBufferLine(image, p, n);
    }
}


/*-------------------------------------------------------------------*/


void RewindImage(IMAGE* image)
/* Rewinds file pointer at the beginning of image data (after header) */
{
  fseek(image->FilePointer, image->HeaderSize, SEEK_SET);
  image->LinesReadTotal = 0;
}


/*-------------------------------------------------------------------*/


void WriteCurrentLine(IMAGE* image)
{
  WriteBuffer(image, image->CurrentLine);
}


/*==================  Image I/O - block buffer	=====================*/


void ReadNewBuffer(IMAGE* image)
{
  int i;

  for(i=1; i<=image->LinesInBuffer; i++)
    {
    if( image->LinesReadTotal < image->ImageSizeY )
      {
      ReadBuffer(image, i);
      }
    else
      {
      MakeEmptyLine(image, i);
      }
    }

}


/*-------------------------------------------------------------------*/


void WriteCurrentBuffer(IMAGE* image)
{
  int i;

  for(i=1; i<=image->LinesInBuffer; i++)
	 {
	 if( image->LinesReadTotal < image->ImageSizeY )
	   {
	   WriteBuffer(image, i);
	   }
	}
}


/*======================  Abstractions  ==========================*/


int SameColor(RGBPIXEL v1, RGBPIXEL v2)
{
  if(v1.Red!=v2.Red)     return(0);
  if(v1.Green!=v2.Green) return(0);
  if(v1.Blue!=v2.Blue)   return(0);
  return(1);
}


