/*-------------------------------------------------------------------*/
/* BINIMAGE.C     Eugene Ageenko                                     */
/*                                                                   */
/* Data structures and interface to binary iamge (PBM file format)   */
/*                                                                   */
/* Image can be opened for INPUT or OUTPUT, it associates with       */
/* given FileName (with extension)                                   */
/*                                                                   */
/* If FileName == TMPNAME then temporary file for image is created   */
/* You may use Rewind() to rewind temporary image (then it change    */
/* its mode to INPUT). Note: after ImageDone() temporary image will  */
/* left!!! - WARNING!!! Doesn't working now.                         */
/*                                                                   */
/* The image can be processed in the following ways:                 */
/*   - Reading whole image into memory (WHOLEIMAGE)                  */
/*   - Using automagic scrolling buffer of a fixed number of lines   */
/*     (BUFACCESS, number of lines is given)                         */
/*   - buffer size is automaticly reduced to fit the image size      */
/*   - AUTOMATIC - buffer size is choosen automaticly, in the case   */
/*     of "Out of memory" problem it is reduced.                     */
/*                                                                   */
/* The image (opened for INPUT) may be rewinded to begin.            */
/*                                                                   */
/* Two image representantion methods are supported:                  */
/*   - one pixel per byte (eptended model)  EXPANDED                 */
/*   - eight pixels per byte (packed model) PACKED                   */
/*                                                                   */
/* For each representation you must use functions                    */
/*   GetImagePixel (PutImagePixel) or                                */
/*   GetImageBitPixel (PutImageBitPixel) respectively                */
/*                                                                   */
/*-------------------------------------------------------------------*/

#define ProgName        "BINIMAGE"
#define VersionNumber   "0.17"
#define LastUpdated     "21.11.97"

/* ------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "file.h"
#include "binimage.h"
#include "pgm.h"
#include "memctrl.h"
#include "error.h"

/* -------------------------- Misc -----------------------------------*/

#if defined(max) || defined(min)
#undef max
#undef min
#endif
#define  max(a,b) ((a) > (b) ? (a) : (b))
#define  min(a,b) ((a) < (b) ? (a) : (b))
#define  ScaleBetween(a,b,c)   max((a), min((b), (c)))

/*----------------- Internal definitiones  -----------------------------------*/

static void ReadLines(IMAGE* image, int n);
static void WriteLines(IMAGE* image, int n);
static void EmptyLine(IMAGE* image, int line);
static void EmptyCurrentLine(IMAGE* image);

/* ====================== internal ==================================*/

int BytesPerLine(int PelsPerLine, int packedbits)

{

   return( (packedbits==PACKED) ? ( 1 + (int) ((PelsPerLine-1)/8) ) : PelsPerLine );

}

/* ==================== ERROR HANDLING ==============================*/

/* In the case of serious internal error special error function is   */
/* called from this module to hadle the error situations, default it */
/* output message (through the Errormsg), breaks the program and     */
/* exit to OS, otherwise (non-serious) error,                        */
/* function returns NO and error status can be obtained from the     */
/* BIErrorStatus()                                                   */

BINIMAGE_ERROR BIErrorStatus_ = 0;

BINIMAGE_ERROR BIErrorStatus() {return BIErrorStatus_;}

#if defined(ERROR_)
#undef ERROR_
#endif
#if defined(ERROR_V)
#undef ERROR_V
#endif

#define  ERROR_(x) bi_error((x),0,ProgName)
#define  ERROR_V(x,value) bi_error((x),(value),ProgName)

void bi_error(int code, int value, char* module)

{
  BIErrorStatus_ = (BINIMAGE_ERROR) code;
  switch (code)
  {
    case BI_NOMEMORY:   ErrorMsg ("(%s) Not enough memory to complete operation, aborted.",module); break;
    case IMSIZE:        ErrorMsg ("(%s) Image size must be positive number.",module); break;
    case IMBOUNDS:      ErrorMsg ("(%s) Bounds can not exseed image size.",module); break;
    case TRBLOCK:       ErrorMsg ("(%s) Block can not be transpose.",module); break;
    case UBIMAGEMODEL:  ErrorMsg ("(%s) Unknow image representation model.",module); break;
    case LINEPASSED:    ErrorMsg ("(%s) Line %i already left buffer and can't be observed.",module,value); break;
    case LINESTORED:    ErrorMsg ("(%s) Line %i already stored and can't be observed.",module,value); break;
    case UFMODE:        ErrorMsg ("(%s) Unknow file open mode.",module); break;
    case UBUFFERING:    ErrorMsg ("(%s) Unknow image buffering method.",module); break;
    case LARGEBUFFER:   ErrorMsg ("(%s) Module supports buffer size up to %i lines only.",module,value); break;
    default:            ErrorMsg ("(%s) Internal error - Call program vendor",module);
  }
}


/* ==================  Internal I/O routines  ========================*/

void EmptyLine(IMAGE* image, int line)
/* fill line with DeafultColor - now you may use PutPixel ;) */
{
  int x;
  BYTE filler;

  if (image->ModType == OBB)
     filler = (image->DefaultColor == WHITE) ? 0 : 1;
  else
     filler = (image->DefaultColor == WHITE) ? 0 : 0xFF;

  for(x=0; x<image->ImageBytesX; x++)
  {
        image->data[line][x]=filler;
  }
}

/* ------------------------------------------------------------------*/

void EmptyCurrentLine(IMAGE* image)
/* fill current line with DeafultColor */
{
  int line;
  line=image->BufCount;
  EmptyLine(image,line);
}

/* ------------------------------------------------------------------*/

YESNO EmptyWholeImage(IMAGE* image)
/* fill current WHOLEIMAGE image with DefaultColor */

{
  int y;

  if (image->Buffering == WHOLEIMAGE)
  {
     for(y=0; y<image->LinesInBuffer; y++)
     {
        EmptyLine(image,y);
     }
     return (YES);
  }
  else return (NO);
}

/* ------------------------------------------------------------------*/

void ReadLines(IMAGE* image, int n)
{
  int k,tmp;

  for (k=0; k<n; k++)
  {
     switch (image->ModType) {
     case OBB: ReadPbmLine(image->FilePointer, image->data[image->BufCount],
                           image->ImageBytesX); break;
     default: /* case EBB */
               ReadPbmBitLine(image->FilePointer, image->data[image->BufCount],
                              image->ImageBytesX); break;
     } /* switch */
     /* INPUT mode doesn't cheked. Be carefull */

     image->BufCount++;
     if ( (tmp=image->BufCount - image->LinesInBuffer) >= 0 )
        image->BufCount = tmp;

  }

  image->FirstBufLine += n;
  image->LastBufLine  += n;
  image->LinesReadTotal += n;

}

/* ------------------------------------------------------------------*/


void WriteLines(IMAGE* image, int n)
{
  int k,tmp;

  for (k=0; k<n; k++)
  {
     switch (image->ModType)
     {
        case OBB: WritePbmLine(image->FilePointer, image->data[image->BufCount],
                               image->ImageBytesX); break;
        default: /* case EBB */
                  WritePbmBitLine(image->FilePointer, image->data[image->BufCount],
                                  image->ImageBytesX); break;
     } /* switch */
     EmptyCurrentLine(image);
     /* OUTPUT mode doesn't cheked. Be carefull */
     image->BufCount++;
     if ( (tmp=image->BufCount - image->LinesInBuffer) >= 0 )
        image->BufCount = tmp;

  }

  image->FirstBufLine += n;
  image->LastBufLine  += n;
  image->LinesReadTotal += n;

}


/*==================== Image-Interface ==============================*/
/*-------------------- Get Pixel ------------------------------------*/

int GetImagePixel(IMAGE* image, int x, int y)
{
  int by,tmp;

  if ( x<image->Bound.Left || x>image->Bound.Right ||
       y<image->Bound.Upper || y>image->Bound.Lower )
  return(image->DefaultColor);

  if (y < image->FirstBufLine)
  {
     ERROR_V(LINEPASSED,y);
     exit(-1);
  }

  if ( (tmp=y-image->LastBufLine) > 0)
  {
     tmp += min(image->ScrollSize,image->ImageSizeY-y);
     if (image->FileMode==OUTPUT) WriteLines (image, tmp);
     else ReadLines (image, tmp);
  }

  by= (y - image->FirstBufLine + image->BufCount);
  if ( (tmp=by-image->LinesInBuffer) >= 0) by=tmp;

/* by= (y - image->FirstBufLine + image->BufCount) % image->LinesInBuffer */

  return( (int) (image->data[by][x-1] ) );
}


/*-------------------------------------------------------------------*/


int GetImageBitPixel(IMAGE* image, int x, int y)
{
  int by,bx,tmp;
  BYTE value,mask;

  if ( x<image->Bound.Left || x>image->Bound.Right ||
       y<image->Bound.Upper || y>image->Bound.Lower )
  return(image->DefaultColor);

  if (y < image->FirstBufLine)
  {
     ERROR_V(LINEPASSED,y);
     exit(-1);

  }

  if ( (tmp=y-image->LastBufLine) > 0)
  {
     tmp += min(image->ScrollSize,image->ImageSizeY-y);
     if (image->FileMode==OUTPUT) WriteLines (image, tmp);
     else ReadLines (image, tmp);
  }

  by= (y - image->FirstBufLine + image->BufCount);
  if ( (tmp=by-image->LinesInBuffer) >= 0) by=tmp;

/* by= (y - image->FirstBufLine + image->BufCount) % image->LinesInBuffer */

  bx= (x-1)>>3;
  mask= 0x80;
  value= image->data[by][bx];
  if ( (tmp= ((x-1) & 7)) > 0 )
     mask>>=tmp;

  return ( (value & mask) != 0 ? BLACK : WHITE );
}


/*-------------------------------------------------------------------*/


int GetUnboundedImageBitPixel(IMAGE* image, int x, int y)
{
  int by,bx,tmp;
  BYTE value,mask;

  if ( x<1 || x>image->ImageSizeX ||
       y<1 || y>image->ImageSizeY )
  return(image->DefaultColor);

  if (y < image->FirstBufLine)
  {
     ERROR_V(LINEPASSED,y);
     exit(-1);

  }

  if ( (tmp=y-image->LastBufLine) > 0)
  {
     tmp += min(image->ScrollSize,image->ImageSizeY-y);
     if (image->FileMode==OUTPUT) WriteLines (image, tmp);
     else ReadLines (image, tmp);
  }

  by= (y - image->FirstBufLine + image->BufCount);
  if ( (tmp=by-image->LinesInBuffer) >= 0) by=tmp;

/* by= (y - image->FirstBufLine + image->BufCount) % image->LinesInBuffer */

  bx= (x-1)>>3;
  mask= 0x80;
  value= image->data[by][bx];
  if ( (tmp= ((x-1) & 7)) > 0 )
     mask>>=tmp;

  return ( (value & mask) != 0 ? BLACK : WHITE );
}


/*-------------------- Get Relative Pixel ---------------------------*/


int GetRelativeImageBitPixel(IMAGE* image, int x, int y, int deltax, int deltay)
{
  int x1 = x + deltax;
  int y1 = y + deltay;

/*  if ( ( x==image->Bound.Left || x==image->Bound.Right ) &&
       ( ( deltay >= 0 ) || ( y==image->Bound.Upper || y==image->Bound.Lower ) ) )
     return (image->DefaultColor);

  if ( ( y==image->Bound.Upper || y==image->Bound.Lower ) &&
       ( deltax >= 0 ) )
     return(image->DefaultColor);
*/

/*  if ( ( x==image->Bound.Left || x==image->Bound.Right ) ||
       ( y==image->Bound.Upper || y==image->Bound.Lower ) )
     return(image->DefaultColor);
*/

/*
  x1 = ScaleBetween(image->Bound.Left,x1,image->Bound.Right);
  y1 = ScaleBetween(image->Bound.Upper,y1,image->Bound.Lower);
*/
  return(GetImageBitPixel(image,x1,y1));

}


/*-------------------- Put Pixel ------------------------------------*/


void PutImagePixel(IMAGE* image, int x, int y, int value)
{

  int by,tmp;

  if ( x<image->Bound.Left || x>image->Bound.Right ||
       y<image->Bound.Upper || y>image->Bound.Lower )
  return;

  if (y < image->FirstBufLine)
  {
     ERROR_V(LINESTORED,y);
     exit(-1);
  }

  if ( (tmp=y-image->LastBufLine) > 0)
  {
     tmp += min(image->ScrollSize,image->ImageSizeY-y);
     WriteLines (image, tmp);
  }

  by= (y - image->FirstBufLine + image->BufCount);
  if ( (tmp=by-image->LinesInBuffer) >= 0) by=tmp;

/* by= (y - image->FirstBufLine + image->BufCount) % image->LinesInBuffer */

  image->data[by][x-1]= (BYTE) value;

}


/*-------------------------------------------------------------------*/


void PutImageBitPixel(IMAGE* image, int x, int y, int value)
{

  int by,bx,tmp;
  BYTE mask;

  if ( x<image->Bound.Left || x>image->Bound.Right ||
       y<image->Bound.Upper || y>image->Bound.Lower )
  return;

  if (y < image->FirstBufLine)
  {
     ERROR_V(LINESTORED,y);
     exit(-1);
  }

  if ( (tmp=y-image->LastBufLine) > 0)
  {
     WriteLines (image, tmp);
     tmp += min(image->ScrollSize,image->ImageSizeY-y);
  }

  by= (y - image->FirstBufLine + image->BufCount);
  if ( (tmp=by-image->LinesInBuffer) >= 0) by=tmp;

/* by= (y - image->FirstBufLine + image->BufCount) % image->LinesInBuffer */

  mask=0x80;
  bx= (x-1)>>3;
  if ( (tmp= ((x-1) & 7)) > 0 )
     mask>>=tmp;
  if (value!= WHITE) image->data[by][bx] |= mask;    /* set bit */
  else               image->data[by][bx] &= ~mask;   /* kill bit */

}


/*=======================  General routines  ========================*/


YESNO SetImageSize(IMAGE* image, int x, int y)
{
  x=max(1,x); y=max(1,y);
  image->ImageSizeX = x;
  image->ImageSizeY = y;
  return (YES); /* values are rounded automaticly */
}


/*-------------------------------------------------------------------*/


YESNO SetImageBounds(IMAGE* image, int bx, int by, int cx, int cy)
{
  int ex,ey;

  if ( (cx<=0) || (cy<=0) )
  {
    ERROR_(TRBLOCK);
    return(NO);
  }
  ex = bx + cx - 1;
  ey = by + cy - 1;

  bx = max(bx,1);
  by = max(by,1);
  ex = min(ex,image->ImageSizeX);
  ey = min(ey,image->ImageSizeY);

  image->Bound.Left  = bx;
  image->Bound.Right = ex;
  image->Bound.Upper = by;
  image->Bound.Lower = ey;

  return (YES); /* values are rounded automaticly */
}


/*-------------------------------------------------------------------*/


void ImageFlush(IMAGE* image)
{
  int  n,rest,k;

/* --- flush buffer (for OUTPUT) --- */

  if (image->FileMode == OUTPUT)
  {

    n=image->LastBufLine - image->FirstBufLine +1;
    rest=image->ImageSizeY - image->LastBufLine;

    WriteLines (image,n); /* flush buffer */

    if (rest>0)
    /* write rest of lines as default color */

     {
       /* make dummy empty line */
       EmptyCurrentLine (image);
       switch (image->ModType) {

       /* store dummy empty line *rest* times */
       case OBB:  for (k=1;k<=rest;k++)
                  {
                   WritePbmLine(image->FilePointer,
                                image->data[image->BufCount],
                                image->ImageBytesX);
                  }
                  break;
       default: /* case EBB */
                  for (k=1;k<=rest;k++)
                  {
                   WritePbmBitLine(image->FilePointer,
                                   image->data[image->BufCount],
                                   image->ImageBytesX);
                  }
                  break;
       } /* switch */
     } /* rest lines */
     fflush (image->FilePointer);
    }

}


/*-------------------------------------------------------------------*/


YESNO ImageRewind(IMAGE* image)
{
  int n;

  if (image->TargetType == DISK)
  {

     ImageFlush(image);
     fseek (image->FilePointer,0,SEEK_SET);

     image->FileMode       = INPUT;

     ReadPbmHeader(image->FilePointer, &(image->ImageSizeX), &(image->ImageSizeY));
     image->ImageBytesX = BytesPerLine(image->ImageSizeX,image->ModType);
     image->Bound.Left  = 1;
     image->Bound.Right = image->ImageSizeX;
     image->Bound.Upper = 1;
     image->Bound.Lower = image->ImageSizeY;

     n=min(image->LinesInBuffer,image->ImageSizeY);

     image->FirstBufLine   = 1;
     image->LastBufLine    = n;
     image->BufCount       = 0;
     ReadLines(image, n);
     image->LinesReadTotal = n;
     image->FirstBufLine   = 1;
     image->LastBufLine    = n;
     image->BufCount       = 0;
     return (YES);
  }
  else if (image->TargetType == MEMORY)
  {
     image->Bound.Left  = 1;
     image->Bound.Right = image->ImageSizeX;
     image->Bound.Upper = 1;
     image->Bound.Lower = image->ImageSizeY;
     image->BufCount       = 0;
     return(YES);
  }
  else return (NO);

}

/*-------------------------------------------------------------------*/

YESNO ImageInitScrolled (char InputName[], IMAGE* image,
                FILEMODE mode, int AllowOverWrite,
                BUFFERINGYESNO Bufferingyesno, MODETYPE ModType,
                int BufferSize, int ScrollSize, int DefaultColor)
{

int LineAmount,i,n;

  image->TargetType = DISK;
  image->ScrollSize = ScrollSize;

  if (InputName[0]!=0)
  {
     if ((image->FilePointer = FileOpen(InputName, mode, AllowOverWrite))==NULL)
        {
           return (NO);
        }
  }
  else {  image->TargetType = MEMORY;  }
  image->FileType    = PBM;

  if (image->TargetType != MEMORY)
     switch (mode) {
     case INPUT:  ReadPbmHeader(image->FilePointer, &(image->ImageSizeX), &(image->ImageSizeY));
                  break;
     default: /* OUTPUT: */
                  WritePbmHeader(image->FilePointer, image->ImageSizeX, image->ImageSizeY);
                  break;
  } /* if-switch */

  image->FileMode       = mode == INPUT ? INPUT : OUTPUT;
  image->Buffering      = Bufferingyesno;
  image->ModType        = ModType == EXPANDED ? EXPANDED : PACKED;
  image->DefaultColor   = DefaultColor == WHITE ? WHITE : BLACK;

  if (BufferSize == 0) BufferSize = LinesAutoMax;
  switch(Bufferingyesno)
  {
      case WHOLEIMAGE:    image->LinesInBuffer  = image->ImageSizeY;
                          break;
      case BUFACCESS:     image->LinesInBuffer  =
                          min (image->ImageSizeY, BufferSize);
                          break;
      default: /* case AUTOMATIC: */
                          image->LinesInBuffer  =
                          min (image->ImageSizeY, BufferSize);
                          break;
  }

  image->ImageBytesX = BytesPerLine(image->ImageSizeX,image->ModType);
  image->Bound.Left  = 1;
  image->Bound.Right = image->ImageSizeX;
  image->Bound.Upper = 1;
  image->Bound.Lower = image->ImageSizeY;

  /* --- allocate memory --- */

  if( image->LinesInBuffer >= LinesMax )
  {
      ERROR_V(LARGEBUFFER,LinesMax);
      return (NO);
  }

  if ( (image->data = (BINLINE*) allocate(
                      (image->LinesInBuffer +1) * sizeof(BINLINE) )
       ) == NULL)

  {
     ERROR_(BI_NOMEMORY);
     exit(-1);
  }

  LineAmount  = (image->ImageBytesX) * sizeof(BYTE);

  for(i=0; i<(image->LinesInBuffer); i++)
  {
     image->data[i] = (BYTE*) allocate( LineAmount );
     if ( image->data[i] == NULL )
     {
         if ( (i>=LinesMin) && (image->Buffering==AUTOMATIC) &&
              (image->FileMode==INPUT) )
         {
            image->LinesInBuffer = i;
            image->Buffering=BUFACCESS; /* switch to buffering mode */
            break;
         }
         else
         {
            ERROR_(BI_NOMEMORY);
            exit(-1);
         }
     }
  }

   /* for AUTOMATIC mode: if initialization done -> switch to BUFACCESS */
   if (image->Buffering == AUTOMATIC) image->Buffering=BUFACCESS;

   n=min(image->LinesInBuffer,image->ImageSizeY);

   image->LinesReadTotal = 0;
   image->FirstBufLine   = 1;
   image->LastBufLine    = n;
   image->BufCount       = 0;

   switch (mode) {

      case INPUT: /* --- fill buffer (for INPUT) --- */
             ReadLines(image, n);
             image->LinesReadTotal = n;
             image->FirstBufLine   = 1;
             image->LastBufLine    = n;
             image->BufCount       = 0;
             break;

      default: /* case OUTPUT: */
               /* --- fill buffer with default color (for OUTPUT) --- */
              for (i=0;i<n;i++)
              {
                 EmptyCurrentLine(image);
                 image->BufCount++;
              }
              image->FirstBufLine   = 1;
              image->LastBufLine    = n;
              image->BufCount       = 0;
              image->LinesReadTotal = 0;
              break;

      } /* switch */

   return (YES);

}

/*-------------------------------------------------------------------*/

YESNO ImageInit (char InputName[], IMAGE* image,
                FILEMODE mode, int AllowOverWrite,
                BUFFERINGYESNO Bufferingyesno, MODETYPE ModType,
                int BufferSize, int DefaultColor)
{
   return(ImageInitScrolled(InputName, image, mode, AllowOverWrite, Bufferingyesno,
             ModType, BufferSize, 1, DefaultColor));

}

/*----------------------- Memory Image Init ---------------------*/

YESNO MemImageInit(IMAGE* image, int DefaultColor)

{
   int value = NO;

   value = ImageInit("", image, OUTPUT, NO, WHOLEIMAGE, PACKED, 0, DefaultColor);
   return (value); /* NO in the case of "image to large to fit the buffer" */

}

/*------------------------  Free Image  -------------------------*/


/* HUOM! HUOM! HUOM! Before use this function, be sure that following IMAGE */
/* fields are initialized at least to NULL value:                           */
/*        TargetType  = NULL or real value                                  */
/*        FilePointer = NULL or real value                                  */
/*        TargetType  = MEMORY or real value                                */
/* These fields are initialized if ImageInit functions was called.          */
/* Otherwise these fields must be initialized manually in order to allove   */
/* GPF-free deinitialization procedure (e.g. deinitialization calld in      */
/* the case of some internal errors)                                        */

YESNO ImageDone(IMAGE* image)
{
  int y;

  if (image==NULL) return (YES);

  switch (image->TargetType)
  {
     case DISK:      ImageFlush(image);
                     if (image->FilePointer != NULL)
                        fclose(image->FilePointer);
                     break;
     case TEMPORARY:
                     if (image->FilePointer != NULL)
                        fclose(image->FilePointer);
                     /* erase file ... */
                     break;
     default:        break;
  }

  /* --- deallocate memory --- */

  if (image->data != NULL)
  {
     for(y=0; y<image->LinesInBuffer; y++)
     {
         if (image->data[y] != NULL) deallocate(image->data[y]);
     }
     deallocate(image->data);
  }

  return(YES);

}

