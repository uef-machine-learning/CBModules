/*-------------------------------------------------------------------*/
/* PGM.C           Pasi Franti                                       */
/*                 Timo Kaukoranta                                   */
/*                 Eugene Ageenko                                    */
/*                                                                   */
/* - PGM interface.                                                  */
/* - PBM interface.                                                  */
/* - PPM interface.                                                  */
/*                                                                   */
/*-------------------------------------------------------------------*/

#define ProgName        "PGM"
#define VersionNumber   "Version 0.37"
#define LastUpdated     "8.10.98"

/* ----------------------------------------------------------------- */

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "interfc.h"
#include "file.h"
#include "owntypes.h"
#include "pgm.h"

/*-------------------------------------------------------------------*/

#define OK YES
#define FAIL NO

#define  FLUSHNOW  999   /* PBM bit buffer control symbol */


/*======================= P G M - I N P U T ===========================*/


YESNO ReadPgmHeader(FILE* f, int* width, int* height, int* maxgray)
{
  int  ch, Id1, Id2;

  Id1 = getc(f);
  Id2 = getc(f);
  if( Id1 != 'P' || Id2 != '5' )
    {
    ErrorMessage("PGM-identification flag not found. (P5 was expected).\n");
    return(NO);
    }

  /* Read the image size and maximum gray scale value */
  fscanf(f, "%i %i", width, height);
  fscanf(f, "%i", maxgray);

  /* Check validity */
  if( (*maxgray) != 255 )
    {
    ErrorMessage("This module supports 256 gray scales only."
           "The current image file has %i levels.\n", (*maxgray)+1 );
    return(NO);
    }

  /* Skip ONE white space, see 'man pgm' */
  ch = getc(f);

  return( (YESNO) OK );
}


/*-------------------------------------------------------------------*/


YESNO ReadASCIIPgmHeader(FILE* f, int* width, int* height, int* maxgray)
{
  int  ch, Id1, Id2;

  Id1 = getc(f);
  Id2 = getc(f);
  if( Id1 != 'P' || Id2 != '2' )
    {
    ErrorMessage("PGM-identification flag not found. (P2 was expected).\n");
    return(NO);
    }

  /* Read the image size and maximum gray scale value */
  fscanf(f, "%i %i", width, height);
  fscanf(f, "%i", maxgray);

  /* Check validity */
  if( (*maxgray) != 255 )
    {
    ErrorMessage("This module supports 256 gray scales only."
           "The current image file has %i levels.\n", (*maxgray)+1 );
    return(NO);
    }

  /* Skip ONE white space, see 'man pgm' */
  ch = getc(f);

  return( (YESNO) OK );
}


/*-------------------------------------------------------------------*/


YESNO ReadPgmLine(FILE* f, BYTE* data, int  x)
{
  int   i;

  for(i=1; i<=x; i++)
    {
    (*data++) = getc(f);
    }
  return( (YESNO) OK );
}


/*-------------------------------------------------------------------*/


YESNO ReadASCIIPgmLine(FILE* f, BYTE* data, int  x)
{
  int   i;
  char  s[10];

  for(i=1; i<=x; i++)
    {
    fscanf(f, "%s", s);
    (*data++) = atoi(s);
    }
  return( (YESNO) OK );
}


/*======================= P G M - O U P U T ===========================*/


YESNO WritePgmHeader(FILE* f, int width, int height, int maxgray)
{
  /* P5 identification flag */
  putc('P', f);
  putc('5', f);
  putc(10, f);

  /* Image size and maximum gray scale value */
  fprintf(f, "%i %i", width, height);  putc(10, f);
  fprintf(f, "%i", maxgray);           putc(10,f);

  /* No problems at all. */
  return( (YESNO) OK );
}


/*-------------------------------------------------------------------*/


YESNO WritePgmLine(FILE* f, BYTE* data, int x)
{
  int   i;

  for(i=1; i<=x; i++)
    {
    putc(*data++, f);
    }
  return( (YESNO) OK );
}


/*======================= P P M - I N P U T ===========================*/


YESNO ReadPpmHeader(FILE* f, int* width, int* height, int* maxvalue)
{
  int  ch, Id1, Id2;

  Id1 = getc(f);
  Id2 = getc(f);
  if( Id1 != 'P' || Id2 != '6' )
    {
    ErrorMessage("PPM-identification flag not found. (P6 was expected).\n");
    return(NO);
    }

  /* Read the image size and maximum gray scale value */
  fscanf(f, "%i %i", width, height);
  fscanf(f, "%i", maxvalue);

  /* Check validity */
  if( (*maxvalue) != 255 )
    {
    ErrorMessage("This module supports 8-bit color components only."
           "The current image file has %i levels.\n", (*maxvalue)+1 );
    return(NO);
    }

  /* Skip ONE white space, see 'man pgm' */
  ch = getc(f);

  return( (YESNO) OK );
}


/*-------------------------------------------------------------------*/


YESNO ReadPpmLine(FILE*  f,
                       BYTE*  dataR,
                       BYTE*  dataG,
                       BYTE*  dataB,
                       int    x)
{
  int   i;

  for(i=1; i<=x; i++)
    {
    (*dataR++) = getc(f);
    (*dataG++) = getc(f);
    (*dataB++) = getc(f);
    }
  return( (YESNO) OK );
}


/*-------------------------------------------------------------------*/


YESNO ReadPpmPixel(FILE* f, BYTE* R, BYTE* G, BYTE* B)
{
  (*R) = getc(f);
  (*G) = getc(f);
  (*B) = getc(f);
  return( (YESNO) OK );
}


/*======================= P P M - O U P U T ===========================*/


YESNO WritePpmHeader(FILE* f, int width, int height, int maxvalue)
{
  /* P6 identification flad */
  putc('P', f);
  putc('6', f);
  putc(10, f);

  /* Image size and maximum color value */
  fprintf(f, "%i %i", width, height);  putc(10, f);
  fprintf(f, "%i", maxvalue);           putc(10,f);

  /* No problems at all. */
  return( (YESNO) OK );
}


/*-------------------------------------------------------------------*/


YESNO WritePpmLine(FILE*  f,
                        BYTE*  dataR,
                        BYTE*  dataG,
                        BYTE*  dataB,
                        int    x)
{
  int   i;

  for(i=1; i<=x; i++)
    {
    putc(*dataR++, f);
    putc(*dataG++, f);
    putc(*dataB++, f);
    }
  return( (YESNO) OK );
}


/*-------------------------------------------------------------------*/


YESNO WritePpmPixel(FILE* f, BYTE R, BYTE G, BYTE B)
{
  putc(R, f);
  putc(G, f);
  putc(B, f);
  return( (YESNO) OK );
}


/*======================= P B M - I N P U T ===========================*/


YESNO ReadPbmHeader(FILE* f, int* width, int* height)
{
  int  ch, Id1, Id2;

  Id1 = getc(f);
  Id2 = getc(f);
  if( Id1 != 'P' || Id2 != '4' )
    {
    ErrorMessage("PBM-identification flag not found. (P4 was expected).\n");
    return(NO);
    }

  /* Read the image size */
  fscanf(f, "%i %i", width, height);

  /* Skip white space; only one character can appear! */
  ch = getc(f);

  return( (YESNO) OK );
}


/*-------------------------------------------------------------------*/


static int InputPBMBit(FILE* f, int ControlCode)
{
  static int  CurrentBitMask = 0;
  static int  BitQueue       = 0;
  int bit;

  /* Flush bit buffer; only at the end of line situations. */
  if( ControlCode==FLUSHNOW )
    {
    CurrentBitMask = 0;
    return(OK);
    }

  if( CurrentBitMask == 0 )
    {
    BitQueue = getc(f);
    if( BitQueue == EOF )
      {
      ErrorMessage("\nUnexpected EOF in input file.\n");
      return( EOF );
      }
    CurrentBitMask = 0x80;
    }
  bit = ((BitQueue & CurrentBitMask) != 0);
  CurrentBitMask >>= 1;
  return( bit );
}


/*-------------------------------------------------------------------*/


int FlushInputBitBuffer(FILE* f)
{
  return( InputPBMBit(f,FLUSHNOW) );
}


/*-------------------------------------------------------------------*/


YESNO ReadPbmLine(FILE* f, BYTE* data, int  x)
{
  int   i;

  for(i=1; i<=x; i++)
    {
    (*data++) = (BYTE) InputPBMBit(f,OK);
    }
  FlushInputBitBuffer(f);
  return( (YESNO) OK );
}


/*-------------------------------------------------------------------*/


YESNO ReadPbmBitLine(FILE* f, BYTE* data, int  x)
{
  int   i;

  for(i=1; i<=x; i++)
    {
    (*data++) = getc(f);
    }
  return( (YESNO) OK );
}


/*======================= P B M - O U P U T ===========================*/


YESNO WritePbmHeader(FILE* f, int width, int height)
{
  /* P4 identification flad */
  putc('P', f);
  putc('4', f);
  putc(10, f);

  /* Image size */
  fprintf(f, "%i %i", width, height);  putc(10, f);

  /* No problems at all. */
  return( (YESNO) OK );
}


/*-------------------------------------------------------------------*/


static void OutputPBMBit(FILE* f, int bit)
{
  static int  CurrentBitMask = 0x80;
  static int  BitQueue    = 0;

  /* Flush bit buffer; only at the end of line situations. */
  if( bit==FLUSHNOW )
    {
    putc(BitQueue, f);
    CurrentBitMask = 0x80;
    BitQueue = 0;
    return;
    }

  /* Normal bit buffering */
  assert( bit == 0 || bit == 1 );
  if (CurrentBitMask == 0)
    {
    putc(BitQueue, f);
    BitQueue = 0;
    CurrentBitMask = 0x80;
    }
  if (bit) BitQueue |= CurrentBitMask;
  CurrentBitMask >>= 1;

}


/*-------------------------------------------------------------------*/


YESNO FlushOutputBitBuffer(FILE* f)
{
  OutputPBMBit(f, FLUSHNOW);
  return((YESNO) OK);
}


/*-------------------------------------------------------------------*/


YESNO WritePbmLine(FILE* f, BYTE* data, int x)
{
  int   i;

  for(i=1; i<=x; i++)
    {
    OutputPBMBit(f, *data++);
    }
  FlushOutputBitBuffer(f);
  return( (YESNO) OK );
}


/*-------------------------------------------------------------------*/


YESNO WritePbmBitLine(FILE* f, BYTE* data, int  x)
{
  int   i;

  for(i=1; i<=x; i++)
    {
    putc(*data++, f);
    }
  return( (YESNO) OK );
}


/*================== P N M - General PPM/PGM - I N P U T ==============*/


YESNO ReadPnmHeader(FILE* f, int* width, int* height, int* maxvalue,
                    FILETYPE* ft)
{
  int  ch, Id1, Id2;

  Id1 = getc(f);
  Id2 = getc(f);

  if( Id1 != 'P' ||
      (Id2 != '2' && Id2 != '4' && Id2 !='5' && Id2 !='6') )
    {
    ErrorMessage("PPM-identification flag not supported. (P2, P4, P5 or P6 expected).\n");
    return(NO);
    }

  /* Read the image size */
  fscanf(f, "%i %i", width, height);

  switch(Id2)
    {
    case '2': (*ft) = PGMASCII; break;
    case '4': (*ft) = PBM; break;
    case '5': (*ft) = PGM; break;
    case '6': (*ft) = PPM; break;
    default:  (*ft) = PGM; break; /* Error??? */
    }

  if( Id2=='4' ) /* *ft == PBM */
    {
    *maxvalue = 1;
    }
  else
    {
    /* Read the maximum gray scale value */
    fscanf(f, "%i", maxvalue);

    if( *ft == PGM )
      {
      if( *maxvalue == 255 )
        {
        /* OK, filetype is PGM */
        *ft = PGM; /* Redundant but who cares. */
        }
      else
        {
        if( 256 <= *maxvalue && *maxvalue <= 65535 )
          {
          *ft = PGMEXTENDED;
          }
        else
          {
          ErrorMessage("This module supports from 8 to 16 bit gray components for PGM."
                   "The current image file has %i levels.\n", (*maxvalue)+1 );
          return(NO);
          }
        }
      }
    else
      {
      /* Check validity */
      if( (*maxvalue) != 255 )
        {
        ErrorMessage("This module supports 8-bit color components only."
               "The current image file has %i levels.\n", (*maxvalue)+1 );
        return(NO);
        }
      }
    }

  /* Skip ONE white space. See 'man pgm'. */
  ch = getc(f); 

  return( (YESNO) OK );
}


/*-------------------------------------------------------------------*/


YESNO ReadPnmLine(FILE*     f,
                  FILETYPE  filetype,
                  BYTE*     dataR,
                  BYTE*     dataG,
                  BYTE*     dataB,
                  int       x)
{
  int   i;
  char  s[10];

  for( i=1; i<=x; i++ )
    {
    switch(filetype)
      {
      case PBM: (*dataR++) = (BYTE) InputPBMBit(f,OK); break;
      case PGM: (*dataR++) = getc(f); break;
      case PPM: (*dataR++) = getc(f);
                (*dataG++) = getc(f);
                (*dataB++) = getc(f); break;
      case PGMASCII: fscanf(f, "%s", s);
                (*dataR++) = atoi(s); break;
      case PGMEXTENDED:
                (*dataR++) = getc(f);
                (*dataG++) = getc(f); break;
      default:  break;
      }
    }
  FlushInputBitBuffer(f); /* For PGM and PPM also??? */
  return( (YESNO) OK );
}


/*-------------------------------------------------------------------*/


YESNO ReadPnmPixel(FILE* f, FILETYPE filetype, BYTE* R, BYTE* G, BYTE* B)
{
  char  s[10];

  switch(filetype)
    {
    case PBM: (*R) = (BYTE) InputPBMBit(f,OK); break;
    case PGM: (*R) = getc(f); break;
    case PPM: (*R) = getc(f);
              (*G) = getc(f);
              (*B) = getc(f); break;
    case PGMASCII: fscanf(f, "%s", s);
              (*R) = atoi(s); break;
    case PGMEXTENDED:
              (*R) = getc(f);
              (*G) = getc(f); break;
    default:  break;
    }
  return( OK );
}


/*================= P N M - General PPM/PGM - O U T P U T =============*/


YESNO WritePnmHeader(FILE* f, int width, int height, int maxvalue,
                     FILETYPE ft)
{
  /* P4, P5 or P6 identification flag */

  switch(ft)
    {
    case PBM: putc('P',f); putc('4',f); putc(10,f); break;
    case PGM: putc('P',f); putc('5',f); putc(10,f); break;
    case PPM: putc('P',f); putc('6',f); putc(10,f); break;
    case PGMASCII:    putc('P',f); putc('2',f); putc(10,f); break;
    case PGMEXTENDED: putc('P',f); putc('5',f); putc(10,f); break;
    default:  break;
    }

  /* Image size and maximum color value */
  fprintf(f, "%i %i", width, height);  putc(10,f);
  if( ft != PBM)
    {
    fprintf(f, "%i", maxvalue); putc(10,f);
    }

  /* No problems at all. */
  return( (YESNO) OK );
}


/*-------------------------------------------------------------------*/


YESNO WritePnmLine(FILE*     f,
                   FILETYPE  filetype,
                   BYTE*     dataR,
                   BYTE*     dataG,
                   BYTE*     dataB,
                   int       x)
{
  int i;

  for(i=1; i<=x; i++)
    {
    switch(filetype)
      {
      case PBM: OutputPBMBit(f, *dataR++);  break;
      case PGM: putc(*dataR++, f); break;
      case PPM: putc(*dataR++, f);
                putc(*dataG++, f);
                putc(*dataB++, f); break;
      case PGMASCII:
                fprintf(f, "%i ", *dataR++); break;
      case PGMEXTENDED:
                putc(*dataR++, f);
                putc(*dataG++, f); break;
      default:  break;
      }
    }
  return( (YESNO) OK );
}


/*-------------------------------------------------------------------*/


YESNO WritePnmPixel(FILE* f, FILETYPE filetype, BYTE R, BYTE G, BYTE B)
{
  switch(filetype)
    {
    case PBM: OutputPBMBit(f, R);                 break;
    case PGM: putc(R, f);                         break;
    case PPM: putc(R, f); putc(G, f); putc(B, f); break;
    case PGMASCII:
              fprintf(f, "%i ", R); break;
    case PGMEXTENDED:
              putc(R, f); putc(G, f); break;
    default:  break;
    }
  return( (YESNO) OK );
}


/*===================== File type determination  ====================*/


FILETYPE AnalyzeFileType(char* FileName)
{
  FILE*    f;
  FILETYPE ft;
  int      Id1,Id2;
  int      width, height, maxvalue;

  if( !ExistFile(FileName) )
    {
    return(UNKNOWN);
    }

  f = FileOpen(FileName, INPUT, NO);
  Id1 = getc(f);
  Id2 = getc(f);

  if( Id1!='P' )
    {
    ft = UNKNOWN;
    }
  else
    {
    switch(Id2)
      {
      case '2': ft = PGMASCII; break;
      case '4': ft = PBM;      break;
      case '5': fscanf(f, "%i %i %i", &width, &height, &maxvalue);
                if( maxvalue == 255 )
                  {
                  ft = PGM;
                  }
                else
                  {
                  if( 256 <= maxvalue && maxvalue <= 65535 )
                    {
                    ft = PGMEXTENDED;
                    }
                  else
                    {
                    ErrorMessage("This module supports from 8 to 16 bit gray components for PGM."
                             "The current image file has %i levels.\n", (maxvalue)+1 );
                    return(NO);
                    }
                  }
                break;
      case '6': ft = PPM;      break;
      default:  ft = UNKNOWN;  break;
      }
    }

  fclose(f);

  return( ft );
}


/*-------------------------------------------------------------------*/


char* FileFormatName(FILETYPE filetype)
{
  switch(filetype)
    {
    case PBM: return(FormatNamePBM);
    case PGMASCII:
    case PGMEXTENDED:
    case PGM: return(FormatNamePGM);
    case PPM: return(FormatNamePPM);
    default:  return(FormatNamePGM);
    }
  return(FormatNamePGM);
}


/* ---------------------- Filetype determination ------------------------- */


FILETYPE DetermineFileFormat(char* FileName)
/* Determines the type of the input file (PBM/PGM/PPM) and returns
   the result. If no extension given in the file name, the routine
   tries each possibity and chooses the one that matches. */
{
  char      s[128];
  FILETYPE  filetype;

  /* Without extension */
  strcpy(s, FileName);
  filetype = AnalyzeFileType(s);
  if(filetype!=UNKNOWN)
    {
    return(filetype);
    }

  /* Try PBM-extension */
  strcpy(s, FileName);
  CheckFileName(s, FileFormatName(PBM));
  filetype = AnalyzeFileType(s);
  if(filetype!=UNKNOWN)
    {
    return(filetype);
    }

  /* Try PGM-extension */
  strcpy(s, FileName);
  CheckFileName(s, FileFormatName(PGM));
  filetype = AnalyzeFileType(s);
  if(filetype!=UNKNOWN)
    {
    return(filetype);
    }

  /* Try PPM-extension */
  strcpy(s, FileName);
  CheckFileName(s, FileFormatName(PPM));
  filetype = AnalyzeFileType(s);
  if(filetype!=UNKNOWN)
    {
    return(filetype);
    }

  /* No luck */
  printf("--- %s ---\n",FileName);
  return(UNKNOWN);
}
