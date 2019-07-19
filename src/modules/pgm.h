#if ! defined(__PGM_H)
#define __PGM_H

#include "owntypes.h"

/* ----------------------- Type definitions ------------------------------ */

typedef  enum { UNKNOWN=0, PBM=1, PGM=2, PPM=3, PGMASCII=4, PGMEXTENDED=5 }
         FILETYPE;
typedef  enum { Gray=0, Red=0, Green=1, Blue=2 } COLORTYPE;
#define  FormatNamePBM    "pbm"
#define  FormatNamePGM    "pgm"
#define  FormatNamePPM    "ppm"

/* ---------------------- Filetype determination ------------------------- */

char*    FileFormatName(FILETYPE filetype);
FILETYPE AnalyzeFileType(char* FileName);
FILETYPE DetermineFileFormat(char* FileName);

/* -------------------- PGM - Gray-scale images -------------------------- */

YESNO ReadPgmHeader(FILE*  f,
                    int*   width,
                    int*   height,
                    int*   maxgray);
YESNO WritePgmHeader(FILE* f,
                     int   width,
                     int   height,
                     int   maxgray);
YESNO ReadASCIIPgmHeader(FILE*  f,
                         int*   width,
                         int*   height,
                         int*   maxgray);
YESNO ReadPgmLine(FILE* f,  unsigned char *data,  int  x_size);
YESNO WritePgmLine(FILE* f, unsigned char *data,  int  x_size);
YESNO ReadASCIIPgmLine(FILE* f,  unsigned char *data,  int  x_size);

/* ------------------- PBM - black & white images ------------------------ */

YESNO ReadPbmHeader(FILE*  f,
                    int*   width,
                    int*   height);
YESNO WritePbmHeader(FILE* f,
                     int   width,
                     int   height);
YESNO ReadPbmLine(FILE* f,  unsigned char *data,  int  x_size);
YESNO WritePbmLine(FILE* f, unsigned char *data,  int  x_size);
int   FlushInputBitBuffer(FILE* f);
YESNO FlushOutputBitBuffer(FILE* f);
YESNO ReadPbmBitLine(FILE* f, unsigned char *data, int  x_size);
YESNO WritePbmBitLine(FILE* f, unsigned char *data, int  x_size);

/* ----------------------- PPM - Color images ---------------------------- */

YESNO ReadPpmHeader(FILE*  f,
                    int*   width,
                    int*   height,
                    int*   maxvalue);
YESNO WritePpmHeader(FILE* f,
                     int   width,
                     int   height,
                     int   maxvalue);
YESNO ReadPpmLine(FILE*          f,
                  unsigned char *dataR,
                  unsigned char *dataG,
                  unsigned char *dataB,
                  int           x_size);
YESNO WritePpmLine(FILE*          f,
                   unsigned char *dataR,
                   unsigned char *dataG,
                   unsigned char *dataB,
                   int  x_size);
YESNO ReadPpmPixel(FILE*          f,
                   unsigned char *R,
                   unsigned char *G,
                   unsigned char *B);
YESNO WritePpmPixel(FILE*        f,
                    unsigned char R,
                    unsigned char G,
                    unsigned char B);

/* -------------------------- PNM - General ------------------------------ */

YESNO ReadPnmHeader(FILE*      f,
                    int*       width,
                    int*       height,
                    int*       maxvalue,
                    FILETYPE*  filetype);
YESNO WritePnmHeader(FILE*     f,
                     int       width,
                     int       height,
                     int       maxvalue,
                     FILETYPE  filetype);
YESNO ReadPnmLine(FILE*          f,
                  FILETYPE       filetype,
                  unsigned char *dataR,
                  unsigned char *dataG,
                  unsigned char *dataB,
                  int            x_size);
YESNO WritePnmLine(FILE*          f,
                   FILETYPE       filetype,
                   unsigned char *dataR,
                   unsigned char *dataG,
                   unsigned char *dataB,
                   int            x_size);
YESNO ReadPnmPixel(FILE*          f,
                   FILETYPE       filetype,
                   unsigned char  *R,
                   unsigned char  *G,
                   unsigned char  *B);
YESNO WritePnmPixel(FILE*        f,
                    FILETYPE      filetype,
                    unsigned char R,
                    unsigned char G,
                    unsigned char B);


#endif /* __PGM_H */
