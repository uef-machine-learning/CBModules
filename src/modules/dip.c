/*-------------------------------------------------------------------*/
/* DIP.C          Pasi Fr„nti                                        */
/*                                                                   */
/* Digital Image Processing - basic methods                          */
/*                                                                   */
/*-------------------------------------------------------------------*/

#define ProgName        "DIP"
#define VersionNumber   "Version 0.04"
#define LastUpdated     "12.4.99"

/* ----------------------------------------------------------------- */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "image.h"
#include "interfc.h"
#include "owntypes.h"
#include "random.h"

/*-----------------------------  M i s c  ---------------------------*/

#if defined(max) || defined(min)
#undef max
#undef min
#endif
#define  max(a,b) ((a) > (b) ? (a) : (b))
#define  min(a,b) ((a) < (b) ? (a) : (b))
#define  ScaleBetween(a,b,c)  max(a, min(b, c))
#define  Dice(p)			  (frand()<p ? YES : NO )
#define  round(a)             ((int) (a+0.5))


/* ===================== Histogram operations ====================== */


void CalculateHistogram(IMAGE* Image, long Histogram[3][256])
{
  int  i,x,y,v,color;

  for(color=0; color<ColorBands(Image); color++)
     for(i=0; i<=Image->MaxValue; i++)
		Histogram[color][i] = 0;

  for(y=1; y<=Image->ImageSizeY; y++)
    {
	ScrollBuffer(Image);
	for(x=1; x<=Image->ImageSizeX; x++)
	   for(color=0; color<ColorBands(Image); color++)
		  {
          v = GetBufferPixel(Image, x, 0, color);
		  Histogram[color][v]++;
		  }
    }
  RewindImage(Image);
}


/*-------------------------------------------------------------------*/


void EqualizeHistogram(IMAGE* Image, long Histogram[3][256], int Mapping[3][256])
{
  int	 color,i;
  double new;
  long   cumulative;

  for(color=0; color<ColorBands(Image); color++)
     for(i=0,cumulative=0; i<=Image->MaxValue; i++)
       {
       cumulative += Histogram[color][i];
       new = (256.0*cumulative)/ImageSize(Image) - 1;
       Mapping[color][i] = max( 0, round(new) );
	   }
}


/*-------------------------------------------------------------------*/


void PrintHistogram(IMAGE* Image, long Histogram[3][256])
{
  int i,color;

  for(i=0; i<=Image->MaxValue; i++)
	 {
     PrintMessage("%3i: ",i);
	 for(color=0; color<ColorBands(Image); color++)
        PrintMessage("%5li", Histogram[color][i]);
     PrintMessage("\n");
	 }
}


/* ======================= Color conversions ======================= */


int SeparateColor(int r, int g, int b, int color)
{
  int    v;
  double x;

  switch(color)  /* R-G-B  Y-U-V  H-S-I */
	 {
	 case 1: v = r; break;
	 case 2: v = g; break;
	 case 3: v = b; break;
	 case 4: v = round(0.6*r + 0.6*g + 0.1*b); break;
	 case 5: v = round( (0.9*b - 0.3*r - 0.6*g + 229.5) / 1.8 ); break;
	 case 6: v = round( (0.7*r - 0.6*g - 0.1*b + 178.5) / 1.4 ); break;
     case 7: x = (r-g)*(r-g)+(r-b)*(g-b);
             x = x>0 ? 0.5*(r-g+r-b) / sqrt(x) : 1;
             x = b>g ? 1-acos(x)/6.28318 : acos(x)/6.28318;
             /* x = x+0.5; if(x>1) x-=1; */
             v = round(255*x); break;
     case 8: v = (r+b+g);
             x = v>0 ? 3.0*min(min(r,g),b)/v : 1;
             v = round(255 * (1-x)); break;
     case 9: v = round((r+g+b)/3); break;
	 default: v = 0; break;
	 }
/* PrintMessage("(%3i,%3i,%3i) -> %3i \n",r,g,b,v); */
  return( ScaleBetween(0,v,255) );
}


/*-------------------------------------------------------------------*/


void CalculatePseudoColor(int v, int* r, int* g, int* b, int mapping)
{
  switch( mapping )
    {
	case 1: (*r) = 3 * v;
			(*g) = 3 * (v-85);
			(*b) = 3 * (v-171);
			break;
	case 2: (*r) = 3 * (v-171);
			(*g) = v;
			(*b) = 63;
			break;
	case 3: (*r) = 1.5 * (v-85);
			(*g) = 63;
			(*b) = 255 - 1.5*v;
			break;
    }
  (*r) = ScaleBetween(0, *r, 255);
  (*g) = ScaleBetween(0, *g, 255);
  (*b) = ScaleBetween(0, *b, 255);
}


/* ===================== Pixelwise operations ====================== */


int QuantizePixel(int v, int bits)
{
  int UpperBits[9]	= { 0x00, 0x80, 0xc0, 0xe0, 0xf0, 0xf8, 0xfc, 0xfe, 0xff };
  int BreakPoint[9] = { 0x80, 0x40, 0x20, 0x10, 0x08, 0x04, 0x02, 0x01, 0x00 };

  if( bits==1 )  return( v>127 ? 255 : 0 );
  else           return( (v & UpperBits[bits]) + BreakPoint[bits] );
}


/*-------------------------------------------------------------------*/


int ThresholdPixel(int v, int threshold)
{
  return( v<=threshold ? 1 : 0 );
}


/*-------------------------------------------------------------------*/


BYTE GenerateSaltPepperNoise(int v, float p)
{
  if(Dice(p)) return(Dice(0.5) ? 255 : 0);
  return(v);
}


/*-------------------------------------------------------------------*/


BYTE GenerateAdditiveNoise(int v, float p)
{
  /* Not yet properly implemented */
  v = v + (Dice(0.5) ? 1 : -1) * frand()*p*255;
  v = ScaleBetween(0, v, 255);
  return(v);
}


/*-------------------------------------------------------------------*/


float BlockMean(IMAGE* Image, int x, int y, int color)
{
  int p1 = GetBufferPixel(Image, x-1, y-1, color);
  int p2 = GetBufferPixel(Image, x,   y-1, color);
  int p3 = GetBufferPixel(Image, x+1, y-1, color);
  int p4 = GetBufferPixel(Image, x-1, y,   color);
  int p5 = GetBufferPixel(Image, x,   y,   color);
  int p6 = GetBufferPixel(Image, x+1, y,   color);
  int p7 = GetBufferPixel(Image, x-1, y+1, color);
  int p8 = GetBufferPixel(Image, x,   y+1, color);
  int p9 = GetBufferPixel(Image, x+1, y+1, color);

  return( (p1+p2+p3+p4+p5+p6+p7+p8+p9) / 9.0 );
}


/*-------------------------------------------------------------------*/


float BlockVariance(IMAGE* Image, int x, int y, int color)
{
  long  sqsum;
  float mean, sqmean;

  int p1 = GetBufferPixel(Image, x-1, y-1, color);
  int p2 = GetBufferPixel(Image, x,   y-1, color);
  int p3 = GetBufferPixel(Image, x+1, y-1, color);
  int p4 = GetBufferPixel(Image, x-1, y,   color);
  int p5 = GetBufferPixel(Image, x,   y,   color);
  int p6 = GetBufferPixel(Image, x+1, y,   color);
  int p7 = GetBufferPixel(Image, x-1, y+1, color);
  int p8 = GetBufferPixel(Image, x,   y+1, color);
  int p9 = GetBufferPixel(Image, x+1, y+1, color);

  mean = (p1+p2+p3+p4+p5+p6+p7+p8+p9) / 9.0;
  sqsum = (long) p1*p1 + (long) p2*p2 + (long) p3*p3 +
		  (long) p4*p4 + (long) p5*p5 + (long) p6*p6 +
		  (long) p7*p7 + (long) p8*p8 + (long) p9*p9;
  sqmean = sqsum / 9.0;
  return( sqmean - mean*mean );
}


/*-------------------------------------------------------------------*/


int BlockMedian(IMAGE* Image, int x, int y, int color)
{
  int  p[10];
  int  i,j,temp;

  p[1] = GetBufferPixel(Image, x-1, y-1, color);
  p[2] = GetBufferPixel(Image, x,   y-1, color);
  p[3] = GetBufferPixel(Image, x+1, y-1, color);
  p[4] = GetBufferPixel(Image, x-1, y,   color);
  p[5] = GetBufferPixel(Image, x,   y,   color);
  p[6] = GetBufferPixel(Image, x+1, y,   color);
  p[7] = GetBufferPixel(Image, x-1, y+1, color);
  p[8] = GetBufferPixel(Image, x,   y+1, color);
  p[9] = GetBufferPixel(Image, x+1, y+1, color);

  /* sort values */
  for(i=1; i<=8; i++)
	 for(j=i+1; j<=9; j++)
        if(p[i]>p[j]) { temp=p[i]; p[i]=p[j]; p[j]=temp; }

  return( p[5] );
}


/*-------------------------------------------------------------------*/


int Weight[5][9] = {{ 0, 0, 0,   0, 0, 0,   0, 0, 0 },
                    { 1, 1, 1,   1, 1, 1,   1, 1, 1 },
                    {-1,-1,-1,  -1, 9,-1,  -1,-1,-1 },
                    {-1,-2,-1,   0, 0, 0,  +1,+2,+1 },
                    {-1, 0,+1,  -2, 0,+2,  -1, 0,+1 }};


/*-------------------------------------------------------------------*/


long MaskOperation(IMAGE* Image, int x, int y, int color, int mask)
{
  long sum;

  sum = GetBufferPixel(Image, x-1, y-1, color) * Weight[mask][0] +
		GetBufferPixel(Image, x,   y-1, color) * Weight[mask][1] +
		GetBufferPixel(Image, x+1, y-1, color) * Weight[mask][2] +
		GetBufferPixel(Image, x-1, y,	color) * Weight[mask][3] +
		GetBufferPixel(Image, x,   y,	color) * Weight[mask][4] +
		GetBufferPixel(Image, x+1, y,	color) * Weight[mask][5] +
		GetBufferPixel(Image, x-1, y+1, color) * Weight[mask][6] +
		GetBufferPixel(Image, x,   y+1, color) * Weight[mask][7] +
		GetBufferPixel(Image, x+1, y+1, color) * Weight[mask][8];

  return(sum);
}


