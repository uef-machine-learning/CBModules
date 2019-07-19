#if ! defined(__DIP_H)
#define __DIP_H

#include "image.h"

/* ----------------------------------------------------------------- */

void  CalculateHistogram(IMAGE* Image, long Histogram[3][256]);
void  EqualizeHistogram(IMAGE* Image, long Histogram[3][256], int Mapping[3][256]);
void  PrintHistogram(IMAGE* Image, long Histogram[3][256]);
int   SeparateColor(int r, int g, int b, int color);
void  CalculatePseudoColor(int v, int* r, int* g, int* b, int mapping);
int   QuantizePixel(int v, int bits);
int   ThresholdPixel(int v, int threshold);
int   GetBit(int v, int bit);
BYTE  GenerateSaltPepperNoise(int v, float p);
BYTE  GenerateAdditiveNoise(int v, float p);
float BlockMean(IMAGE* Image, int x, int y, int color);
float BlockVariance(IMAGE* Image, int x, int y, int color);
int   BlockMedian(IMAGE* Image, int x, int y, int color);
long  MaskOperation(IMAGE* Image, int x, int y, int color, int mask);

/* ----------------------------------------------------------------- */

#endif /* __DIP_H */

