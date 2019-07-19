/*
$Revision: 1.2 $
$Date: 2005/06/27 11:18:17 $
$Author: mtuonone $
$Name:  $
$Id: kmeans.h,v 1.2 2005/06/27 11:18:17 mtuonone Exp $
*/

#if ! defined(__KMEANS_H)
#define __KMEANS_H

int PerformKMeans(TRAININGSET *pTS, CODEBOOK *pCB, PARTITIONING *pP,
		  int clus, int repeats, int InitMethod, int InitMethodVar, int quietLevel,
          int useInitial, int MaxIter);

char* KMeansInfo(void);

void InitializeSolutions(TRAININGSET *pTS, CODEBOOK *pCB,
                                PARTITIONING *pP, CODEBOOK *pCBnew,
                                PARTITIONING *pPnew, CODEBOOK *pCBinit,
                                PARTITIONING *pPinit, llong *distanceInit,
                                int clus, int initial);

void GenerateSolution(TRAININGSET *pTS, CODEBOOK *pCBnew,
                             PARTITIONING *pPnew, CODEBOOK *pCBinit, PARTITIONING *pPinit,
			     llong *distance, llong *distanceInit,
			     int InitMethod, int InitMethodVar, int initial);



void KMeansIterate(TRAININGSET *pTS, CODEBOOK *pCBnew,
PARTITIONING *pPnew, llong *distance, int quietLevel, int i, int *iter,
double time, double *error, int initial, int MaxIter);

//Moved from rs.c:
void CalculateDistances(TRAININGSET *pTS, CODEBOOK *pCB, PARTITIONING *pP,
    llong *distance);
void OptimalRepresentatives(PARTITIONING *pP, TRAININGSET *pTS,
    CODEBOOK *pCB, int *active, llong *cdist, int *activeCount, int travellerSearch);
void OptimalPartition(CODEBOOK *pCB, TRAININGSET *pTS, PARTITIONING *pP,
int *active, llong *cdist, int activeCount, llong *distance, int quietLevel);
int BinarySearch(int *arr, int size, int key);

int KMeansDensityPeaks(TRAININGSET *pTS, CODEBOOK *pCB, PARTITIONING *pP,
        int clus, int repeats, int quietLevel);


typedef     struct {
    int densityPeaksFilterMethod; //0=delta, 1=gamma=delta*density
    int densityPeaksDeltaMethod; //0=normal, 1=DDDE experimental
    int densityMethod; //0=knnGraph,1=DDDE
    int QuietLevel; //0=knnGraph,1=DDDE
    float sample; // From 0 to 1
} KmeansOpt;

extern KmeansOpt kmOpt;

#endif /* __KMEANS_H */


