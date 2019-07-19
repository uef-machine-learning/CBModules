
#ifndef INITS_H
#define INITS_H
void LuxburgInitialCentroids(TRAININGSET *pTS,CODEBOOK *pCB);
void ProjectionInitialCentroids(TRAININGSET *pTS,CODEBOOK *pCB,int variant);

void DensityInitialCentroids(TRAININGSET *pTS,CODEBOOK *pCB,int variant);
double* calcDeltas( TRAININGSET *pTS, double* density, int* densityOrder );
double* calcDeltasPerDim( TRAININGSET *pTS, double* density, int* densityOrder );
double* calcDeltasFromWindow( TRAININGSET *pTS, double* density, int* densityOrder);

void hierarchicalSplittingCodebook(TRAININGSET *TS,CODEBOOK *pCBnew,int kmIters);
int SplitCluster(TRAININGSET *TS,CODEBOOK *CBleft,CODEBOOK *CBright,TRAININGSET *TSleft,TRAININGSET *TSright,int kmIters);
float dotprod(int *a,int *b,int size);
int *vecminus(int *a,int *b,int size);
void RandomPartitionInitialCentroids(TRAININGSET *pTS,CODEBOOK *pCB,int NumTries);
void BradleyInitialCentroids(TRAININGSET *pTS,CODEBOOK *pCB,int J);
void SortingHeuristic(TRAININGSET *pTS,CODEBOOK *pCB,int J);
void MaxMinInitialCentroids(TRAININGSET *pTS,CODEBOOK *pCB,int InitialPoint);
void CreateMaxMinCodebook(TRAININGSET *pTS,CODEBOOK *sourceCB,CODEBOOK *targetCB,int k);
void SelectRandomWeightedRepresentatives(TRAININGSET *pTS,CODEBOOK *pCB, int numLocalTries);
llong *GetLongVector(int dim);
void SelectRandomRepresentativesbyMarko(TRAININGSET *pTS,CODEBOOK *pCB);
void RandomCodebook(TRAININGSET *pTS,CODEBOOK *pCB);
void SelectRandomRepresentatives2(TRAININGSET *pTS,CODEBOOK *pCB);
void SelectRandomRepresentatives(TRAININGSET *pTS,CODEBOOK *pCB);
/*double* knnGraphDensity(TRAININGSET *pTS, int knnK);*/
void  decreaseOverlapInit(TRAININGSET *pTS, CODEBOOK *pCB, int variant);


#endif
