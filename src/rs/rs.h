#if ! defined(__RS_H)
#define __RS_H

int PerformRS(TRAININGSET *pTS, CODEBOOK *pCB, PARTITIONING *pP, int iter, 
    int kmIter, int deterministic, int travellerSearch, int quietLevel, 
    int useInitialCB, int monitoring);

void SelectRandomRepresentatives(TRAININGSET *pTS, CODEBOOK *pCB);
void SelectRandomRepresentatives2(TRAININGSET *pTS, CODEBOOK *pCB);
void SelectRandomRepresentativesbyMarko(TRAININGSET *pTS, CODEBOOK *pCB);
void SelectRandomWeightedRepresentatives(TRAININGSET *pTS, CODEBOOK *pCB);
void LuxburgInitialCentroids(TRAININGSET *pTS, CODEBOOK *pCB); 

void CalculateDistances(TRAININGSET *pTS, CODEBOOK *pCB, PARTITIONING *pP,
    llong *distance);

void OptimalRepresentatives(PARTITIONING *pP, TRAININGSET *pTS, CODEBOOK *pCB, 
     int *active, llong *cdist, int *activeCount, int travellerSearch);

void OptimalPartition(CODEBOOK *pCB, TRAININGSET *pTS, PARTITIONING *pP,
    int *active, llong *cdist, int activeCount, llong *distance, int travellerSearch);

char* RSInfo(void);

#endif /* __RS_H */
