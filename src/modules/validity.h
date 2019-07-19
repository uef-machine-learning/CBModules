#if ! defined(__CVI_H)
#define __CVI_H

/*--------------------  Basic type definitions -----------------------*/

#include "cb.h"



// the function for getting mean vector of total dataset
VECTORTYPE MeanOfTs(TRAININGSET *pTS, CODEBOOK *pCB);
// SSW and SSB for sum-of-squares
double ValidityMSE(TRAININGSET *pTS, CODEBOOK *pCB, PARTITIONING *pP);
double ValiditySSB(TRAININGSET *pTS, CODEBOOK *pCB, PARTITIONING *pP);

// commonly used index 
double ValidityDunn(TRAININGSET *pTS, CODEBOOK *pCB, PARTITIONING *pP);
double ValidityDBI(TRAININGSET *pTS, CODEBOOK *pCB, PARTITIONING *pP);
double ValidityFratio(TRAININGSET *pTS, CODEBOOK *pCB, PARTITIONING *pP);
//another implementation of F-ratio i.e. f-test from ISMO
double FTest(TRAININGSET *pTS, CODEBOOK *pCB, PARTITIONING *pP);
double Silhouette(TRAININGSET *pTS, CODEBOOK *pCB, PARTITIONING *pP);
double SilhouetteCI(TRAININGSET *pTS, CODEBOOK *pCB, PARTITIONING *pP);
double RMSSTD(TRAININGSET *pTS, CODEBOOK *pCB, PARTITIONING *pP);
double RSquare(TRAININGSET *pTS, CODEBOOK *pCB, PARTITIONING *pP);
double SD_Scatter(TRAININGSET *pTS, CODEBOOK *pCB, PARTITIONING *pP);
double SD_Dis(TRAININGSET *pTS, CODEBOOK *pCB, PARTITIONING *pP);

double Scatter(TRAININGSET *pTS, CODEBOOK *pCB, PARTITIONING *pP);
double variance_cluster(TRAININGSET *pTS, CODEBOOK *pCB, PARTITIONING *pP);
double S_Dbw(TRAININGSET *pTS, CODEBOOK *pCB, PARTITIONING *pP);
// Xie-Beni in partition-based clustering
double Xie_Beni(TRAININGSET *pTS, CODEBOOK *pCB, PARTITIONING *pP);

// for fuzzy clustering 
//double Xie_Beni_Fuzzy(TRAININGSET *pTS, CODEBOOK *pCB, FuzzyPartitioning *pP);

// Baysian Information Criteria in partition-based clustering
double ValidityBIC(TRAININGSET *pTS, CODEBOOK *pCB, PARTITIONING *pP);

//External index
double RandIndex(TRAININGSET *pTS, PARTITIONING *pP, PARTITIONING *pPtruth);
double CorrectedRI(TRAININGSET *pTS, PARTITIONING *pP, PARTITIONING *pPtruth);
//Mutual information function
double NMI(TRAININGSET *pTS, CODEBOOK *pCB, PARTITIONING *pP, PARTITIONING *pPtruth);


/* ----------------------------------------------------------------- */

#endif /* __CVI_H */

