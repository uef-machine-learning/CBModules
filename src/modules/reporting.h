#if ! defined(__REPORTING_H)
#define __REPORTING_H

#define MAXFILENAME     40    /* maximum length of filename */

double SetClock(double* start);

double GetClock(double start);

void PrintHeader(int quietLevel);

void PrintIterationKM(int quietLevel, int i, int iter, double error,
     double time);

void PrintIterationActivity(double time, int iter, int activeCount,
     int CBCount, int quietLevel);

void PrintIterationKMSummary(double KMtime, double KMinitTime);

void PrintIterationRS(int quietLevel, int iter, double error, int ci, 
     double time, int better);

void PrintRepeat(int quietLevel, int repeats, int i, int iter, double error,
     double time, int better);

void PrintXM(int quietLevel, int i, int iter, double error,
     double time, int better);

void PrintFooterKM(int quietLevel, double error, int repeats,
     double totalTime, int totalIter);
	
void PrintFooterXM(int quietLevel, double error, int clusters, int score, int ScoreType,
     double totalTime, int totalIter);

void PrintFooterDBSCAN(int quietLevel, double time, int nCluster, 
     int nNoisePoint);

void PrintFooterRS(int quietLevel, int iter, double error, double time);

int DetermineFileName(char *name);

TRAININGSET CheckParameters(char *TSName, char *CBName, char *PAName, 
    char *InName, int clus, int ow);

int ReadInitialCBorPA(char *InName, int clus, TRAININGSET *pTS, 
    CODEBOOK *pCB, PARTITIONING *pP);

#endif /* __REPORTING_H */

