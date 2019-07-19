/*-------------------------------------------------------------------*/
/* REPORTING.C     Marko Tuononen + Radu Mariescu-Istodor            */
/*                                                                   */
/* Provides reporting (time and error) facilities.                   */
/*                                                                   */
/*-------------------------------------------------------------------*/

#define ProgName       "REPORTING"
#define VersionNumber  "Version 0.25" // PF
#define LastUpdated    "7.9.2016"

/* ----------------------------------------------------------------- */
/* Changelog:                                                        */
/*                                                                   */
/* PF  7.9.16   0.25    Added PrintFooterDBSCAN                      */
/* PF  5.7.16   0.24    Modified time profilation for K-means        */
/* PF 24.6.16   0.23    CI-value to PrintIterationRS (redefined)     */
/* PF 23.6.16   0.22    Unifying; RSL->RS; Minor editing;            */
/* RMI          0.21    Printing functions for XMeans                */
/*                                                                   */
/* ----------------------------------------------------------------- */

#include <sys/time.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "interfc.h"
#include "cb.h"
#include "file.h"
#include "reporting.h"
#include "textfile.h"


/* =================== FUNCTIONS ==================================== */


double SetClock(double* start)
{
  struct timeval tv;
  gettimeofday(&tv, NULL);

  return( *start = (double)tv.tv_sec + ((double)tv.tv_usec /1e6) );
}


/* ------------------------------------------------------------------ */


double GetClock(double start)
{
  struct timeval tv;
  double elapsed;

  gettimeofday(&tv, NULL);
  elapsed = (double)tv.tv_sec + ((double)tv.tv_usec /1e6) - start;

  return elapsed;
}


/* ------------------------------------------------------------------ */


void PrintHeader(int quietLevel)
{
  if ((quietLevel >= 2))
  {
    PrintMessage("\n");
    PrintMessage("Iteration\tMSE\t\tTime\n");
  }
}


/* ------------------------------------------------------------------ */


void PrintIterationKM(int quietLevel, int i, int iter, double error,
double time)
{
if ((quietLevel >= 2) && (i == 0))
  {
  PrintMessage("%d\t\t%f\t%f\n", iter, error, time);
  }
}


/* ------------------------------------------------------------------ */


void PrintIterationActivity(double time, int iter, int activeCount,
int CBCount, int quietLevel)
{
  if (quietLevel >= 3)
    {
    PrintMessage("KM-time(%i)= %f   ", iter, time);
    }
  if (quietLevel >= 4)
    {
    PrintMessage("\t\t  K-means iter=%d\t  \t%d/%d Active centroids\n",
                 iter, activeCount, CBCount);
    }
}


/* ------------------------------------------------------------------ */


void PrintIterationKMSummary(double KMtime, double KMinitTime)
{
PrintMessage("\t\t\t\t\t <%f>\t<onInit %f>\n\n", KMtime, KMinitTime);
}



/* ------------------------------------------------------------------ */


void PrintIterationRS(int quietLevel, int iter, double error, int ci,
double time, int better)
{
  if ((quietLevel >= 3) || ((quietLevel >= 2) && better))
    {
    PrintMessage("%d\t\t%f\t%f", iter, error, time);
    if (better)  PrintMessage(" *");
    if (ci)      PrintMessage(" CI=%i",ci);
    PrintMessage("\n");
    }
}


/* ------------------------------------------------------------------ */


void PrintRepeat(int quietLevel, int repeats, int i, int iter, double error,
double time, int better)
{
  char str[127];

  sprintf(str, "Repeat %d:\t%f\t%f\t(%d iterations)", i+1, error, time, iter);

  if ((quietLevel >= 2) && (repeats > 1))
  {
    if (i == 0)
    {
      PrintMessage("%s *\n\n", str);
    }
    else if (better)
    {
      PrintMessage("%s *\n", str);
    }
    else if (quietLevel >= 3)
    {
      PrintMessage("%s\n", str);
    }
  }
}


/* ------------------------------------------------------------------ */


void PrintXM(int quietLevel, int i, int iter, double error,
double time, int better)
{
  char str[127];

  sprintf(str, "Clusters %d:\t%f\t%f\t(%d iterations)", i, error, time, iter);

  if (quietLevel >= 2)
  {
    if (i == 0)
    {
      PrintMessage("%s *\n\n", str);
    }
    else if (better)
    {
      PrintMessage("%s *\n", str);
    }
    else if (quietLevel >= 3)
    {
      PrintMessage("%s\n", str);
    }
  }
}


/* ------------------------------------------------------------------ */


void PrintFooterKM(int quietLevel, double error, int repeats,
double totalTime, int totalIter)
{
  if (quietLevel >= 2)
  {
    PrintMessage("\nmse = %f  time = %f  %d repeats  %d iterations\n\n", error, totalTime, repeats, totalIter);
  }
  else if (quietLevel == 1)
  {
    PrintMessage("%f\t%f\n", error, totalTime);
  }
}


/* ------------------------------------------------------------------ */


void PrintFooterXM(int quietLevel, double error, int clusters, int score, int ScoreType,
double totalTime, int totalIter)
{
  if (quietLevel >= 2)
  {
    if(ScoreType==1)
		PrintMessage("\nmse = %f  time = %f    %d clusters:(BIC score = %d)\n\n", error, totalTime, clusters, score);
	else
	    PrintMessage("\nmse = %f  time = %f    %d clusters:(AIC score = %d)\n\n", error, totalTime, clusters, score);
  }
  else if (quietLevel == 1)
  {
    PrintMessage("%f\t%f\n", error, totalTime);
  }
}


/* ------------------------------------------------------------------ */


void PrintFooterDBSCAN(int quietLevel, double time, int nCluster,
                       int nNoisePoint)
{
  if(quietLevel >= 1) // always display these
    {
    PrintMessage("\n Number of clusters = %d,  Number of Noisy points = "
    "%d,   Time = %f  \n\n", nCluster, nNoisePoint, time);
    }
}


/* ------------------------------------------------------------------ */


void PrintFooterRS(int quietLevel, int iter, double error, double time)
{
  if (quietLevel >= 2)
  {
    PrintMessage("\nmse = %f  time = %f  %d iterations\n\n", error, time, iter);
  }
  else if (quietLevel == 1)
  {
    PrintMessage("%f\t%f\n", error, time);
  }
}


/* ------------------------------------------------------------------ */


int DetermineFileName(char *name)
{
  char newName[MAXFILENAME], suffix[MAXFILENAME];
  int  i;

  /* Without extension */
  if (ExistFile(name))
  {
    return 1;
  }

  for (i = 0; i < 3; i++)
  {
    if (i == 0)      /* Try TS-file extension */
      strcpy(suffix, FormatNameTS);

    else if (i == 1)  /* Try CB-file extension */
      strcpy(suffix, FormatNameCB);

    else             /* Try PA-file extension */
      strcpy(suffix, FormatNamePA);

    if (strlen(name) < MAXFILENAME-strlen(suffix)-1)
    {
      strcpy(newName, name);
      CheckFileName(newName, suffix);
      if (ExistFile(newName))
      {
        strcpy(name, newName);
        return 1;
      }
    }
  }

  /* No luck this time */
  return 0;
}


/* ------------------------------------------------------------------ */


TRAININGSET CheckParameters(char *TSName, char *CBName, char *PAName,
char *InName, int clus, int ow) {
  TRAININGSET TS;
  char *fileExtension;

  /* input training set doesn't exist */
  if (!ExistFile(TSName))
  {
    ErrorMessage("\nERROR: Input training set doesn't exist: "
        "%s\n\n", TSName);
    ExitProcessing(FATAL_ERROR);
  }

  /* result codebook file exists and we are told not to overwrite */
  if (ExistFile(CBName) && !ow)
  {
    ErrorMessage("\nERROR: Result codebook already exists: "
        "%s\n\n", CBName);
    ExitProcessing(FATAL_ERROR);
  }

  /* result partitioning file exists and we are told not to overwrite */
  if (*PAName && ExistFile(PAName) && !ow)
  {
    ErrorMessage("\nERROR: Result partitioning already exists: "
        "%s\n\n", PAName);
    ExitProcessing(FATAL_ERROR);
  }

  /* initial codebook / partitioning doesn't exist */
  if (*InName && !DetermineFileName(InName))
  {
    ErrorMessage("\nERROR: Initial codebook/partitioning doesn't exist: %s\n\n", InName);
    ExitProcessing(FATAL_ERROR);
  }

  ReadTrainingSet(TSName, &TS);

  /* result codebook cannot contain more vectors than training set */
  if (BookSize(&TS) < clus)
  {
      ErrorMessage("\nERROR: Number of vectors in training set ");
      ErrorMessage("(%d) < number of clusters ", BookSize(&TS));
      ErrorMessage("(%d)!\n\n", clus);
      FreeCodebook(&TS);
      ExitProcessing(FATAL_ERROR);
  }

  return TS;
}


/* ------------------------------------------------------------------ */


int ReadInitialCBorPA(char *InName, int clus, TRAININGSET *pTS,
CODEBOOK *pCB, PARTITIONING *pP)
{
  int useInitial = 0;

  if (*InName)  /* we use initial codebook/partitioning */
  {
    switch (DetermineCBFileType(InName)) {
      case TSFILE: case CBFILE:
        ReadCodebook(InName, pCB);
        useInitial = 1;

        if (BookSize(pCB) != clus)
        {
          ErrorMessage("\nERROR: Number of vectors in initial codebook ");
          ErrorMessage("(%d) <> number of clusters ", BookSize(pCB));
          ErrorMessage("(%d)!\n\n", clus);
          FreeCodebook(pTS);
          FreeCodebook(pCB);
          ExitProcessing(FATAL_ERROR);
        }

        CreateNewPartitioning(pP, pTS, clus);
        break;

      case PAFILE:
        ReadPartitioning(InName, pP, pTS);
        useInitial = 2;

        if (PartitionCount(pP) != clus)
        {
          ErrorMessage("\nERROR: Number of partitions in initial partitioning ");
          ErrorMessage("(%d) <> number of clusters ", PartitionCount(pP));
          ErrorMessage("(%d)!\n\n", clus);
          FreeCodebook(pTS);
          FreePartitioning(pP);
          ExitProcessing(FATAL_ERROR);
        }

        CreateNewCodebook(pCB, clus, pTS);
        break;

      case NOTFOUND:
        ErrorMessage("\nERROR: Type of initial codebook/partitioning file "
            "%s is unidentified!\n\n", InName);
        FreeCodebook(pTS);
        ExitProcessing(FATAL_ERROR);
        break;
    }
  }
  else  /* we don't use initial codebook/partitioning */
  {
    CreateNewCodebook(pCB, clus, pTS);
    CreateNewPartitioning(pP, pTS, clus);
    useInitial = 0;
  }

  return useInitial;
}

/* ------------------------------------------------------------------ */


