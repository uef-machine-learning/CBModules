/*-------------------------------------------------------------------*/
/* CBIMS.C        Timo Kaukoranta                                    */
/*                                                                   */
/*                                                                   */
/* - Iterative Merge and Split                                       */
/*                                                                   */
/*-------------------------------------------------------------------*/

#define ProgName "CBIMS"
#define VersionNumber "V.0.15"
#define LastUpdated "16.6.09" /* PF */

/* ----------------------------------------------------------------- */

#include <string.h>
#include <time.h>
#include <values.h>

#define FACTFILE "cbims.fac"
#include "parametr.c"

#include "file.h"
#include "ims.h"
#include "memctrl.h"
#include "pnn.h"
#include "random.h"
#include "sa.h"
#include "sort.h"
#include "sortcb.h"
#include "split.h"

/*-----------------------------  M i s c  ----------------------------*/

#define MAXERROR (MAXLONG >> 1)

/* =================  T I M E   F U N C T I O N S  ================= */

long SetWatch(long *start) {
  time_t tmp;
  return (*start = (long)time(&tmp));
}

/* ----------------------------------------------------------------- */

void PrintTime(long start) {
  time_t tmp;
  long elapsed = (long)time(&tmp) - start;

  switch (Value(ShowProgress)) {
  case 0:
    break;
  case 1:
    printf("%6li ", elapsed);
    break;
  default:
    printf("Time:%6li\n", elapsed);
    break;
  }
}

/*============================  P R I N T  ===============================*/

void PrintInfo(void) {
  printf("%s %s %s\n", ProgName, VersionNumber, LastUpdated);
  printf("This program improves a codebook by Iterated Merge and Split method.\n");
  printf("Usage: %s [%coption] <training-set> [initial-codebook] <result-codebook>\n", ProgName,
         OPTION_SYMBOL);
  printf("For example: %s set_of_4 genetic gen1\n", ProgName);
  printf("\n  Options:\n");
  PrintOptions();
  printf("%s %s %s\n", ProgName, VersionNumber, LastUpdated);
}

/*-------------------------------------------------------------------*/

static void PrintOperatingInfo(char TSetName[], char InitialCBName[], char ResultCBName[]) {
  if (Value(ShowProgress) >= 2) {
    printf("\n%s %s %s\n", ProgName, VersionNumber, LastUpdated);
    printf("Training Set:     %s\n", TSetName);
    printf("Initial Codebook: %s\n", (InitialCBName[0] == 0x00 ? "TRAINING SET" : InitialCBName));
    printf("Result Codebook:  %s\n", ResultCBName);
    PrintSelectedOptions();
    printf("\n");
  }
}

/*=========================  CODEBOOK HANDLING  ============================*/

static void GenerateRandomCodebook(TRAININGSET *TS, CODEBOOK *CB) {
  int i, j;
  int founddupl;
  int picked;
  int candidate;

  for (i = 0; i < BookSize(CB); i++) {
    do {
      founddupl = NO;
      picked = irand(1, TotalFreq(TS));
      candidate = 0;
      picked -= VectorFreq(TS, candidate);
      while (picked > 0) {
        candidate++;
        picked -= VectorFreq(TS, candidate);
      }
      for (j = 0; j < i; j++) {
        if (EqualVectors(Vector(CB, j), Vector(TS, candidate), VectorSize(CB))) {
          founddupl = YES;
        }
      }
    } while (founddupl);
    CopyNode(&Node(TS, candidate), &Node(CB, i), VectorSize(TS));
  }
}

/* ----------------------------------------------------------------- */

static void GenerateInitialCodebook(char *InitialCBName, CBFILETYPE InitFT, TRAININGSET *TS,
                                    CODEBOOK *CB, PARTITIONING *P) {
  if (Value(ShowProgress) >= 4) {
    printf("Gen.Init.CB: FT=%i\n", InitFT);
  }

  switch (InitFT) {
  case NOTFOUND: {
    /* Create random initial codebook from TS. */
    if (BookSize(TS) < Value(CodebookSize)) {
      printf("ERROR: Training set size (%i) < requested codebook size (%i).\n", BookSize(TS),
             Value(CodebookSize));
      exit(-1);
    }
    CreateNewCodebook(CB, Value(CodebookSize), TS);
    CreateNewPartitioning(P, TS, Value(CodebookSize));
    GenerateRandomCodebook(TS, CB);
    GenerateOptimalPartitioning(TS, CB, P);
    AddGenerationMethod(CB, "IMS(random)");
    break;
  }
  case CBFILE: {
    /* Read initial codebook. */
    ReadCodebook(InitialCBName, CB);
    CreateNewPartitioning(P, TS, BookSize(CB));
    GenerateOptimalPartitioning(TS, CB, P);
    AddGenerationMethod(CB, "IMS");
    break;
  }
  case PAFILE: {
    /* Read initial partitioning and generate initial codebook
       based on it. */
    ReadPartitioning(InitialCBName, P, TS);
    CreateNewCodebook(CB, P->PartitionCount, TS);
    GenerateOptimalCodebook(TS, CB, P);
    AddGenerationMethod(CB, "IMS");
    break;
  }
  default: {
    printf("ERROR: Unknown initial filetype=%i\n", InitFT);
    exit(-1);
  }
  }
}

/* ----------------------------------------------------------------- */

static void SaveResultCodebook(char *ResultCBName, CODEBOOK *CB) {
  int ndups;

  if ((ndups = DuplicatesInCodebook(CB)) > 0) {
    printf("WARNING: Final CB contains %i duplicates.\n", ndups);
  }
  SortCodebook(CB, FREQ_DESCENDING);
  WriteCodebook(ResultCBName, CB, Value(OverWrite));
}

/* ================================================================= */

static void Initialize(CODEBOOK *CB, SASchedule *SAS) {
  int ndups;

  if ((ndups = DuplicatesInCodebook(CB)) > 0) {
    if (Value(ShowProgress) >= 1) {
      printf("WARNING: Init. CB does contain %i duplicates.\n", ndups);
    }
  }

  SetIMSParameters(Value(ShowProgress), Value(CodebookSize), Value(AlternatingOrder),
                   Value(AlternatingAmount), Value(MinCodebookSize), Value(MaxCodebookSize),
                   Value(StepSizeChange), Value(Iterations), Value(IMSGLAIterations),
                   Value(RequireDistortionDecrease));

  SetPNNParameters(Value(ShowProgress) - 2, 0,  /* Value(MergeErrorType)==Equitz */
                   0,                           /* Value(CentroidCalculation)==CodevectorBased */
                   0,                           /* Value(MergeMethod)==Kissing */
                   NO,                          /* Value(PartitionRemapping) because of Equitz */
                   Value(PNNGLAIterations), 50, /* Value(MergePercent) */
                   8,                           /* Value(MaxBucket) */
                   4);                          /* Value(MinBucket) */

  SetSplitParameters(Value(ShowProgress) - 2, Value(PartitionErrorType),
                     7,   /* Value(NewVectors) == MeansOfCurrentAndFurthest */
                     YES, /* Value(NewVectorsHeuristically) == YES */
                     YES, /* Value(Hyperplane) */
                     Value(HyperplanePivot), Value(HyperplaneLine), Value(PartitionRemapping),
                     NO, /* Value(LocalGLA), */
                     0,  /* Value(LocalGLAIterations), */
                     Value(SplitGLAIterations));

  InitializeSASchedule(
      Value(TemperatureDistribution), Value(TemperatureFunction), Value(TemperatureChange),
      (float)Value(TemperatureInitial) * (float)(CB->MaxValue) / 100.0,
      (Value(VectorRandomizing) == RandomCentroid || Value(VectorRandomizing) == RandomBoth),
      (Value(VectorRandomizing) == RandomTrainingSet || Value(VectorRandomizing) == RandomBoth),
      SAS);
}

/* ================================================================= */

static void ShutUp(TRAININGSET *TS, CODEBOOK *CB, PARTITIONING *P) {
  FreeCodebook(TS);
  FreeCodebook(CB);
  FreePartitioning(P);
}

/*===========================  M A I N  ==============================*/

int main(int argc, char *argv[]) {
  CODEBOOK CB;
  TRAININGSET TS;
  SASchedule SAS;
  PARTITIONING P;
  char TSetName[MAXFILENAME] = {0x00};
  char InitialCBName[MAXFILENAME] = {0x00};
  char ResultCBName[MAXFILENAME] = {0x00};
  CBFILETYPE InitFT;
  double FinalError;
  int iteration;
  long watch;

  ParameterInfo paraminfo[3] = {{TSetName, FormatNameTS, 0, INFILE},
                                {InitialCBName, FormatNameCB, 1, INFILE},
                                {ResultCBName, FormatNameCB, 0, OUTFILE}};

  ParseParameters(argc, argv, 3, paraminfo);
  InitFT = DetermineCBFileTypeConsideringOrder(InitialCBName, 23);
  initrandom(Value(RandomSeed));
  SetWatch(&watch);

  ReadTrainingSet(TSetName, &TS);
  GenerateInitialCodebook(InitialCBName, InitFT, &TS, &CB, &P);

  Initialize(&CB, &SAS);
  PrintOperatingInfo(TSetName, InitialCBName, ResultCBName);

  IMS(&TS, &CB, &SAS, &P, &iteration);

  PrintTime(watch);

  if (TS.InputFormat == TXT) {
    // If input data is in txt format and scaled to integers,
    // no sense to save in binary format
    SaveCB2TXT(&CB, ResultCBName, TS.MinMax, NULL, NULL);
  } else {
    WriteCodebook(ResultCBName, &CB, Value(OverWrite));
  }

  FinalError = AverageErrorCBFast(&TS, &CB, &P, MSE);
  if (Value(ShowProgress) >= 2) {
    printf("Final dist. = %9.4f\n", PrintableError(FinalError, &CB));
    printf("Iterations  = %3i\n", iteration);
  } else {
    printf("%9.4f ", PrintableError(FinalError, &CB));
    printf("%3i ", iteration);
  }

  ShutUp(&TS, &CB, &P);
  checkmemory();

  return (0);
}
