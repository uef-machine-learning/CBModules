/*-------------------------------------------------------------------*/
/* CBKM.C         Marko Tuononen                                     */
/*                                                                   */
/* Console interface for repeated K-means algorithm                  */
/*                                                                   */
/*                                                                   */
/*-------------------------------------------------------------------*/

#define ProgName       "CBKM"
#define VersionNumber  "Version 0.65" /* SS */
#define LastUpdated    "4.4.2017"
#define FACTFILE       "cbkm.fac"

/* ----------------------------------------------------------------- */

#include "parametr.c"
#include "cb.h"
#include "file.h"
#include "interfc.h"
#include "memctrl.h"
#include "random.h"
#include "reporting.h"
#include "kmeans.h"
#include "textfile.h"


/* ======================== PRINT ROUTINES =========================== */


void PrintInfo(void)
{
  PrintMessage("%s\t%s\t%s\n"
        "Repeated K-means algorithm.\n"
        "Usage: %s [%coption] <training set> [initial cb/pa] "
        "<result codebook>\n"
        "For example: %s bridge initial tmp\n\n  Options:\n",
        ProgName, VersionNumber, LastUpdated, ProgName,
        OPTION_SYMBOL, ProgName);
  PrintOptions();
  PrintMessage("\n");
}


/* ------------------------------------------------------------------ */

static char* PrintInitialData(char* TSName, char* InName,
                              char* OutCBName, char* OutPAName,
                             int useInitial)
{
  char* str;
  str = KMeansInfo();

  if( Value(QuietLevel) >= 2 )
    {
    printf("\n%s %s %s\n", ProgName, VersionNumber, LastUpdated);
    PrintMessage("Training Set:     %s\n", TSName);
    PrintMessage("Codebook:         %s\n", OutCBName);
    if (useInitial)
      {
      PrintMessage("Initial %s",
                  ((useInitial == 1) ? "codebook:         "
                                     : "partition:        "));
      PrintMessage("%s\n", InName);
      }
    if (*OutPAName)
      {
      PrintMessage("Result partition: %s\n", OutPAName);
      }
    printf("\n");
    PrintSelectedOptions();
    printf("\n");
    }
  return str;
}


/* ===========================  MAIN  ================================ */


int main(int argc, char* argv[])
{
  TRAININGSET   TS;
  CODEBOOK      CB;
  PARTITIONING  P;
  char          TSName[MAXFILENAME] = {'\0'};
  char          InName[MAXFILENAME] = {'\0'};
  char          OutCBName[MAXFILENAME] = {'\0'};
  char          OutPAName[MAXFILENAME] = {'\0'};
  int           useInitial = 0;
  char*         genMethod;
  ParameterInfo paraminfo[3] = { { TSName,  FormatNameTS, 0, INFILE },
                                 { InName,  FormatNameCB, 1, INFILE },
                                 { OutCBName,  FormatNameCB, 0, OUTFILE } };


  ParseParameters(argc, argv, 3, paraminfo);
  initrandom(Value(RandomSeed));


  if (Value(SavePartition))
  {
      PickFileName(OutCBName, OutPAName);
      CheckFileName(OutPAName, FormatNamePA);
  }

  TS = CheckParameters(TSName, OutCBName, OutPAName, InName,
       Value(Clusters), Value(OverWrite));


  //float* density = dimBasedDensity(&TS, 10);
  //return;

  //kNNDensityForDim(&TS, 100, 0);
  //return;

  useInitial = ReadInitialCBorPA(InName, Value(Clusters), &TS, &CB, &P);

  int initMeth = Value(InitMethod);
  int initMethVar = Value(InitMethodVar);

  if(Value(HybridMethod) == KMPeaks ) {

      printf("KMPEAKS\n");
      //Value(InitMethod) = 101;
      //ModuleParameters.InitMethod=101;
      initMeth=101;
      initMethVar = Value(HybridMethodVar);
      //KMeansDensityPeaks(&TS, &CB, &P, Value(Clusters), Value(Repeats), Value(QuietLevel));
      //return 0 ;
  }
  kmOpt.sample=Value(Sample)/10000.0f;

  kmOpt.densityPeaksFilterMethod=Value(InitMethodOpt);
  kmOpt.densityMethod=Value(InitMethodThirdOpt);
  kmOpt.QuietLevel=Value(QuietLevel);

  genMethod = PrintInitialData(TSName, InName, OutCBName,
              OutPAName, useInitial);



  if (PerformKMeans(&TS, &CB, &P, Value(Clusters), Value(Repeats),
              initMeth,initMethVar, Value(QuietLevel), useInitial,Value(Iterations)))
    {
    ErrorMessage("\nERROR: Clustering failed!\n\n");
    FreeCodebook(&TS);
    FreeCodebook(&CB);
    FreePartitioning(&P);
    free(genMethod);
    ExitProcessing(FATAL_ERROR);
    }

  AddGenerationMethod(&CB, genMethod);

  if(TS.InputFormat==TXT) {
      // If input data is in txt format and scaled to integers,
      // no sense to save in binary format
      SaveCB2TXT(&CB, OutCBName, TS.MinMax, NULL, NULL);
  }
  else {
      WriteCodebook(OutCBName, &CB, Value(OverWrite));
  }



  if (Value(SavePartition))
    {
    WritePartitioning(OutPAName, &P, &TS, Value(OverWrite));
    }

  FreeCodebook(&TS);
  FreeCodebook(&CB);
  FreePartitioning(&P);
  free(genMethod);

  return EVERYTHING_OK;
}


/* ----------------------------------------------------------------- */
