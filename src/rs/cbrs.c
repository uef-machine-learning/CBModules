/*-------------------------------------------------------------------*/
/* CBRS.C         Pasi Franti                                        */
/*                                                                   */
/* Random Swap (RS) algorithm                                        */
/*                                                                   */
/* ChangeLog:                                                        */
/*                                                                   */
/* 0.64: 25.6.16  PF: Monitor progress (CI-value)                    */
/* 0.63: 23.6.16  PF: Added Random Seed + Text modifications         */
/* 0.62: 28.2.10  AH: Traveller search                               */
/*                                                                   */
/*-------------------------------------------------------------------*/

#define ProgName        "CBRS"
#define VersionNumber   "Version 0.64"  /* PF */
#define LastUpdated     "25.6.2016"
#define FACTFILE        "cbrs.fac"

/* ------------------------------------------------------------------- */

#include "parametr.c"
#include "cb.h"
#include "file.h"
#include "interfc.h"
#include "memctrl.h"
#include "random.h"
#include "reporting.h"
#include "rs.h"
#include "textfile.h"


/* ======================== PRINT ROUTINES =========================== */


void PrintInfo(void)
{
  PrintMessage("%s\t%s\t%s\n\n"
        "Random Swap algorithm.\n"
        "Use: %s [%coption] <dataset> [initial cb/pa] <codebook>\n"
        "For example: %s bridge initial tmp\n\n  Options:\n",
        ProgName, VersionNumber, LastUpdated, ProgName, OPTION_SYMBOL,
        ProgName);
  PrintOptions();
  PrintMessage("\n");
}


/* ------------------------------------------------------------------ */


static char* PrintInitialData(char *TSName, char *InName,
             char *OutCBName, char* OutPAName, int useInitial)
{
  char* str;
  str = RSInfo();

  if (Value(QuietLevel) >= 2)
    {
    PrintMessage("\n%s\n\n", str);
    PrintMessage("Dataset                   = %s \n", TSName);

    if (useInitial)
      {
      PrintMessage("Initial ");
      if (useInitial == 1)  PrintMessage("codebook");
      else                  PrintMessage("partitioning");
      PrintMessage("          = %s \n", InName);
      if(Value(MonitorProgress))
        {
        PrintMessage("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! \n");
        PrintMessage("!!! Progress Monitor mode selected:       !!! \n");
        PrintMessage("!!! Initial codebook is used as REFERENCE !!! \n");
        PrintMessage("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! \n");
        }
      }

    PrintMessage("Codebook                  = %s\n", OutCBName);

    if (*OutPAName)
      {
      PrintMessage("Partition file            = %s\n", OutPAName);
      }

    if (Value(Iterations)==0)
       {
       PrintMessage("Number of iterations      = AUTOMATIC\n");
       }

    PrintSelectedOptions();
    PrintMessage("\n");
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

    useInitial = ReadInitialCBorPA(InName, Value(Clusters), &TS, &CB, &P);

    genMethod = PrintInitialData(TSName, InName, OutCBName,
            OutPAName, useInitial);

    if (PerformRS(&TS, &CB, &P, Value(Iterations),
                Value(KMeansIterations), Value(Deterministic), Value(TravellerSearch),
                Value(QuietLevel), useInitial, Value(MonitorProgress)))
    {
        ErrorMessage("ERROR: Clustering failed!\n");
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
