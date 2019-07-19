/* ----------------------------------------------------------------- */
/* SA.C            Timo Kaukoranta                                   */
/*                                                                   */
/* - Temperature module of Simulated Annealing                       */
/*                                                                   */
/* ----------------------------------------------------------------- */

#define ProgName        "SA"
#define VersionNumber   "Version 0.05"
#define LastUpdated     "20.4.98"

/* In this module enum constants are based on the following definitions.
   Copy these Fact macros into the *.fac file when using this module
   in new application.
*/

/* -------------------------------------------------------------------------- */
/*    Parameter identifier                                                    */
/*      key    no   type value  min  max   default   save to file             */
/*      enum of values                                                        */
/*      name of values                                                        */
/*      printing                                                              */
/* -------------------------------------------------------------------------- */
/*
Fact( VectorRandomizing,
      "Vector Randomizing",
      'R',    1,  ENUM,   0,    0,   3,      0,      NO,
      NoRandomizing, RandomCentroid, RandomTrainingSet, RandomBoth, vr5, vr6, vr7, vr8, vr9, vr10,
      "No Randomizing", "Random Centroid", "Random Training Set", "Random Both", 0, 0, 0, 0, 0, 0,
      YES )

Fact( TemperatureDistribution,
      "Temperature Distribution",
      'R',    2,  ENUM,   0,    0,   1,      0,      NO,
      Even, Gaussian, td3, td4, td5, td6, td7, td8, td9, td10,
      "Even Temperature", "Gaussian Temperature", 0, 0, 0, 0, 0, 0, 0, 0,
      Value(VectorRandomizing) != NoRandomizing )

Fact( TemperatureFunction,
      "Temperature Decrease Function",
      'R',    3,  ENUM,   0,    0,   1,      1,      NO,
      TFLinear, TFlogarithmical, tf3, tf4, tf5, tf6, tf7, tf8, tf9, tf10,
      "Linear Temperature", "Logarithmical Temperature", 0, 0, 0, 0, 0, 0, 0, 0,
      Value(VectorRandomizing) != NoRandomizing )

Fact( TemperatureInitial,
      "Temperature Initial (% of MaxValue)",
      'R',    4,  INT,   0,    0,   100,    20,      NO,
      ti1, ti2, ti3, ti4, ti5, ti6, ti7, ti8, ti9, ti10,
      0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
      Value(VectorRandomizing) != NoRandomizing )

Fact( TemperatureChange,
      "Temperature Change (% of current)",
      'R',    5,  INT,   0,    0, 10000,   100,      NO,
      tc1, tc2, tc3, tc4, tc5, tc6, tc7, tc8, tc9, tc10,
      0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
      Value(VectorRandomizing) != NoRandomizing )
*/

#define TEMPERATUREMAXCHANGE 10000.0


/* ----------------------------------------------------------------- */

#include <stdlib.h>

#include "random.h"
#include "cb.h"
#include "memctrl.h"
#include "sa.h"

/* ----------------------------------------------------------------- */

#define round(a)      ((a) + 0.5)

/*===========================  RANDOM VECTORS ==============================*/


void InitializeSASchedule(int         Distribution,
                          int         DecreaseFunction,
                          double      Change,
                          double      initial,
                          int         ToCodevectors,
                          int         ToTrainingvectors,
                          SASchedule* t)
{
  if( t != NULL )
    {
    t->Distribution           = Distribution;
    t->DecreaseFunction       = DecreaseFunction;
    t->Change                 = Change;
    t->CurrentTemperature     = initial;
    t->ApplyToCodevectors     = ToCodevectors;
    t->ApplyToTrainingvectors = ToTrainingvectors;
    }
}


/* ----------------------------------------------------------------- */


static double RandomNoise(SASchedule* t)
{
  double result = 0.0;

  switch( t->Distribution )
    {
    case 0: /* Even */
      {
      result = (drand() * 2 * SASTemp(t)) - SASTemp(t);
      break;
      }
    case 1: /* Gaussian */
      {
      break;
      }
    default:
      {
      break;
      }
    }
  return( result );
}


/* ----------------------------------------------------------------- */


void DecreaseTemperature(SASchedule* t)
{
  if( t != NULL )
    {
    switch( t->DecreaseFunction )
      {
      case 0: /* TFLinear */
        {
        t->CurrentTemperature -= t->Change;
        break;
        }
      case 1: /* TFlogarithmical */
        {
        t->CurrentTemperature -= (t->Change * t->CurrentTemperature) / TEMPERATUREMAXCHANGE;
        break;
        }
      default:
        break;
      }
    if( t->CurrentTemperature < 0.0 )
      {
      t->CurrentTemperature = 0.0;
      }
    }
}


/* ----------------------------------------------------------------- */


BOOKNODE* RandomizeVectorBySA(SASchedule* t,
                              BOOKNODE*   source,
                              BOOKNODE*   dest,
                              int         Vsize,
                              int         maxvalue)
{
  int i;
  int elem;

  if( t != NULL )
    {

/* WE DO NOT APPLY THIS YET!!!

    if( maxvalue == 1 )
      {
      int  nchanges = (int)round(SASTemp(t) * Vsize);
      int* order    = allocate(Vsize * sizeof(int));

      for( i = 0; i < Vsize; i++ )
        {
        order[i] = i;
        }

      for( i = 0; i < nchanges; i++ )
        {
        elem = order[irand(i, Vsize-1)];
        source->vector[elem] = 1 - source->vector[elem];
        order[elem] = order[i];
        }
      deallocate(order);
      }
    else
*/
         /* Non-binary vector */
      {
      for( i = 0; i < Vsize; i++ )
        {
        elem = (int)round((double)source->vector[i] + RandomNoise(t));
        elem = ( elem < 0 ? 0 :
               ( elem > maxvalue ? maxvalue : elem ) );
        dest->vector[i] = (VECTORELEMENT)elem;
        }
      }
    dest->freq = source->freq;
    }
  else
    {
    CopyNode(source, dest, Vsize);
    }

  return( dest );
}

