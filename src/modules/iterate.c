/*-------------------------------------------------------------------*/
/* ITERATE.C    Juha Kivijärvi                                       */
/*                                                                   */
/* Codebook iterator                                                 */
/*                                                                   */
/*-------------------------------------------------------------------*/

#include "cb.h"
#include "iterate.h"

#define ProgName        "ITERATE"
#define VersionNumber   "Version 0.02"
#define LastUpdated     "5.9.96"

/* ----------------------------------------------------------------- */

void InitializeCombinationIterator(int CBsize, int *Array)
{
int i;

for (i=0; i<CBsize; i++)
    Array[i] = i;
    
}

/* ----------------------------------------------------------------- */

int IterateNextCombination(int TSsize, int CBsize, int *Array)
{
int i;

if (Array[0] == TSsize-CBsize) return 0;

i=CBsize-1;

while (Array[i] == TSsize-CBsize+i && i>0)
    i--;
    
Array[i]++;
i++;

while (i<CBsize)
    {
    Array[i]=Array[i-1]+1;
    i++;
    }    

return 1;
}

/* ----------------------------------------------------------------- */

void GenerateIteratedCodebook(TRAININGSET *TS, CODEBOOK *CB, int *Array)
{
int i;

for (i=0; i<BookSize(CB); i++)
    CopyVector(Vector(TS,Array[i]),Vector(CB,i),VectorSize(TS));
    
}


/* ----------------------------------------------------------------- */

void InitializeStirlingIterator(int TSsize, int CBsize, int *Array)
{
int i;

for (i=0; i<TSsize; i++)
    Array[i] = 0;

for (i=0; i<CBsize; i++)
    Array[TSsize-1-i] = CBsize-1-i;

}


/* ----------------------------------------------------------------- */

int IterateNextStirling(int TSsize, int CBsize, int *Array)
   /* With this function one can walk through all the partitionings.
      Only partitionings with no empty partitions are returned.      */

{
int i,j,k;
int used[CBsize];
int maxsofar[TSsize];
int OK;

do
    {
    for (i=0; i<CBsize; i++)
        used[i] = 0;
    maxsofar[0] = Array[0];
    for (i=1; i<TSsize; i++)
        maxsofar[i] = (Array[i] > maxsofar[i-1]) ? Array[i] : maxsofar[i-1];

    i=TSsize-1;

    while ( (Array[i] == CBsize-1 || Array[i] == maxsofar[i-1]+1) && i>0 )
        i--;
    
    if (i==0) return 0;
    Array[i]++;

    for (j=i+1; j<TSsize; j++)
        Array[j]=0;
        
    for (j=0; j<TSsize; j++)
        used[Array[j]]++;
        
    k=TSsize-1;

    for (j=CBsize-1; j>=0; j--)
        if (used[j]==0 && k>i)
            {
            used[Array[k]]--;
            used[j]++;
            Array[k] = j;
            k--;
            }

    OK=1;
    for (i=0; i<CBsize; i++)
        if (used[i]==0) OK=0;
        
    }
while (!(OK));    

return 1;
}
