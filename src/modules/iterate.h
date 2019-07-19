#if ! defined(__ITERATE_H)
#define __ITERATE_H


/* ----------------------------------------------------------------- */

void InitializeCombinationIterator(int CBsize, int *Array);
int IterateNextCombination(int TSsize, int CBsize, int *Array);
void GenerateIteratedCodebook(TRAININGSET *TS, CODEBOOK *CB, int *Array);
void InitializeStirlingIterator(int TSsize, int CBsize, int *Array);
int IterateNextStirling(int TSsize, int CBsize, int *Array);

/* ----------------------------------------------------------------- */


#endif /* __ITERATE_H */
