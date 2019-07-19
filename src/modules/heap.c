/*-------------------------------------------------------------------*/
/* HEAP.C          Timo Kaukoranta                                   */
/*                                                                   */
/* General data structure of heap                                    */
/*                                                                   */
/*-------------------------------------------------------------------*/

#define ProgName        "HEAP"
#define VersionNumber   "Version 0.02"
#define LastUpdated     "27.3.98"

/* ----------------------------------------------------------------- */


#include <assert.h>
#include <string.h>
#include <stdio.h>

#include "heap.h"
#include "memctrl.h"
#include "owntypes.h"


/* ----------------------------------------------------------------- */

void Heap_init(HEAP* H,
               int   allocatedsize,
               int (*compf)(void* elem1, void* elem2, void* info))
{
  assert( 0 <= allocatedsize );

  H->allocatedsize = allocatedsize+1;
  H->data    = (void**)allocate(H->allocatedsize * sizeof(void*));
  H->index   = (int**)allocate(H->allocatedsize * sizeof(int*));
  H->compare = compf;
  H->size    = 0;
}

/* ----------------------------------------------------------------- */

void Heap_free(HEAP* H)
{
  deallocate(H->data);
  deallocate(H->index);
}

/* ----------------------------------------------------------------- */

static void Heap_siftdown(HEAP* H, int target, void* info)
{
  int   i;
  YESNO flag = YES;
  int*  tmpi;
  void* tmpd;

  while( flag && target > 1 )
    {
    i = target / 2;
    if( H->compare(H->data[target], H->data[i], info) >= 0 )
      {
      flag = NO;
      }
    else
      {
      tmpd = H->data[i];
      H->data[i] = H->data[target];
      H->data[target] = tmpd;

      tmpi = H->index[i];
      H->index[i] = H->index[target];
      H->index[target] = tmpi;

      *(H->index[i]) = i;
      *(H->index[target]) = target;

      target = i;
      }
    }
}

/* ----------------------------------------------------------------- */

static void Heap_siftup(HEAP* H, int target, void* info)
{
  int   j;
  void* tmpd;
  int*  tmpi;

  while( 2*target <= Heap_size(H) )
    {
    j = 2*target;
    if( j < Heap_size(H) )
      {
      /* Select smaller of two sons. */
      if( H->compare(H->data[j], H->data[j+1], info) > 0 )
        {
        j++;
        }
      }
    if( H->compare(H->data[target], H->data[j], info) > 0 )
      {
      tmpd = H->data[j];
      H->data[j] = H->data[target];
      H->data[target] = tmpd;

      tmpi = H->index[j];
      H->index[j] = H->index[target];
      H->index[target] = tmpi;

      *(H->index[j]) = j;
      *(H->index[target]) = target;

      target = j;
      }
    else
      {
      target = Heap_size(H) + 1;
      }
    }
}

/* ----------------------------------------------------------------- */

void Heap_update(HEAP* H, int target, void* info)
{

  if( target > 1 &&
      H->compare(H->data[target/2], H->data[target], info) > 0 )
    {
    Heap_siftdown(H, target, info);
    }
  else
    {
    Heap_siftup(H, target, info);
    }
}

/* ----------------------------------------------------------------- */

void Heap_insert(HEAP* H, void* elem, void* info, int* index)
{
  int**  tmpi;
  void** tmpd;

  if( H->allocatedsize == Heap_size(H) )
    {
    tmpd = H->data;
    tmpi = H->index;
    H->allocatedsize = H->allocatedsize * 2 + 1;

    H->data  = (void**)allocate(H->allocatedsize * sizeof(void*));
    H->index = (int**)allocate(H->allocatedsize * sizeof(int*));
    int ret1,ret2;
    ret1 = memcpy(H->data, tmpd, (Heap_size(H)+1) * sizeof(void*));
    ret2 = memcpy(H->index, tmpi, (Heap_size(H)+1) * sizeof(int*));

    if(!ret1 || !ret2) {
        printf("memcpy fail\n");
    }

    if(!tmpd || !tmpi) {
        printf("memcpy fail\n");
    }



    deallocate(tmpd);
    deallocate(tmpi);
    }

  Heap_size(H)++;         /* A new element */
  H->data[Heap_size(H)]  = elem;
  H->index[Heap_size(H)] = index;
  *(H->index[Heap_size(H)]) = Heap_size(H);

  Heap_siftdown(H, Heap_size(H), info);
}

/* ----------------------------------------------------------------- */

void* Heap_remove(HEAP* H, int heapindex, void* info)
{
  void* elem;

  assert( ! Heap_isempty(H) );
  assert( 1 <= heapindex );
  assert( heapindex <= Heap_size(H) );

  elem = H->data[heapindex];
  *(H->index[heapindex]) = 0; /* Outside of 1..size */

  H->data[heapindex] = H->data[Heap_size(H)];
  H->index[heapindex] = H->index[Heap_size(H)];
  *(H->index[heapindex]) = heapindex;

  Heap_size(H)--;

  Heap_siftup(H, heapindex, info);

  return( elem );
}

/* ----------------------------------------------------------------- */

void Heap_moveelem(HEAP* H, int target, void* elem, int* index)
{
  H->data[target]  = elem;
  H->index[target] = index;
  *(H->index[target]) = target;
}
