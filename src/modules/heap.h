#ifndef __HEAP_H
#define __HEAP_H


/*----------------------------  Heap  -------------------------------*/

typedef struct
          {
          void** data;
          int**  index;
          int  (*compare)(void* a, void* b, void* info);
          int    size;
          int    allocatedsize;
          } HEAP;

/* ----------------------------------------------------------------- */

#define Heap_size(H)        ((H)->size)
#define Heap_elem(H,i)      ((H)->data[i])
#define Heap_isempty(H)     (Heap_size(H) == 0)

void    Heap_init(HEAP* H, int allocatedsize,
                  int (*compf)(void* elem1, void* elem2, void* info));
void    Heap_free(HEAP* H);
void    Heap_update(HEAP* H, int target, void* info);
void    Heap_insert(HEAP* H, void* elem, void* info, int* index);
void*   Heap_remove(HEAP* H, int heapindex, void* info);
#define Heap_removeroot(H,info)   Heap_remove((H), 1, (info))
void    Heap_moveelem(HEAP* H, int target, void* elem, int* index);

/* ----------------------------------------------------------------- */

#endif /* __HEAP_H */


