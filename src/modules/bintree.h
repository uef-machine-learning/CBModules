#if ! defined(__BINTREE_H)
#define __BINTREE_H

/*--------------------------  Node Stack  -----------------------------*/


struct  STACKSTRUCT { void*               pointer;
                      struct STACKSTRUCT* next;
                    };

typedef struct STACKSTRUCT  STACKNODE;
typedef STACKNODE*          STACK;


/*--------------------------  Binary Tree  ----------------------------*/


struct  NODESTRUCT  { void*              data;
                      struct NODESTRUCT* left;
                      struct NODESTRUCT* right;
                    };

typedef struct NODESTRUCT   NODETYPE;
typedef NODETYPE*           NODE;

struct  ROOTSTRUCT  { NODE   first;
                      int  (*compare)(void* a, void* b, void* info);
                      int    nnodes;
                    };

typedef struct ROOTSTRUCT   BINTREE;


/* ----------------------------------------------------------------- */

#define BintreeSize(tree) ((tree)->nnodes)


void    InitBintree(BINTREE* tree,
                    int (*compf)(void* a, void* b, void* info));
void*   DeleteMinimumFromBintree(BINTREE* tree);
void*   DeleteMaximumFromBintree(BINTREE* tree);
void*   DeleteNodeFromBintree(BINTREE* tree, void* d, void* info);
void    InsertToBintree(BINTREE* tree, void* d, void* info);
void*   InsertToBintreeNoDuplicates(BINTREE* tree, void* d, void* info);
void*   FindFromBintree(BINTREE* tree, void* d, void* info);
void    FreeBintree(BINTREE* tree);
int     CheckBintree(BINTREE* tree);

void    ClearBintreeIterator(STACK* s);
void    InitPreOrderBintree(BINTREE* tree, STACK* s);
void*   PreOrderBintree(STACK* s);

void    InitInOrderBintree(BINTREE* tree, STACK* s);
void*   InOrderBintree(STACK* s);

typedef enum {
	PREORDER = 1,
	INORDER,
	POSTORDER
} BINTREE_ORDER;

void    IterateBintreeWithCallback(BINTREE* tree,
                                   BINTREE_ORDER order,
                                   int (*callbackf)(void*));

/* ----------------------------------------------------------------- */

#endif /* __BINTREE_H */
