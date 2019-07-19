/*-------------------------------------------------------------------*/
/* BINTREE.C      Timo Kaukoranta                                    */
/*                                                                   */
/*                                                                   */
/* - Implementation of stack and binary tree                         */
/*                                                                   */
/*-------------------------------------------------------------------*/

#define ProgName        "BINTREE"
#define VersionNumber   "Version 0.09a"
#define LastUpdated     "28.7.99"

/* ----------------------------------------------------------------- */

#include <assert.h>
#include <stdio.h>

#include "memctrl.h"
#include "bintree.h"
#include "interfc.h"


/*====================  N O D E   S T A C K  =========================*/


#define StackInit(s)     ( (*(s)) = NULL )
#define StackEmpty(s)    ( (*(s)) == NULL )


/* ----------------------------------------------------------------- */


void Push(void* p, STACK* s)
{
  STACK tmp;

  tmp = (STACK)allocate(sizeof(STACKNODE));
  if( tmp == NULL )
    {
    ErrorMessage("\nERROR: allocate Push\n");
    ExitProcessing( -1 );
    }
  tmp->pointer = p;
  tmp->next    = *s;
  *s           = tmp;
}


/*-------------------------------------------------------------------*/


void* Pop(STACK* s)
{
  void* p;
  STACK tmp;

  if( StackEmpty(s) )
    {
    ErrorMessage("\nERROR: Stack empty, cannot POP.\n");
    ExitProcessing(-1);
    }
  tmp   = (*s);
  p     = (*s)->pointer;
  *s    = (*s)->next;
  deallocate(tmp);

  return( p );
}


/* ----------------------------------------------------------------- */


void ClearStack(STACK* s)
{
  while( ! StackEmpty(s) )
    {
    Pop(s);
    }
}


/*===================  B I N A R Y   T R E E   ========================*/


void InitBintree(BINTREE* tree,
                 int     (*compf)(void* a, void* b, void* info))
{
  tree->first   = NULL;
  tree->compare = compf;
  tree->nnodes  = 0;
}


/* ----------------------------------------------------------------- */


static NODE CreateNode(void* d)
{
  NODE t;

  /* Allocate memory for new node */
  t = (NODE)allocate( sizeof(NODETYPE) );
  if( t == NULL )
    {
    ErrorMessage("\nERROR: allocate CreateNode.\n");
    ExitProcessing(-1);
    }

  /* Initialize */
  t->data  = d;
  t->left  = NULL;
  t->right = NULL;

  return( t );
}


/* ----------------------------------------------------------------- */


static void SearchMinAndFather(NODE* father, NODE* minimum)
{
  if( *minimum != NULL )
    {
    while( (*minimum)->left != NULL )
      {
      *father = *minimum;
      *minimum = (*minimum)->left;
      }
    }
}


/* ----------------------------------------------------------------- */


static void DeleteMinimumBintreeNode(NODE father, NODE target)
{
  NODE rminfather;
  NODE rmintarget;

  assert( target != NULL );
  assert( target->left == NULL );
  assert( father != NULL );
  assert( father->left == target || father->right == target );

  rminfather = target;
  rmintarget = target->right;
  SearchMinAndFather(&rminfather, &rmintarget);

  if( rmintarget == NULL ) 
    {
    /* Target has not a right subtree.
       Put (possible) targets' left child to
       fathers' left or right child.
    */
    if( father->left == target )
      {
      /* Target is the left child of its father. */
      father->left = target->left;
      }
    else
      {
      /* Target is the right child of its father. */
      father->right = target->left;
      }
    deallocate(target);
    }
  else
    {
    target->data = rmintarget->data;

    /* Target has a right subtree.
    */
    if( rminfather->left == rmintarget )
      {
      rminfather->left = rmintarget->right;
      }
    else
      {
      rminfather->right = rmintarget->right;
      }
    deallocate(rmintarget);
    }
}



/*-------------------------------------------------------------------*/


void* DeleteMinimumFromBintree(BINTREE* tree)
{
  void* d;
  NODE  father = NULL;
  NODE  target = tree->first;

  if( tree->first != NULL )
    {
    SearchMinAndFather(&father, &target);
    d = target->data;
    if( father == NULL ) /* Root is minimum. */
      {
      tree->first = target->right;
      deallocate(target);
      }
    else
      {
      DeleteMinimumBintreeNode(father, target);
      }
    tree->nnodes--;
/*
ErrorMessage("BINTREE: DELMIN nodes=%i first=%p\n", tree->nnodes, tree->first);
*/
    return( d );
    }

  return( NULL );
}


/* ----------------------------------------------------------------- */


static void SearchMaxAndFather(NODE* father, NODE* maximum)
{
  if( *maximum != NULL )
    {
    while( (*maximum)->right != NULL )
      {
      *father = *maximum;
      *maximum = (*maximum)->right;
      }
    }
}


/* ----------------------------------------------------------------- */


static void DeleteMaximumBintreeNode(NODE father, NODE target)
{
  NODE lmaxfather;
  NODE lmaxtarget;

  assert( target != NULL );
  assert( target->right == NULL );
  assert( father != NULL );
  assert( father->left == target || father->right == target );

  lmaxfather = target;
  lmaxtarget = target->left;
  SearchMaxAndFather(&lmaxfather, &lmaxtarget);

  if( lmaxtarget == NULL )
    {
    /* Target has not a left subtree.
       Put targets' (possible) right child to
       fathers' right or left child.
    */
    if( father->right == target )
      {
      /* Target is the left child of its father. */
      father->right = target->right;
      }
    else
      {
      /* Target is the left child of its father. */
      father->left = target->right;
      }
    deallocate(target);
    }
  else
    {
    target->data = lmaxtarget->data;

    /* Target has a left subtree.
    */
    if( lmaxfather->right == lmaxtarget )
      {
      lmaxfather->right = lmaxtarget->left;
      }
    else
      {
      lmaxfather->left = lmaxtarget->left;
      }
    deallocate(lmaxtarget);
    }
}



/*-------------------------------------------------------------------*/


void* DeleteMaximumFromBintree(BINTREE* tree)
{
  void* d;
  NODE  father = NULL;
  NODE  target = tree->first;

  if( tree->first != NULL )
    {
    SearchMaxAndFather(&father, &target);
    d = target->data;
    if( father == NULL ) /* Root is maximum. */
      {
      tree->first = target->left;
      deallocate(target);
      }
    else
      {
      DeleteMaximumBintreeNode(father, target);
      }
    tree->nnodes--;
/*
ErrorMessage("BINTREE: DELMAX nodes=%i first=%p\n", tree->nnodes, tree->first);
*/
    return( d );
    }

  return( NULL );
}


/* ================================================================= */


static void SearchNodeAndFather(BINTREE* tree,
                                void*    d,
                                void*    info,
                                NODE*    father,
                                NODE*    t)
{
  while( *t != NULL )
    {
    switch( tree->compare(d, (*t)->data, info) )
      {
      case -1: *father = *t;
               *t = (*t)->left;
               break;
      case  0: /* Found. */
               return;
      case  1: *father = *t;
               (*t) = (*t)->right;
               break;
      default: ErrorMessage("ERROR: SearchNodeAndFather undefined compare value.\n");
               exit( -1 );
      }
    }
ErrorMessage("ERROR: NOT FOUND\n");
  /* Not found. */
  *father = NULL;
  return;
}


/* ----------------------------------------------------------------- */


void* DeleteNodeFromBintree(BINTREE* tree, void* d, void* info)
{
  NODE  father = NULL;
  NODE  target = tree->first;

  if( tree->first != NULL )
    {
    SearchNodeAndFather(tree, d, info, &father, &target);
    if( target != NULL ) /* 'd' is in the tree. */
      {
      if( target->left == NULL )
        {
        if( father == NULL ) /* Root. */
          {
          tree->first = target->right;
          }
        else
          {
          if( father->left == target )
            {
            father->left = target->right;
            }
          else
            {
            father->right = target->right;
            }
          }
        deallocate(target);
        }
      else /* target->left != NULL */
        {
        NODE lfather = target;
        NODE ltarget = target->left;
        SearchMaxAndFather(&lfather, &ltarget);
        if( lfather == target )
          {
          ltarget->right = target->right;
          if( father == NULL )
            {
            tree->first = target->left;
            }
          else
            {
            if( father->left == target )
              {
              father->left = target->left;
              }
            else
              {
              father->right = target->left;
              }
            }
          deallocate(target);
          }
        else
          {
          target->data = ltarget->data;
          DeleteMaximumBintreeNode(lfather, ltarget);
          }
        }
      tree->nnodes--;

      return( d );
      }
    }

  return( NULL );
}



/*===================================================================*/


void InsertToBintree(BINTREE* tree, void* d, void* info)
/*  This insert allows duplicates in the bintree.
 */
{
  NODE t = tree->first;
  int  inserted;

  if( t == NULL )
    {
    tree->first = CreateNode(d);
/* ErrorMessage("ROOT\n"); */
    }
  else
    {
    inserted = 0;
    while( !inserted )
      {
      switch( tree->compare(d, t->data, info) )
        {
        case -1: if( t->left == NULL )
                   {
                   t->left = CreateNode(d);
                   inserted = 1;
/* ErrorMessage("LEFT INSERT\n"); */
                   }
                 else
                   {
                   t = t->left;
/* ErrorMessage("LEFT "); */
                   }
                 break;
        case  0: /* Equals to the right child. */
        case  1: if( t->right == NULL )
                   {
                   t->right = CreateNode(d);
                   inserted = 1;
/* ErrorMessage("RIGHT INSERT\n"); */
                   }
                 else
                   {
                   t = t->right;
/* ErrorMessage("RIGHT "); */
                   }
                 break;
        default: ErrorMessage("ERROR: InsertToBintree undefined compare value.\n");
                 exit( -1 );
        }
      }
    }
  tree->nnodes++;
}


/*-------------------------------------------------------------------*/


void* InsertToBintreeNoDuplicates(BINTREE* tree, void* d, void* info)
/*  Returns: NULL   : 'd' is inserted into the tree.
             ! NULL : 'd' is already in the tree, returns pointer to
                      an existing node.
 */
{
  NODE t = tree->first;
  int  finish = 0;
  int  inserted = 0;

  if( t == NULL )
    {
    tree->first = CreateNode(d);
    inserted = 1;
    }
  else
    {
    while( !finish && !inserted )
      {
      switch( tree->compare(d, t->data, info) )
        {
        case  0: finish = 1;
                 break;
        case -1: if( t->left == NULL )
                   {
                   t->left = CreateNode(d);
                   inserted = 1;
                   }
                 else
                   {
                   t = t->left;
                   }
                 break;
        case  1: if( t->right == NULL )
                   {
                   t->right = CreateNode(d);
                   inserted = 1;
                   }
                 else
                   {
                   t = t->right;
                   }
                 break;
        default: ErrorMessage("ERROR: InsertToBintree undefined compare value.\n");
                 exit( -1 );
        }
      }
    }

  if( inserted )
    {
    tree->nnodes++;
    return( NULL );
    }
  else
    {
    return( t->data );
    }
}


/* ----------------------------------------------------------------- */


void* FindFromBintree(BINTREE* tree, void* d, void* info) {
  NODE  father = NULL;
  NODE  target = tree->first;
  if( tree->first != NULL ) {
    SearchNodeAndFather(tree, d, info, &father, &target);
    if (target) return target->data;
  }
  return NULL;
}


/* ----------------------------------------------------------------- */


void FreeBintree(BINTREE* tree)
{
  NODE  t = tree->first;
  STACK s;

  StackInit(&s);

  /* Push root of the tree to stack */
  if( t != NULL )
    {
    Push(t, &s);

    /* Traverse the tree */
    while( ! StackEmpty(&s) )
      {
      t = (NODE)Pop(&s);
      if( t->left != NULL )
        {
        Push(t->left, &s);
        }
      if( t->right != NULL )
        {
        Push(t->right, &s);
        }
      deallocate(t);
      }
    }

  /* Clear the root */
  tree->first   = NULL;
/*   tree->compare = NULL; */
  tree->nnodes  = 0;
}


/* ==========================  BINTREE SECURITY  =================== */


int CheckBintree(BINTREE* tree)
/* Returns: == 0 OK.
 *          != 0 Difference between reported and counted nodes.
 */
{
  NODE  t = tree->first;
  STACK s;
  int   count = 0;

  StackInit(&s);

  /* Push root of the tree to stack */
  if( t != NULL )
    {
    Push(t, &s);
ErrorMessage("S ");

    /* Traverse the tree */
    while( ! StackEmpty(&s) )
      {
      t = (NODE)Pop(&s);
      if( t->left != NULL )
        {
        Push(t->left, &s);
ErrorMessage("L ");
        }
      if( t->right != NULL )
        {
        Push(t->right, &s);
ErrorMessage("R ");
        }

      count++;
      if( count > tree->nnodes + 20 )
        {
        break;
        }
      }
    }

  return( count - tree->nnodes );
}


/* =========================  BINTREE ITERATORS  =================== */


void ClearBintreeIterator(STACK* s)
{
  ClearStack(s);
}


/* ----------------------------------------------------------------- */


void InitPreOrderBintree(BINTREE* tree, STACK* s)
{
  StackInit(s);
  if( tree->first != NULL )
    {
    Push(tree->first, s);
    }
}


/* ----------------------------------------------------------------- */


void* PreOrderBintree(STACK* s)
{
  NODE t;

  if( StackEmpty(s) )
    {
    return( NULL );
    }
  else
    {
    t = (NODE)Pop(s);
    if( t->left )
      {
      Push(t->left, s);
      }
    if( t->right )
      {
      Push(t->right, s);
      }
    return( t->data );
    }
}


/* ----------------------------------------------------------------- */


void InitInOrderBintree(BINTREE* tree, STACK* s)
{
  NODE t= tree->first;

  StackInit(s);
  while( t != NULL )
    {
    Push(t, s);
    t = t->left;
    }
}


/* ----------------------------------------------------------------- */


void* InOrderBintree(STACK* s)
{
  NODE t;
  NODE t2;

  if( StackEmpty(s) )
    {
    return( NULL );
    }
  else
    {
    t = (NODE)Pop(s);
    t2 = t->right;
    while( t2 != NULL )
      {
      Push(t2, s);
      t2 = t2->left;
      }
    return( t->data );
    }
}


/* ----------------------------------------------------------------- */

static int IterateWithCallback(NODE t, BINTREE_ORDER order,
	int (*callbackf)(void*))
{
	int stop = 0;
	if (order == PREORDER)
		stop = (*callbackf)(t->data);
	if (t->left && !stop)
		stop = IterateWithCallback(t->left, order, callbackf);
	if ((order == INORDER) && !stop)
		stop = (*callbackf)(t->data);
	if (t->right && !stop)
		stop = IterateWithCallback(t->right, order, callbackf);
	if ((order == POSTORDER) && !stop)
		stop = (*callbackf)(t->data);
	return stop;
}


void IterateBintreeWithCallback(BINTREE* tree, BINTREE_ORDER order,
	int (*callbackf)(void*))
{
	NODE t = tree->first;
    if (t) IterateWithCallback(t, order, callbackf);
}

