
#ifndef KNNGRAPH_H_
#define KNNGRAPH_H_

#include <string.h> // For memmove
#include <stdio.h>
#include <stdlib.h>

#define FLOAT float
/*#define UINT unsigned int;*/
#define UINT int
#define bool int
#define MAX_FLOAT 1.0e30; //TODO
#define false 0
#define true 1

//TODO
#define safemalloc malloc

typedef struct {
    UINT id;
    FLOAT dist;
    bool new_item;
} kNNItem;

typedef struct  {
    // List of <size> number of nearest neighbors
    kNNItem * items;
    FLOAT max_dist;
    unsigned int size;
    // <id> of point which nearest neighbors this represents
    UINT id;
    bool is_exact;
} kNNList;

// format

#define AUTODETECT 0
#define SAMPLED_BRUTEFORCE 2
#define RANDOM_SAMPLED_BRUTEFORCE 3
#define KGRAPH 4

typedef struct  {
    int size;
    int k;
    int format;
    kNNList * list;
} kNNGraph;

kNNGraph * init_kNNGraph(int N, int K, int maxK);
void free_kNNGraph(kNNGraph * kNN);

/*kNNGraph * init_kNNGraph(int N, int K ) {*/
/*return init_kNNGraph(N, K, K);*/
/*}*/

void debug_graph(kNNGraph* knng);
int updatekNN(kNNGraph* kNN, UINT p1, UINT p2, FLOAT dist);
/*inline */
UINT get_kNN_item_id(kNNGraph * kNN, int i_list, int i_k);

inline void set_kNN_val(kNNGraph * kNN, int i_list, int i_k, UINT id, FLOAT dist, bool new_item);
inline kNNList* get_kNN_list(kNNGraph * kNN, int i_list);
inline void set_kNN_id(kNNGraph * kNN, UINT i_list, UINT id);
inline kNNItem* get_kNN_item(kNNGraph * kNN, int i_list, int i_k);

#endif
