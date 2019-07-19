
#define bool int
// #include "options.h"
#include "knngraph.h"
//#include "knng_dataset.h"

// struct options g_options;


kNNGraph * init_kNNGraph(int N, int K, int maxK) {
    kNNGraph* kNN = (kNNGraph*) safemalloc(sizeof(kNNGraph));
    kNN->list = (kNNList*) safemalloc(sizeof(kNNList)*N);
    kNN->size = N;
    kNN->k = K;

    kNNList * curlist = kNN->list;
    for (int i = 0; i < N; i++)
    {
        curlist->items = (kNNItem*) safemalloc(sizeof(kNNItem)*maxK);
        /*curlist->dist = (FLOAT*) safemalloc(sizeof(FLOAT)*K);*/
        curlist->size = 0;
        curlist->max_dist = MAX_FLOAT;
        curlist->is_exact = false;
        curlist++;
    }
    return kNN;
}

void free_kNNGraph(kNNGraph * kNN) {
    int i;
    for (i = 0; i < kNN->size; i++)
    {
        free(kNN->list[i].items);
    }
    free(kNN);
    free(kNN->list);
}


/*kNNGraph * init_kNNGraph(int N, int K ) {*/
/*return init_kNNGraph(N, K, K);*/
/*}*/

void debug_graph(kNNGraph* knng) {
    printf("knng->k: %d\n",knng->k);

    for (int i_row=0; i_row < 10; i_row++) {
        for (int j = 0; j < knng->k; ++j) {
            printf("%d ",knng->list[i_row].items[j].id);
        }
        printf("\n");
    }
}

int updatekNN(kNNGraph* kNN, UINT p1, UINT p2, FLOAT dist)
{

#ifdef SANITYCHECK
if(p1 == p2 && kNN->format != RANDOM_SAMPLED_BRUTEFORCE) {
    printf("p1=%u=p2\n",p1);
    terminal_error("p1=p2");
}
if(p1 >= kNN->size) {
    printf("p1 = %u >= kNN->size\n",p1);
    terminal_error(" ");
}
#endif

    kNNList* kl = &kNN->list[p1];
    kNNItem* ki = kl->items;

    if(kl->max_dist > dist || kl->size < kNN->k) {
        int i = 0;
        for(i=0; i <= kl->size; i++) {
            /*if(*datap == p2) { return 0 ;} //TODO:??*/
            if(ki->id == p2) { return 0 ;} //TODO:??
            /*if(ki->id == p2) {num_same_edg++; return 0 ;} //TODO:??*/

            if (ki->dist > dist || i == kl->size) {
                int moveamount = kNN->k - i -1;
                /*int moveamount = kNN->size - i -1;*/
                if(moveamount > 0 ) { //TODO: needed?
                    // Move from ki to ki+1
                    memmove(ki+1,ki,moveamount*sizeof(kNNItem));
                }
                /**datap = p2;*/
                /**distp = dist;*/
                ki->id = p2;
                ki->dist = dist;
                ki->new_item = true; //TODO: not needed in all search types

                if(kl->size < kNN->k) { kl->size++;}
                kl->max_dist = kl->items[kl->size -1].dist; //TODO: optimize?
                break;
            }
            ki++;

        }
        /*kNN->data[p1]->*/

        return 1; // Did update
    }
    else {
        return 0; // Did not update
    }
}

/*inline */
UINT get_kNN_item_id(kNNGraph * kNN, int i_list, int i_k) {
#ifdef SANITYCHECK
    if (i_list < 0 || i_list > kNN->size || i_k >= kNN->k
            || i_k >= kNN->list[i_list].size) {
        printf("%d x %d\n",i_list,kNN->size);
        std::raise(SIGINT);
        terminal_error("get_kNN_item_id: invalid i_list or i_k params");
    }
#endif
    return kNN->list[i_list].items[i_k].id;
}




/*inline void set_kNN_val(kNNGraph * kNN, int i_list, int i_k, UINT id) {*/
/*kNNItem* ki = &kNN->list[i_list].items[i_k];*/
/*ki->id = id;*/
/*}*/

inline void set_kNN_val(kNNGraph * kNN, int i_list, int i_k, UINT id, FLOAT dist, bool new_item) {
    //TODO: &
    kNNItem ki = kNN->list[i_list].items[i_k];
    ki.id = id;
    ki.dist = dist;
    ki.new_item = new_item;
}

inline kNNList* get_kNN_list(kNNGraph * kNN, int i_list) {
#ifdef SANITYCHECK
    if (i_list < 0 || i_list > kNN->size) {
        printf("%d x %d\n",i_list,kNN->size);
        std::raise(SIGINT);
        terminal_error("invalid i_list");
    }
#endif
    return &(kNN->list[i_list]);
}

inline void set_kNN_id(kNNGraph * kNN, UINT i_list, UINT id) {
#ifdef SANITYCHECK
    if (i_list < 0 || i_list > kNN->size ) {
        printf("%d x %d\n",i_list,kNN->size);
        std::raise(SIGINT);
        terminal_error("get_kNN_item_id: invalid i_list or i_k params");
    }
#endif
    kNN->list[i_list].id = id;
}


inline kNNItem* get_kNN_item(kNNGraph * kNN, int i_list, int i_k) {
    return &kNN->list[i_list].items[i_k];
}




