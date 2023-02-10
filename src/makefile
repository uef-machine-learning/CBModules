MODULES = modules/
OBJECTS = modules/
CC=gcc

COMMONDEP = $(OBJECTS)cb.o       \
          $(OBJECTS)file.o       \
          $(OBJECTS)interfc.o    \
          $(OBJECTS)memctrl.o    \
          $(OBJECTS)random.o     \
          $(OBJECTS)reporting.o  \
		  $(OBJECTS)textfile.o \
		  $(OBJECTS)fvec.o \

KMEANSDEP=$(OBJECTS)sortcb.o \
		  $(OBJECTS)sort.o \
		  $(OBJECTS)heap.o \
		  $(OBJECTS)stack.o \
		  $(OBJECTS)graph.o \
		  $(OBJECTS)knngraph.o \
		
		
SMDEP = $(OBJECTS)pnn.o \
	$(OBJECTS)split.o \
	$(OBJECTS)sort.o \
	$(OBJECTS)sortcb.o \
	$(OBJECTS)sa.o \
	$(OBJECTS)bintree.o \
	$(OBJECTS)heap.o \
	$(OBJECTS)stack.o 
	
SPLITDEP = $(OBJECTS)sa.o	\
          $(OBJECTS)sortcb.o \
          $(OBJECTS)sort.o  	\
          $(OBJECTS)split.o 	\
          $(OBJECTS)bintree.o 	\
          
ISDEP = $(OBJECTS)sa.o	\
          $(OBJECTS)sortcb.o \
          $(OBJECTS)sort.o  	\
          $(OBJECTS)solution.o 	\


GADEP= $(OBJECTS)sa.o  $(OBJECTS)sortcb.o $(OBJECTS)sort.o $(OBJECTS)solution.o $(OBJECTS)heap.o $(OBJECTS)stack.o $(OBJECTS)pnn.o

KMEANSOBJ = kmeans/inits.o kmeans/kmeans.o
RSOBJ = rs/rs.o


SMOBJ = splitmerge/ims.o
# splitmerge/cbims.o

BINS=cbkm cbrs cbga cbims cbsplit

OPT     = -O3 $(DEBUG) -std=c99  -Wall  -I. -I$(MODULES) -Isplit -Ikmeans -Irs -Iga -Isplitmerge -Isplit -Iis -Wno-unused-function -Wno-unused-but-set-variable -Wno-unused-variable


.PHONY : all clean
all: $(BINS)


cbis: $(COMMONDEP) $(ISDEP) is/cbis.o
	$(CC) -o $@ $(OPT) $+ -lm

cbims: $(COMMONDEP) $(SMDEP) $(SMOBJ) splitmerge/cbims.o
	$(CC) -o $@ $(OPT) $+ -lm
	
cbsplit: $(COMMONDEP) $(SMDEP) split/cbsplit.o
	$(CC) -o $@ $(OPT) $+ -lm

cbkm: $(COMMONDEP) $(KMEANSDEP) $(KMEANSOBJ) kmeans/cbkm.o
	$(CC) -o $@ $(OPT) $+ -lm

cbrs: $(COMMONDEP) $(KMEANSDEP) $(KMEANSOBJ) $(RSOBJ) rs/cbrs.c
	$(CC) -o $@ $(OPT) $+ -lm

cbga: $(COMMONDEP) $(GADEP)  ga/cbga.c
	$(CC) -o $@ $(OPT) $+ -lm
	

%.o: %.c
	$(CC) $(OPT) -c $(patsubst %.o,%.c,$@) -o $@

clean:
	rm  $(BINS) */*.o

