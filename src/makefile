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

GADEP= $(OBJECTS)sa.o  $(OBJECTS)sortcb.o $(OBJECTS)sort.o $(OBJECTS)solution.o $(OBJECTS)heap.o $(OBJECTS)stack.o $(OBJECTS)pnn.o

KMEANSOBJ = kmeans/inits.o kmeans/kmeans.o
RSOBJ = rs/rs.o

BINS=cbkm cbrs cbga

OPT     = -O2 $(DEBUG) -std=c99  -Wall  -I. -I$(MODULES) -Ikmeans -Irs -Iga 


.PHONY : all clean
all: $(BINS)

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

