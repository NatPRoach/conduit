
AR=ar rc

TARGETS=poa liblpo.a poa_doc libbflag.a

# align_score.c CAN BE USED TO ADD CUSTOMIZED SCORING FUNCTIONS
OBJECTS= \
	align_score.o \
	main.o


LIBOBJECTS= \
	black_flag.o \
	seq_util.o \
	fasta_format.o \
	align_lpo2.o \
	align_lpo_po2.o \
	buildup_lpo.o \
	lpo.o \
	lpo_format.o \
	create_seq.o \
	numeric_data.o \
	stringptr.o


ifeq ($(CC),)
CC := gcc
endif
# CC = gcc
#CFLAGS= -g -ansi-strict -W -Wall -DUSE_WEIGHTED_LINKS -DUSE_PROJECT_HEADER -I.
CFLAGS= -g -O2 -DUSE_WEIGHTED_LINKS -DUSE_PROJECT_HEADER -I.
# -I$(HOME)/lib/include
# -DREPORT_MAX_ALLOC

# NB: LIBRARY MUST FOLLOW OBJECTS OR LINK FAILS WITH UNRESOLVED REFERENCES!!
align_score.o: align_score.c liblpo.a
	$(CC) -c align_score.c -o $@

# poa: $(OBJECTS) liblpo.a
# 	$(CC) -o $@ $(OBJECTS) -lm liblpo.a

clean:
	rm -f $(OBJECTS) $(LIBOBJECTS) $(TARGETS)

liblpo.a: $(LIBOBJECTS)
	rm -f $@
	$(AR) $@ $(LIBOBJECTS)
	ranlib $@

# clean2:
# 	rm -f $(LIBOBJECTS) poa poa_doc libbflag.a

what:
	@echo poa: partial-order based sequence alignment program
	@echo liblpo.a: partial-order alignment and utilities function library
