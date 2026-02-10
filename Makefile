# Usage:
# make        # compile all binary
# make clean  # remove ALL binaries and objects

.PHONY = all clean
CC = gcc
LD = ld
TARGET = PTracker

HDF5_INCLUDE_DIR = /home/gyl/anaconda3/include
HDF5_LIB_DIR = /home/gyl/anaconda3/lib

LIBCONFIG_INCLUDE_DIR = /home/gyl/anaconda3/pkgs/glib-2.69.1-he621ea3_2/lib/glib-2.0/include
LIBCONFIG_LIB_DIR = /home/gyl/anaconda3/pkgs/glib-2.69.1-he621ea3_2/lib/glib-2.0/lib


CFLAGS = -I$(HDF5_INCLUDE_DIR) -I$(LIBCONFIG_INCLUDE_DIR) -g -fopenmp -DDEBUG
LDFLAGS = -L$(HDF5_LIB_DIR) -L$(LIBCONFIG_LIB_DIR) -lm -lconfig -lhdf5 -lhdf5_hl -fopenmp
SRCS = $(wildcard *.c)
BINS = $(SRCS:%.c=%)

direct:
	${CC} *.c $(CFLAGS) ${LDFLAGS} -o ${TARGET}

all: ${BINS}  

show:
	@echo ${SRCS}
	@echo ${BINS}

%: %.c
	@echo "Compiling..."
	${CC} $(CFLAGS) -c $? -o $@.o

link:
	@echo "Linking..."
	${LD} ${LDFLAGS} -o ${TARGET} *.o

clean:
	@echo "Cleaning up..."
	rm -rvf *.o

