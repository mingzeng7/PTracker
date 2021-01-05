# Usage:
# make        # compile all binary
# make clean  # remove ALL binaries and objects

.PHONY = all clean
CC = gcc
LD = ld
TARGET = PTracker
CFLAGS = -I$(HOME)/.local/include -I$(HOME)/.local/hdf5-1.10.7/include
LDFLAGS = -L$(HOME)/.local/lib -L$(HOME)/.local/hdf5-1.10.7/lib -lm -lconfig -lhdf5 -lhdf5_hl
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

