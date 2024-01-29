# C Compiler for final compilation
CC        = gcc

# CUDA Compiler
NVCC      = nvcc

# Flags to configure the field theory code
CODEFLAGS = -DVERSION="\"1.0.0\"" 

# Flags to the C compiler
CFLAGS    = -w -DBUILD_INFO="\"`date +%d-%b-%Y`\""

all: src/dpd.o
	$(CC) $(CFLAGS) $(CODEFLAGS) $(INCS) $(LPATH) -o dpd dpd.o $(LIBS) 
	rm *.o

src/dpd.o:
	$(CC) $(CFLAGS) $(CODEFLAGS) $(INCS) -c src/dpd.c
clean:
	rm -f *.o *.dat *.xyz dpd
