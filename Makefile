C_SRC    = bgzf.c faidx.c index.c kprobaln.c kstring.c razf.c
CPP_SRC  = MuSE.cpp bam_supplement.cpp bam.cpp fet.cpp sample.cpp
BIN      = MuSE
C_OBJ    = $(C_SRC:.c=.o)
CPP_OBJ  = $(CPP_SRC:.cpp=.o)
CPP      = g++
CPPFLAGS = -O3 -w -g -Wall -D_FILE_OFFSET_BITS=64 -D_LARGEFILE64_SOURCE

$(BIN): $(C_OBJ) $(CPP_OBJ)
	$(CPP) $(C_OBJ) $(CPP_OBJ) -o $(BIN) -lm -lz
	
bam.o : bam.cpp
	$(CPP) $(CPPFLAGS) -c bam.cpp -o bam.o
	
MuSE.o : MuSE.cpp
	$(CPP) $(CPPFLAGS) -c MuSE.cpp -o MuSE.o

bam_supplement.o : bam_supplement.cpp
	$(CPP) $(CPPFLAGS) -c bam_supplement.cpp -o bam_supplement.o

bgzf.o : bgzf.c
	$(CPP) $(CPPFLAGS) -c bgzf.c -o bgzf.o
	
faidx.o : faidx.c
	$(CPP) $(CPPFLAGS) -c faidx.c -o faidx.o

kprobaln.o : kprobaln.c
	$(CPP) $(CPPFLAGS) -c kprobaln.c -o kprobaln.o

kstring.o : kstring.c
	$(CPP) $(CPPFLAGS) -c kstring.c -o kstring.o
	
razf.o : razf.c
	$(CPP) $(CPPFLAGS) -c razf.c -o razf.o

sample.o : sample.cpp
	$(CPP) $(CPPFLAGS) -c sample.cpp -o sample.o

fet.o : fet.cpp
	$(CPP) $(CPPFLAGS) -c fet.cpp -o fet.o

index.o : index.c
	$(CPP) $(CPPFLAGS) -c index.c -o index.o

clean: 
	rm -f $(C_OBJ) $(CPP_OBJ) $(BIN)

