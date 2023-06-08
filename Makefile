CPP := g++
CC := gcc
LINK := g++
mkfile_dir := $(dir $(abspath $(lastword $(MAKEFILE_LIST))))

CSOURCES= $(wildcard src/*.c)  
CPPSOURCES= $(wildcard src/*.cpp)
OBJS=$(CSOURCES:.c=.c.o)  $(CPPSOURCES:.cpp=.cpp.o)
COMMONOBJS= lib/libhts.a lib/libboost_iostreams.a lib/libtcmalloc_minimal.a
# Warnings is included in WarningsAsErrors to make sure that the warning is enabled.
Warnings=-Wreturn-type -Warray-bounds -Wmaybe-uninitialized -Waddress
WarningsAsErrors=$(Warnings) -Werror=return-type -Werror=array-bounds -Werror=address
CFLAGS=  $(WarningsAsErrors) -Wno-unused-function
CPPFLAGS=  $(WarningsAsErrors) -Wno-unused-function -std=c++11

RELEASE_FLAGS= -O3 -g

# Includes
INCLUDES = -Iinc/
#
# Common flags
COMMONFLAGS += $(INCLUDES)

CXXFLAGS += $(COMMONFLAGS)
CFLAGS += $(COMMONFLAGS)
CPPFLAGS += $(COMMONFLAGS)
COMMONLIBS= -Llib/ -lz -lm -lpthread -lbz2 -lcurl -lcrypto -llzma -fopenmp

#LIBS += $(COMMONLIBS) -ltcmalloc
LIBS += $(COMMONLIBS)

TARGET = MuSE
LINKLINE = $(LINK)  -O3 -o $(TARGET) $(OBJS) $(COMMONOBJS) $(MATCHOBJS) $(LIBS)

#all:
all: $(TARGET) 
.SUFFIXES: .c .cpp .o

%.c.o: %.c
	$(CC) $(CFLAGS) $(RELEASE_FLAGS) -c $< -o $@

%.cpp.o: %.cpp
	$(CPP) $(CPPFLAGS) $(RELEASE_FLAGS) -c $< -o $@

$(TARGET): $(OBJS) Makefile
	$(LINKLINE)

clean:
	rm -f $(OBJS) $(TARGET)
