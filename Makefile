
#CFLAGS = $(shell gsl-config --cflags) -std=c++0x 
CFLAGS =-O3 -I include
#LIBS   = $(shell gsl-config --libs)

CC     = g++

BDIR   = bin
ODIR   = obj
SDIR   = src
EDIR   = EXEs
OS    := $(shell uname)

ifeq ($(OS),Linux)
	_OBJS = $(shell ls $(SDIR)/*.cpp | xargs -n 1 basename | sed -r 's/(\.cc|.cpp)/.o/')
else
	_OBJS = $(shell ls $(SDIR)/*.cpp | xargs -n 1 basename | sed -E 's/(\.cc|.cpp)/.o/')
endif

_EXECUTABLES = 	rd_trap_rewgt_1D \
		rd_pair_oned_rewgt \
		rd_trap_Events_rewgt1D \

OBJS = $(patsubst %,$(ODIR)/%,$(_OBJS))

EXECUTABLES = $(patsubst %,$(BDIR)/%,$(_EXECUTABLES))


_SOURCES = ${_EXECUTABLES:=.cpp}
SOURCES = $(patsubst %,$(EDIR)/%,$(_SOURCES))

all: dirs $(EXECUTABLES)

$(ODIR)/%.o: $(SDIR)/%.cpp
	@echo "Compiling $<"
	$(CC) $(CFLAGS) $(CFLAGS2) -c $< -o $@  
	@echo "------------"

$(EXECUTABLES): $(OBJS)
	@echo "Compiling $(EDIR)/$(@F).cpp"
	$(CC) $(CFLAGS) $(CFLAGS2) -o $@ $(EDIR)/$(@F).cpp $(OBJS) $(LIBS)
	@echo "------------"

dirs: 
	mkdir -p bin
	mkdir -p obj

clean:
	rm -f $(ODIR)/*.o $(EXECUTABLES)


