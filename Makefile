PROGRAM       = toyFlow

version       = toka
CXX           = g++
CXXFLAGS      = -Ofast -Wall -g -Wno-deprecated -D$(version)

LD            = g++
LDFLAGS       = -Ofast
SOFLAGS       = -shared
CXXFLAGS += $(shell root-config --cflags)
LDFLAGS  = $(shell root-config --libs)

HDRSDICT =  JToyMCTrack.h JHistos.h

HDRS    += $(HDRSDICT) nanoDict.h


SRCS = $(HDRS:.h=.cxx)
OBJS = $(HDRS:.h=.o)

all:            $(PROGRAM) $(PROGRAM2)

$(PROGRAM):     $(OBJS) $(PROGRAM).C
		@echo "Linking $(PROGRAM) ..."
		$(CXX) -lEG -lPhysics -L$(PWD) $(PROGRAM).C $(CXXFLAGS) $(OBJS) $(LDFLAGS) $(FFTWINC) -o $(PROGRAM)
		@echo "done"

%.cxx:

%: %.cxx
#  commands to execute (built-in):
	$(LINK.cc) $^ $(CXXFLAGS) $(LOADLIBES) $(LDLIBS)  -o $@

%.o: %.cxx %.h
#  commands to execute (built-in):
	$(COMPILE.cc) $(OUTPUT_OPTION) $(FFTWINC) $<




clean:
	rm -f $(PROGRAM) *.o nanoDict* *_h* *_C* *_c*

nanoDict.cc: $(HDRSDICT)
		@echo "Generating dictionary ..."
		@rm -f nanoDict.cc nanoDict.hh nanoDict.h
		@rootcint nanoDict.cc -c -D$(version) $(HDRSDICT)
