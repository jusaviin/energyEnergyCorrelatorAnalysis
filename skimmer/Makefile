PROGRAM       = megaSkimmerPp

version       = skimming
CXX           = g++
CXXFLAGS      = -g -Wall -Wno-deprecated -D$(version) 
LD            = g++
LDFLAGS       = -O2
SOFLAGS       = -shared
#############################################
# -bind_at_load helps to remove linker error
############################################
CXXFLAGS += $(shell root-config --cflags)
LDFLAGS  += $(shell root-config --libs)

# Put to hdrsdict header files from all classes that inherit TObject
# HDRSDICT = 
        
# Use the following form if you have classes inherit TObject
# HDRS += $(HDRSDICT) src/Class.h ... nanoDict.h       
HDRS += TrackingEfficiencyInterface.h trackingEfficiency2017pp.h

SRCS = $(HDRS:.h=.cxx)
OBJS = $(HDRS:.h=.o)

all:            $(PROGRAM)

$(PROGRAM):     $(OBJS) $(PROGRAM).cxx
		@echo "Linking $(PROGRAM) ..."
		$(CXX) -lEG -L$(PWD) $(PROGRAM).cxx $(CXXFLAGS) $(OBJS) $(LDFLAGS) -o $(PROGRAM)
		@echo "done"

%.cxx:

%: %.cxx
#  commands to execute (built-in):
	$(LINK.cc) $^ $(CXXFLAGS) $(LOADLIBES) $(LDLIBS) -o $@

%.o: %.cxx %.h
#  commands to execute (built-in):
	$(COMPILE.cc) $(OUTPUT_OPTION) $<

# If dictionaries built, need to clean also them: *Dict*
clean:
		rm -rf $(OBJS) $(PROGRAM).o *.dSYM $(PROGRAM)

cl:  clean $(PROGRAM)

# Dictionary is needed for all classes inheriting TObject from root
# nanoDict.cc: $(HDRSDICT)
#		@echo "Generating dictionary ..."
#		@rm -f nanoDict.cc nanoDict.hh nanoDict.h
#		@rootcling nanoDict.cc -c -D$(version) $(HDRSDICT)
