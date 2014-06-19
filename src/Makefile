ARCH = linux

ROOTCFLAGS := $(shell root-config --cflags)
ROOTLIBS := $(shell root-config --libs)
ROOTGLIBS := $(shell root-config --glibs)

ifeq ($(ARCH),linux)
CXX = g++
CXXFLAGS = -O -Wno-deprecated -fPIC -m32 -fno-inline -Wno-write-strings
LD = g++
LDFLAGS = -O
SOFLAGS = -shared
INCLUDE = $(ROOTSYS)/include/
endif

CXXFLAGS += $(ROOTCFLAGS)
LIBS = $(ROOTLIBS) $(SYSLIBS)
GLIBS = $(ROOTGLIBS) $(SYSLIBS)

#-----------------------------------------------

CXXFILES = RunData13.cxx
HEADERS = RunData13.h

#-----------------------------------------------

All: RunData13.so

RunData13.so: LinkDef.h $(HEADERS) $(CXXFILES)
	@echo "generating dictionary $@..."
	@rootcint -f RunData13Dict.cxx -c $(HEADERS) LinkDef.h
	g++ $(CXXFILES) RunData13Dict.cxx -shared -o RunData13.so $(CXXFLAGS) $(GLIBS)

clean:
	@rm RunData13.so
