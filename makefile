######################
ARCH          = $(shell root-config --arch)

SrcDir        = src
LibDir        = lib
ObjSuf        = o
SrcSuf        = cxx
FcnSuf        = cxx
HdrSuf        = h
ExeSuf        =
ifeq ($(ARCH),macosx)
DllSuf        = dylib
else
DllSuf        = so
endif
OutPutOpt     = -o 

ROOTCFLAGS   := $(shell root-config --cflags)
ROOTLIBS     := $(shell root-config --libs) -lMinuit
ROOTGLIBS    := $(shell root-config --glibs)


CXX           = g++
CXXFLAGS      = -O  -fPIC
LD            = g++
LDFLAGS       = -O  
ifeq ($(ARCH),macosx)
SOFLAGS       = -dynamiclib -single_module -undefined dynamic_lookup
else
SOFLAGS       = -shared
endif


CXXFLAGS     += $(ROOTCFLAGS)
CXXFLAGS     += -I$(SrcDir)
LIBS          = $(ROOTLIBS) $(SYSLIBS)
GLIBS         = $(ROOTGLIBS) $(SYSLIBS)

#------------------------------------------------------------------------------

CLASSES			= Pidrix
FNCS        	= initialization plotting updating

#------------------------------------------------------------------------------

SEP_CLASSESO    = $(addsuffix .$(ObjSuf), $(CLASSES))
SEP_CLASSESH    = $(addprefix $(SrcDir)/, $(addsuffix .$(HdrSuf), $(CLASSES)))
SEP_CLASSESS    = $(addprefix $(SrcDir)/, $(addsuffix .$(SrcSuf), $(CLASSES)))

FNCSS		    = $(addprefix $(SrcDir)/, $(addsuffix .$(FcnSuf), $(FNCS)))
FNCSH		    = $(addprefix $(SrcDir)/, $(addsuffix .$(HdrSuf), $(FNCS)))
FNCSO		    = $(addsuffix .$(ObjSuf), $(FNCS))

CLASSDICTO		= ClassesDict.$(ObjSuf)
CLASSDICTS		= $(SrcDir)/ClassesDict.$(SrcSuf)

PIDRIXSO		= $(LibDir)/libPidrix.$(DllSuf)

OBJS          	= $(SEP_CLASSESO) $(CLASSDICTO) $(FNCSO)
PROGRAMS      	= $(PIDRIXSO)
 
#------------------------------------------------------------------------------

.SUFFIXES: .$(SrcSuf) .$(ObjSuf) .$(DllSuf)
.PHONY:   Classes Functions Pidrix all clean distclean

ifeq ($(ARCH),macosx)
all:            $(PROGRAMS)
		@ln -sf lib/libPidrix.$(DllSuf) \
		        lib/libPidrix.so
else
all:            $(PROGRAMS)
endif

clean:
		@rm -f $(OBJS) $(SrcDir)/*Dict.* core

distclean:      clean
		@rm -f $(PROGRAMS) *Dict.* *.def *.exp \
		   *.root *.ps *.so .def so_locations
		@rm -rf cxx_repository

###

Pidrix:	$(PIDRIXSO)	

$(PIDRIXSO):	$(SEP_CLASSESO) $(CLASSDICTO) $(FNCSO)
	@echo "linking pidrix"
	$(LD) $(SOFLAGS) $(LIBS) $(LDFLAGS) $^ $(OutPutOpt) $@

Classes:	$(SEP_CLASSESO) $(CLASSDICTO)

$(SEP_CLASSESO): 	$(SEP_CLASSESS) $(SEP_CLASSESH) 
	@echo "compiling with $(CXXFLAGS)"
	$(CXX) $(CXXFLAGS) -c $(SEP_CLASSESS) 

$(CLASSDICTO): $(CLASSDICTS)
	@echo "compiling with $(CXXFLAGS)"
	$(CXX) $(CXXFLAGS) -c $(CLASSDICTS) 

$(CLASSDICTS): $(SEP_CLASSESH) $(FNCSH) src/PidrixClassesLinkDef.$(HdrSuf)
	@echo "Generating dictionary $@..."
	@cd $(SrcDir); echo rootcint -f $(@:$(SrcDir)/%=%) -c $(^:$(SrcDir)/%=%)
	@cd $(SrcDir); rootcint -f $(@:$(SrcDir)/%=%) -c $(^:$(SrcDir)/%=%)

Functions:	$(FNCSSO)	

$(FNCSSO):	$(FNCSO)
	@echo "linking functions"
	$(LD) $(SOFLAGS) $(LIBS) $(LDFLAGS) $^ $(OutPutOpt) $@

$(FNCSO): 	$(FNCSS) $(SEP_CLASSESH) 
	@echo "compiling with $(CXXFLAGS)"
	$(CXX) $(CXXFLAGS) -c $(FNCSS) 

.$(SrcSuf).$(ObjSuf):
	$(CXX) $(CXXFLAGS) -c $<
