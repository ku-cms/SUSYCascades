# /cvmfs/cms.cern.ch/slc7_amd64_gcc700/cms/cmssw/CMSSW_10_6_5/external/slc7_amd64_gcc700/
# /cvmfs/cms.cern.ch/el9_amd64_gcc12/cms/cmssw/CMSSW_13_3_1/external/el9_amd64_gcc12/
# /cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/lhapdf/6.2.1-pafccj3/bin/lhapdf-config
# /cvmfs/cms.cern.ch/el9_amd64_gcc12/external/lhapdf/6.4.0-52852f9a177b8e8b5b72e2ae6b1327b6/bin/lhapdf-config
# /cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/lwtnn/2.4-gnimlf3/
# /cvmfs/cms.cern.ch/el9_amd64_gcc12/external/lwtnn/2.14.1-f5d88789e43e4522aa674449efbea1c5
ROOTCFLAGS  = $(shell root-config --cflags)
ROOTGLIBS   = $(shell root-config --glibs)

RFCFLAGS    = $(shell restframes-config --cxxflags)
RFGLIBS     = $(shell restframes-config --libs)

CXX         = g++

CXXFLAGS       = -fPIC $(filter-out -stdlib=libc++ -pthread , $(ROOTCFLAGS))
CXXFLAGS      += $(filter-out -stdlib=libc++ -pthread , $(RFCFLAGS))

GLIBS          = $(filter-out -stdlib=libc++ -pthread , $(ROOTGLIBS))
GLIBS         += $(filter-out -stdlib=libc++ -pthread , $(RFGLIBS))
GLIBS         += -lRooFit -lRooFitCore

SGLIBS          = $(filter-out -stdlib=libc++ -pthread , $(ROOTGLIBS))
SGLIBS         += $(filter-out -stdlib=libc++ -pthread , $(RFGLIBS))
SGLIBS         += -lRooFit -lRooFitCore

INCLUDEDIR       = ./include/
SRCDIR           = ./src/
CXX	            += -I$(INCLUDEDIR) -I.
OUTOBJ	         = ./obj/
OUTOBJ_CMSSW	 = ./obj_cmssw/

INCLUDEDIR_CMSSW  = ./include_cmssw/
SRCDIR_CMSSW      = ./src_cmssw/

SCXX   = $(CXX)

#LHAPDF stuff is cmssw specific
cmssw: LHAPDFCFLAGS = $(shell /cvmfs/cms.cern.ch/el9_amd64_gcc12/external/lhapdf/6.4.0-52852f9a177b8e8b5b72e2ae6b1327b6/bin/lhapdf-config --cflags --ldflags)
cmssw: LHAPDFGLIBS = $(shell /cvmfs/cms.cern.ch/el9_amd64_gcc12/external/lhapdf/6.4.0-52852f9a177b8e8b5b72e2ae6b1327b6/bin/lhapdf-config --libs)
cmssw: CXXFLAGS      += $(filter-out -stdlib=libc++ -pthread , $(LHAPDFCFLAGS))
cmssw: GLIBS         += $(filter-out -stdlib=libc++ -pthread , $(LHAPDFGLIBS))
cmssw: CXX += -I$(INCLUDEDIR_CMSSW)

# Needed for correctionlib
# PYTHON_INCLUDE = $(shell python3 -c "from sysconfig import get_paths as gp; print(gp()['include'])")
# PYTHON_LIB = $(shell python3 -c "import sysconfig; print(sysconfig.get_config_var('LIBDIR'))")
# CORRECTIONLIB_DIR = /cvmfs/cms.cern.ch/el9_amd64_gcc12/external/py3-correctionlib/2.2.2-4f091fd2adcff55f05f1d262b3254c25/lib/python3.9/site-packages/correctionlib/
# INCLUDEDIR_CLIB = -I$(PYTHON_INCLUDE) -I$(CORRECTIONLIB_DIR)/include
# LIBDIR_CLIB = -L$(PYTHON_LIB) -L$(CORRECTIONLIB_DIR)/lib
CORRECTION_FLAGS = $(shell correction config --cflags --ldflags --rpath)

# For correctionlib
# cmssw : CXXFLAGS += $(INCLUDEDIR_CLIB)
# cmssw : GLIBS += $(LIBDIR_CLIB)
cmssw : CXXFLAGS += $(CORRECTION_FLAGS)
cmssw : GLIBS += $(CORRECTION_FLAGS)

CC_FILES := $(wildcard src/*.cc)
HH_FILES := $(wildcard include/*.hh)

CC_FILES_CMSSW := $(wildcard src_cmssw/*.cc)
cmssw: HH_FILES += $(wildcard include_cmssw/*.hh)

OBJ_FILES := $(addprefix $(OUTOBJ),$(notdir $(CC_FILES:.cc=.o)))
OBJ_FILES_CMSSW := $(addprefix $(OUTOBJ_CMSSW),$(notdir $(CC_FILES_CMSSW:.cc=.o)))

SOBJ_FILES = $(filter-out ./obj/AnalysisBase.o ./obj/SVDiscrTool.o ./obj/ReducedNtuple.o ./obj/NtupleBase.o ./obj/LHETool.o, $(OBJ_FILES))

all : GLIBS += -L/cvmfs/cms.cern.ch/el9_amd64_gcc12/external/lwtnn/2.14.1-f5d88789e43e4522aa674449efbea1c5/lib -llwtnn
all : CXX   += -I/cvmfs/cms.cern.ch/el9_amd64_gcc12/external/lwtnn/2.14.1-f5d88789e43e4522aa674449efbea1c5/include/

local : GLIBS += -L/Users/christopherrogan/GitHub/lwtnn/lib -llwtnn
local : CXX   += -I/Users/christopherrogan/GitHub/lwtnn/include 

# cmssw : GLIBS += -L../../lib/slc7_amd64_gcc700 -lCombineHarvesterCombinePdfs -lHiggsAnalysisCombinedLimit -lCombineHarvesterCombineTools
cmssw : GLIBS += -L/cvmfs/cms.cern.ch/el9_amd64_gcc12/external/lwtnn/2.14.1-f5d88789e43e4522aa674449efbea1c5/lib -llwtnn
cmssw : GLIBS += -L/cvmfs/cms.cern.ch/el9_amd64_gcc12/cms/cmssw/CMSSW_13_3_1/external/el9_amd64_gcc12/lib/ -lvdt -lboost_program_options -lboost_filesystem -lboost_regex -lboost_system

cmssw : CXX   += -I../. -I/cvmfs/cms.cern.ch/el9_amd64_gcc12/external/boost/1.80.0-e6c5c62cb3bb17bd5a918a7756ae6c58/include/
cmssw : CXX   += -I/cvmfs/cms.cern.ch/el9_amd64_gcc12/external/lwtnn/2.14.1-f5d88789e43e4522aa674449efbea1c5/include/
# cmssw : CXX   += -I../../src/HiggsAnalysis/CombinedLimit/interface/
cmssw : CXX   += -D_CMSSW_

locallib : GLIBS += -L/Users/christopherrogan/GitHub/lwtnn/lib -llwtnn
locallib : CXX   += -I/Users/christopherrogan/GitHub/lwtnn/include


all: alltargets lib

#cmssw: alltargets lib BuildFit.x
cmssw: alltargets lib

local: alltargets lib

locallib: lib

lib: lib/libKUEWKino.so

#alltargets: MakeReducedNtuple_NANO.x EventCountPlot.x MakeEventCount_NANO.x BuildFitInput.x Condor_Plot_1D_NANO.x BuildPlotInput.x BuildFitShapes.x BuildFitInputCondor.x BuildPlotInputCondor.x BuildFitCondor.x
#alltargets: MakeReducedNtuple_NANO.x MakeEventCount_NANO.x Condor_Plot_1D_NANO.x Submit_Plot_1D_NANO.x BuildFitInput.x BuildFitInputCondor.x
alltargets: MakeReducedNtuple_NANO.x MakeEventCount_NANO.x Condor_Plot_1D_NANO.x Submit_Plot_1D_NANO.x

EventCountPlot.x:  $(SRCDIR)EventCountPlot.C $(OBJ_FILES) $(HH_FILES)
	$(CXX) $(CXXFLAGS) -o EventCountPlot.x $(OUTOBJ)/*.o $(GLIBS) $ $<
	touch EventCountPlot.x

MakeReducedNtuple.x:  $(SRCDIR)MakeReducedNtuple.C $(OBJ_FILES) $(HH_FILES)
	$(CXX) $(CXXFLAGS) -o MakeReducedNtuple.x $(OUTOBJ)/*.o $(GLIBS) $ $<
	touch MakeReducedNtuple.x

MakeEventCount.x:  $(SRCDIR)MakeEventCount.C $(OBJ_FILES) $(HH_FILES)
	$(CXX) $(CXXFLAGS) -o MakeEventCount.x $(OUTOBJ)/*.o $(GLIBS) $ $<
	touch MakeEventCount.x

MakeReducedNtuple_NANO.x:  $(SRCDIR)MakeReducedNtuple_NANO.C $(OBJ_FILES) $(HH_FILES)
	$(CXX) $(CXXFLAGS) -o MakeReducedNtuple_NANO.x $(OUTOBJ)/*.o $(GLIBS) $ $<
	touch MakeReducedNtuple_NANO.x

MakeEventCount_NANO.x:  $(SRCDIR)MakeEventCount_NANO.C $(OBJ_FILES) $(HH_FILES)
	$(CXX) $(CXXFLAGS) -o MakeEventCount_NANO.x $(OUTOBJ)/*.o $(GLIBS) $ $<
	touch MakeEventCount_NANO.x

BuildFit.x:  $(SRCDIR)BuildFit.C $(OBJ_FILES) $(OBJ_FILES_CMSSW) $(HH_FILES)
	$(CXX) $(CXXFLAGS) -o BuildFit.x $(OUTOBJ)/*.o $(OUTOBJ_CMSSW)/*.o $(GLIBS) $ $<
	touch BuildFit.x

BuildFitInput.x:  $(SRCDIR)BuildFitInput.C $(OBJ_FILES) $(HH_FILES)
	$(CXX) $(CXXFLAGS) -o BuildFitInput.x $(OUTOBJ)/*.o $(GLIBS) $ $<
	touch BuildFitInput.x

BuildPlotInput.x:  $(SRCDIR)BuildPlotInput.C $(SOBJ_FILES) $(HH_FILES)
	$(CXX) $(CXXFLAGS) -o BuildPlotInput.x $(SOBJ_FILES) $(SGLIBS) $ $<
	touch BuildPlotInput.x

Condor_Plot_1D_NANO.x:  $(SRCDIR)Condor_Plot_1D_NANO.C $(OBJ_FILES) $(HH_FILES)
	$(CXX) $(CXXFLAGS) -o Condor_Plot_1D_NANO.x $(OUTOBJ)/*.o $(GLIBS) $ $<
	touch Condor_Plot_1D_NANO.x

Submit_Plot_1D_NANO.x:  $(SRCDIR)SubmitCondor_Plot_1D_NANO.C $(OBJ_FILES) $(HH_FILES)
	$(CXX) $(CXXFLAGS) -o Submit_Plot_1D_NANO.x $(OUTOBJ)/*.o $(GLIBS) $ $<
	touch Submit_Plot_1D_NANO.x

BuildFitShapes.x:  $(SRCDIR)BuildFitShapes.C $(OBJ_FILES) $(HH_FILES)
	$(CXX) $(CXXFLAGS) -o BuildFitShapes.x $(OUTOBJ)/*.o $(GLIBS) $ $<
	touch BuildFitShapes.x

BuildFitInputCondor.x:  $(SRCDIR)BuildFitInputCondor.C $(OBJ_FILES) $(HH_FILES)
	$(CXX) $(CXXFLAGS) -o BuildFitInputCondor.x $(OUTOBJ)/*.o $(GLIBS) $ $<
	touch BuildFitInputCondor.x

BuildPlotInputCondor.x:  $(SRCDIR)BuildPlotInputCondor.C $(OBJ_FILES) $(HH_FILES)
	$(CXX) $(CXXFLAGS) -o BuildPlotInputCondor.x $(OUTOBJ)/*.o $(GLIBS) $ $<
	touch BuildPlotInputCondor.x

BuildFitCondor.x:  $(SRCDIR)BuildFitCondor.C $(OBJ_FILES) $(HH_FILES)
	$(CXX) $(CXXFLAGS) -o BuildFitCondor.x $(OUTOBJ)/*.o $(GLIBS) $ $<
	touch BuildFitCondor.x

lib/libKUEWKino.so: $(SOBJ_FILES)
	mkdir -p lib
	$(SCXX) -shared -o lib/libKUEWKino.so $(SOBJ_FILES) $(SGLIBS)
	touch lib/libKUEWKino.so

$(OUTOBJ)%.o: src/%.cc include/%.hh
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(OUTOBJ_CMSSW)%.o: src_cmssw/%.cc include_cmssw/%.hh
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	rm -f $(OUTOBJ)*.o
	rm -f $(OUTOBJ_CMSSW)*.o
	rm -f *.x
	rm -f AutoDict*
	rm -f lib/libKUEWKino.so
