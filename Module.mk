# Module.mk for RooStats/HistFactory module
# Copyright (c) 2008 Rene Brun and Fons Rademakers
#
# Author: Kyle Cranmer

MODNAME      := histfactory
MODDIR       := $(ROOT_SRCDIR)/roofit/$(MODNAME)
MODDIRS      := $(MODDIR)/src
MODDIRI      := $(MODDIR)/inc
MODDIRC      := $(MODDIR)/config

HISTFACTORYDIR  := $(MODDIR)
HISTFACTORYDIRS := $(HISTFACTORYDIR)/src
HISTFACTORYDIRI := $(HISTFACTORYDIR)/inc

### prepareHIstFactory script
HF_PREPAREHISTFACTORY := bin/prepareHistFactory


##### tf_makeworkspace.exe #####

HF_MAKEWORKSPACEEXES   := $(MODDIRS)/hist2workspace.cxx
HF_MAKEWORKSPACEEXEO   := $(call stripsrc,$(HF_MAKEWORKSPACEEXES:.cxx=.o))
HF_MAKEWORKSPACEEXEDEP := $(HF_MAKEWORKSPACEEXEO:.o=.d)

HF_MAKEWORKSPACEEXE    := bin/hist2workspace$(EXEEXT)

#for fast (not working)
#HF_MAKEWORKSPACEEXE2S   := $(MODDIRS)/MakeModelAndMeasurementsFast.cxx
#HF_MAKEWORKSPACEEXE2O   := $(call stripsrc,$(HF_MAKEWORKSPACEEXE2S:.cxx=.o))
#HF_MAKEWORKSPACEEXE2O   := roofit/histfactory/src/ConfigParser.o roofit/histfactory/src/EstimateSummary.o roofit/histfactory/src/Helper.o roofit/histfactory/src/HistoToWorkspaceFactoryFast.o roofit/histfactory/src/LinInterpVar.o  roofit/histfactory/src/MakeModelAndMeasurementsFast.o roofit/histfactory/src/PiecewiseInterpolation.o
#HF_MAKEWORKSPACEEXE2DEP := $(HF_MAKEWORKSPACEEXE2O:.o=.d)

#HF_MAKEWORKSPACEEXE2    := bin/hist2workspaceFast$(EXEEXT)



ifeq ($(PLATFORM),win32)
HF_LIBS = $(HISTFACTORYLIBEXTRA) "$(ROOTSYS)/lib/libHistFactory.lib" 
ifeq ($(BUILDMATHMORE),yes)
HF_LIBS += "$(ROOTSYS)/lib/libMathMore.lib" 
endif
else
#for other platforms HISTFACTORYLIBEXTRA is not defined 
#need to copy from config/Makefile.depend
HF_LIBS = -Llib -lRooFit -lRooFitCore -lTree -lRIO -lMatrix \
          -lHist -lMathCore -lGraf -lGpad -lMinuit -lFoam \
          -lRooStats -lXMLParser
HF_LIBS += -lHistFactory 
ifeq ($(BUILDMATHMORE),yes)
HF_LIBS += -lMathMore
HF_LIBS += $(GSLLIBDIR) $(GSLLIBS)
endif
endif


HF_LIBSDEP := $(HISTFACTORYLIBDEP)


##### libHistFactory #####
HISTFACTORYL    := $(MODDIRI)/LinkDef.h
HISTFACTORYDS   := $(call stripsrc,$(MODDIRS)/G__HistFactory.cxx)
HISTFACTORYDO   := $(HISTFACTORYDS:.cxx=.o)
HISTFACTORYDH   := $(HISTFACTORYDS:.cxx=.h)

HISTFACTORYH    := $(filter-out $(MODDIRI)/LinkDef%,$(wildcard $(MODDIRI)/RooStats/HistFactory/*.h))
HISTFACTORYS    := $(filter-out $(MODDIRS)/G__%,$(wildcard $(MODDIRS)/*.cxx))
HISTFACTORYO    := $(filter-out roofit/histfactory/src/hist2workspace.o, $(call stripsrc,$(HISTFACTORYS:.cxx=.o)))


HISTFACTORYDEP  := $(HISTFACTORYO:.o=.d) $(HISTFACTORYDO:.o=.d)

HISTFACTORYLIB  := $(LPATH)/libHistFactory.$(SOEXT)
HISTFACTORYMAP  := $(HISTFACTORYLIB:.$(SOEXT)=.rootmap)

# used in the main Makefile
ALLHDRS      += $(patsubst $(MODDIRI)/RooStats/HistFactory/%.h,include/RooStats/HistFactory/%.h,$(HISTFACTORYH))
ALLLIBS      += $(HISTFACTORYLIB)
ALLMAPS      += $(HISTFACTORYMAP)

# include all dependency files
INCLUDEFILES += $(HISTFACTORYDEP)

#needed since include are in inc and not inc/RooStats
#HISTFACTORYH_DIC   := $(subst $(MODDIRI),include/RooStats/HistFactory,$(HISTFACTORYH))
##### local rules #####
.PHONY:         all-$(MODNAME) clean-$(MODNAME) distclean-$(MODNAME)

include/RooStats/HistFactory/%.h:    $(HISTFACTORYDIRI)/RooStats/HistFactory/%.h
		@(if [ ! -d "include/RooStats/HistFactory" ]; then    \
		   mkdir -p include/RooStats/HistFactory;             \
		fi)
		cp $< $@

$(HISTFACTORYLIB): $(HISTFACTORYO) $(HISTFACTORYDO) $(ORDER_) $(MAINLIBS) \
                $(HISTFACTORYLIBDEP) 
		echo "histfactorylib"
		@$(MAKELIB) $(PLATFORM) $(LD) "$(LDFLAGS)" \
		   "$(SOFLAGS)" libHistFactory.$(SOEXT) $@ \
		   "$(HISTFACTORYO) $(HISTFACTORYDO)" \
		   "$(HISTFACTORYLIBEXTRA)"

$(HISTFACTORYDS):  $(HISTFACTORYH) $(HISTFACTORYL) $(ROOTCINTTMPDEP)
		$(MAKEDIR)
		@echo "Generating dictionary $@..."
		$(ROOTCINTTMP) -f $@ -c $(HISTFACTORYH) $(HISTFACTORYL)

$(HISTFACTORYMAP): $(RLIBMAP) $(MAKEFILEDEP) $(HISTFACTORYL)
		$(RLIBMAP) -o  $@ -l $(HISTFACTORYLIB) \
		   -d $(HISTFACTORYLIBDEPM) -c $(HISTFACTORYL)

$(HF_MAKEWORKSPACEEXE): $(HF_MAKEWORKSPACEEXEO) $(ROOTLIBSDEP) $(RINTLIB) \
                $(HISTFACTORYLIBDEPM) $(HF_PREPAREHISTFACTORY) $(HISTFACTORYLIB)
		echo "amkeworkspacexe"
		$(LD) $(LDFLAGS) -o $@ $(HF_MAKEWORKSPACEEXEO) $(ROOTICON) \
		   $(ROOTULIBS) $(RPATH) $(ROOTLIBS) $(RINTLIBS) $(HF_LIBS) \
		   $(SYSLIBS)

$(HF_PREPAREHISTFACTORY): $(MODDIRC)/prepareHistFactory
		cp $(MODDIRC)/prepareHistFactory $@
		chmod +x $@

ALLEXECS     += $(HF_MAKEWORKSPACEEXE)

all-$(MODNAME): $(HISTFACTORYLIB) $(HISTFACTORYMAP) $(HF_MAKEWORKSPACEEXE)

clean-$(MODNAME):
		@rm -f $(HISTFACTORYO) $(HISTFACTORYDO)

clean::         clean-$(MODNAME)

distclean-$(MODNAME): clean-$(MODNAME)
		@rm -rf $(HISTFACTORYDEP) $(HISTFACTORYLIB) $(HISTFACTORYMAP) \
		   $(HISTFACTORYDS) $(HISTFACTORYDH)

distclean::     distclean-$(MODNAME)
