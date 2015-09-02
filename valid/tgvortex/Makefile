###############################################################################
#  2006/09/23 - Juan Pedro Mellado
#               Cleaned
#  1996/11/22 - Carlos Pantano
#               Created
###############################################################################
BINS=vtgvortex.x

SRCS=\
     vtgvortex.f

MODS=\
     dns_global.o    \
     thermo_global.o

INCS=\
     types.h

LIBS=\
     libfields.a  \
     libdns.a     \
     libthermo.a  \
     libflow.a    \
     libio.a      \
     libfdm.a     \
     liblinear.a  \
     libavg.a     \
     libutils.a 

INC_PATH=../../include
LIB_PATH=../../lib
BIN_PATH=../../bin

#########################################################################
# No modification beyond this points should be required
#########################################################################
INCLUDES=$(INCS:%=$(INC_PATH)/%)

ROBJS=$(SRCS:%.f=release/%.o)
DOBJS=$(SRCS:%.f=debug/%.o)

RLIB_PATH=$(LIB_PATH:%=%/release)
RLIBS=$(LIBS:%=$(RLIB_PATH)/%)
RLINK=$(LPRE)-L$(RLIB_PATH) $(LIBS:lib%.a=-l%) $(EXTLIB)

DLIB_PATH=$(LIB_PATH:%=%/debug)
DLIBS=$(LIBS:%=$(DLIB_PATH)/%)
DLINK=$(LPRE)-L$(DLIB_PATH) $(LIBS:lib%.a=-l%) $(EXTLIB)

RMODS=$(MODS:%=$(RLIB_PATH)/%)
DMODS=$(MODS:%=$(DLIB_PATH)/%)

#########################################################################
#   Directives
#########################################################################
installrl: release \
           $(BIN_PATH)/release/$(BINS) \
           $(BIN_PATH)/release/$(BINS2)

installdb: debug \
           $(BIN_PATH)/debug/$(BINS) \
           $(BIN_PATH)/debug/$(BINS2)

release debug:
	mkdir $@

$(BIN_PATH)/release/%.x: release/%.x
	@echo "Installing : $*.x release"

$(BIN_PATH)/debug/%.x: debug/%.x
	@echo "Installing : $*.x debug"

#########################################################################
release/$(BINS): $(ROBJS) $(RLIBS) $(RMODS)
	$(LD) release/$(BINS) $(RLDFLAGS) $(ROBJS) $(RLINK) $(RMODS)
debug/$(BINS): $(DOBJS) $(DLIBS) $(DMODS)
	$(LD) debug/$(BINS) $(DLDFLAGS) $(DOBJS) $(DLINK) $(DMODS)

#########################################################################
include ../../common.mak
