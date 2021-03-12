###############################################################################
#  2006/09/23 - Juan Pedro Mellado
#               Cleaned
#  1996/11/22 - Carlos Pantano
#               Created
###############################################################################
BINS=vhelmholtz_fxz.x

SRCS=\
     vhelmholtz_fxz.f 	\

MODS=\
     dns_global.o \
     thermo_global.o \
     dns_mpi.o

INCS=\
     types.h	     \
     dns_bcs_inf.h   \
     dns_bcs_out.h   \
     dns_const.h     \
     avgij_map.h     \
     les_const.h     \
     les_avgij_map.h

LIBS=\
     libfields.a  \
     libflow.a    \
     libdns.a     \
     libio.a      \
     libparallel.a\
     libthermo.a  \
     libpdfs.a    \
     libavg.a     \
     libfdm.a     \
     libfilter.a  \
     liblinear.a  \
     libquad.a    \
     libsplines.a \
     libutils.a

INC_PATH=../../include
LIB_PATH=../../lib
BIN_PATH=../../bin

RMODS=$(MODS:%=$(RLIB_PATH)/%)
DMODS=$(MODS:%=$(DLIB_PATH)/%)

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

#########################################################################
#   Directives
#########################################################################
installrl: release $(BIN_PATH)/release/$(BINS)
installdb: debug $(BIN_PATH)/debug/$(BINS)

release debug:
	mkdir $@

$(BIN_PATH)/release/$(BINS): release/$(BINS)
	@echo "Installing : $(BINS) release"

$(BIN_PATH)/debug/$(BINS): debug/$(BINS)
	@echo "Installing : $(BINS) debug"

release/$(BINS): $(ROBJS) $(RLIBS) $(RMODS)
	$(LD) release/$(BINS) $(RLDFLAGS) $(ROBJS) $(RLINK) $(RMODS)

debug/$(BINS): $(DOBJS) $(DLIBS) $(DMODS)
	$(LD) debug/$(BINS) $(DLDFLAGS) $(DOBJS) $(DLINK) $(DMODS)

#########################################################################
include ../../common.mak
