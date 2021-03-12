###############################################################################
#  2006/09/23 - Juan Pedro Mellado
#               Cleaned
#  1996/11/22 - Carlos Pantano
#               Created
###############################################################################
BINS=vnonad.x

SRCS=\
     vnonad.f

INCS=\
     types.h

LIBS=\
     libdns.a   \
     libio.a    \
     libfdm.a  \
     liblinear.a\
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

release/$(BINS): $(ROBJS) $(RLIBS)
	$(LD) release/$(BINS) $(RLDFLAGS) $(ROBJS) $(RLINK)

debug/$(BINS): $(DOBJS) $(DLIBS)
	$(LD) debug/$(BINS) $(DLDFLAGS) $(DOBJS) $(DLINK)

#########################################################################
include ../../common.mak
