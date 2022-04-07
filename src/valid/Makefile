###############################################################################
#  2006/09/23 - Juan Pedro Mellado
#               Cleaned
#  1996/11/22 - Carlos Pantano
#               Created
###############################################################################
PRJS=\
     partial    \
     efilter    \
     poisson    \
     pentad     \
     heptad     \
     nonad      \
     matrix     \
     integral   \
     helmholtz_fde  \
     helmholtz_fxz  \
     tgvortex   \
     diffusion
#     superlayer \

PASSVARS = FF="$(FF)"               \
           LD="$(LD)"               \
	   CPPFLAGS="$(CPPFLAGS)"   \
           FFLAGS="$(FFLAGS)"       \
           EXTLIB="$(EXTLIB)"       \
           CPP="$(CPP)"             \
           LPRE="$(LPRE)"           \
	   DCPPFLAGS="$(DCPPFLAGS)" \
           DFFLAGS="$(DFFLAGS)"     \
	   DLDFLAGS="$(DLDFLAGS)"   \
           RCPPFLAGS="$(RCPPFLAGS)" \
           RFFLAGS="$(RFFLAGS)"     \
	   RLDFLAGS="$(RLDFLAGS)"

##########################################################################
#   Directives
##########################################################################
installrl:
	@(for prj in $(PRJS); do (cd $$prj; make installrl $(PASSVARS) ) done )

installdb: 
	@(for prj in $(PRJS); do (cd $$prj; make installdb $(PASSVARS) ) done )

##########################################################################
clean: 
	@(for prj in $(PRJS); do (cd $$prj; make clean SRC="$(SRC)" ) done )

proper:
	-/bin/rm -rf *~
	@(for prj in $(PRJS); do (cd $$prj; make proper SRC="$(SRC)" ) done )

cleandb: 
	@(for prj in $(PRJS); do (cd $$prj; make cleandb SRC="$(SRC)" ) done )
