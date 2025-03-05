.SUFFIXES: .f .F .F90 .f90 .o .mod
.SHELL: /bin/sh

# PATH options
objdir = obj
srcdir = src

# Command-line options at make call
debug  ?= 0  
openmp ?= 1

defaulttarget: usage

## COMPILER CONFIGURATION ##
# (should be loaded from config directory)

<COMPILER_CONFIGURATION>

# Determine whether to use normal flags (default), 
# debugging flags (debug=1) or profiling flags (debug=2)
DFLAGS   = $(DFLAGS_NODEBUG)
ifeq ($(debug), 1)
	DFLAGS   = $(DFLAGS_DEBUG)
endif
ifeq ($(debug), 2)
	DFLAGS   = $(DFLAGS_PROFILE)
endif

# Determine whether to use openmp flags 
ifeq ($(openmp), 1)
	FFLAGS_FULL := $(FFLAGS_FULL) $(FFLAGS_OPENMP)
	FFLAGS_CLIM := $(FFLAGS_CLIM) $(FFLAGS_OPENMP)
endif

###############################################
##							
## List of rules and source files
##
###############################################

include config/Makefile_climber.mk

###############################################
##
## Compilation of complete programs
##
###############################################

# Compile the static libraries using their own Makefiles
coord-static: 
	$(MAKE) -C $(COORDROOT) coord-static openmp=$(openmp)

yelmo-static: 
	$(MAKE) -C $(YELMOROOT) yelmo-static openmp=$(openmp)
#	cp $(YELMOROOT)/git_yelmo.txt ./

climber-clim: FFLAGS = $(FFLAGS_CLIM)
climber-clim: LFLAGS = $(LFLAGS_CLIM)
climber-clim: CPPFLAGS = $(CPPFLAGS_CLIM)
climber-clim: LDFLAGS = $(FFLAGS) $(DFLAGS) $(CPPFLAGS)
climber-clim: coord-static \
				utils bnd geo-clim \
				atm ocn sic lnd \
				bgc-dummy co2 ch4 \
				ice-dummy smb-dummy bmb-dummy main-clim
	$(FC) $(LDFLAGS) -o climber.x $(climber_clim_obj) $(LFLAGS)
	@echo " "
	@echo "    climber.x is ready."
	@echo "    [clim]"
	@echo " "

climber-clim-bgc: FFLAGS = $(FFLAGS_CLIM)
climber-clim-bgc: LFLAGS = $(LFLAGS_CLIM)
climber-clim-bgc: CPPFLAGS = $(CPPFLAGS_CLIM)
climber-clim-bgc: LDFLAGS = $(FFLAGS) $(DFLAGS) $(CPPFLAGS)
climber-clim-bgc: coord-static \
				utils bnd geo-clim \
				atm ocn sic lnd \
				bgc co2 ch4 \
				ice-dummy smb-dummy bmb-dummy main-clim-bgc
	$(FC) $(LDFLAGS) -o climber.x $(climber_clim_bgc_obj) $(LFLAGS)
	@echo " "
	@echo "    climber.x is ready."
	@echo "    [clim,bgc]"
	@echo " "

climber-clim-ice: FFLAGS = $(FFLAGS_FULL)
climber-clim-ice: LFLAGS = $(LFLAGS_FULL)
climber-clim-ice: CPPFLAGS = $(CPPFLAGS_FULL)
climber-clim-ice: LDFLAGS = $(FFLAGS) $(DFLAGS) $(CPPFLAGS)
climber-clim-ice: coord-static yelmo-static \
				utils bnd geo \
				atm ocn sic lnd \
				bgc-dummy co2 ch4 \
				ice smb bmb main-clim-ice
	$(FC) $(LDFLAGS) -o climber.x $(climber_clim_ice_obj) $(LFLAGS)
	@echo " "
	@echo "    climber.x is ready."
	@echo "    [clim,ice]"
	@echo " "

climber-clim-bgc-ice: FFLAGS = $(FFLAGS_FULL)
climber-clim-bgc-ice: LFLAGS = $(LFLAGS_FULL)
climber-clim-bgc-ice: CPPFLAGS = $(CPPFLAGS_FULL)
climber-clim-bgc-ice: LDFLAGS = $(FFLAGS) $(DFLAGS) $(CPPFLAGS)
climber-clim-bgc-ice: coord-static yelmo-static \
				utils bnd geo \
				atm ocn sic lnd \
				bgc co2 ch4 \
				ice smb bmb main
	$(FC) $(LDFLAGS) -o climber.x $(climber_clim_bgc_ice_obj) $(LFLAGS)
	@echo " "
	@echo "    climber.x is ready."
	@echo "    [clim,bgc,ice]"
	@echo " "

climber: climber-clim-bgc-ice

check-args: 
	@echo $(openmp) $(debug)
	@echo $(FFLAGS)
	@echo $(CPPFLAGS)
	@echo $(DFLAGS)
	@echo $(LDFLAGS)

###############################################
##
## Compilation of individual modules
##
###############################################

atm: FFLAGS = $(FFLAGS_CLIM)
atm: LFLAGS = $(LFLAGS_CLIM)
atm: CPPFLAGS = $(CPPFLAGS_CLIM)
atm: LDFLAGS = $(FFLAGS) $(DFLAGS) $(CPPFLAGS)
atm: $(obj_atm)
	@echo ""
	@echo "atm files compiled successfully."
	@echo ""

atm-out: FFLAGS = $(FFLAGS_CLIM)
atm-out: LFLAGS = $(LFLAGS_CLIM)
atm-out: CPPFLAGS = $(CPPFLAGS_CLIM)
atm-out: LDFLAGS = $(FFLAGS) $(DFLAGS) $(CPPFLAGS)
atm-out: $(objdir)/atm_out.o
	@echo ""
	@echo "atm-out files compiled successfully."
	@echo ""

bgc: FFLAGS = $(FFLAGS_CLIM)
bgc: LFLAGS = $(LFLAGS_CLIM)
bgc: CPPFLAGS = $(CPPFLAGS_CLIM)
bgc: LDFLAGS = $(FFLAGS) $(DFLAGS) $(CPPFLAGS)
bgc: $(obj_bgc)
	@echo ""
	@echo "bgc files compiled successfully."
	@echo ""

bgc-dummy: FFLAGS = $(FFLAGS_CLIM)
bgc-dummy: LFLAGS = $(LFLAGS_CLIM)
bgc-dummy: CPPFLAGS = $(CPPFLAGS_CLIM)
bgc-dummy: LDFLAGS = $(FFLAGS) $(DFLAGS) $(CPPFLAGS)
bgc-dummy: $(obj_bgc_dummy)
	@echo ""
	@echo "bgc-dummy files compiled successfully."
	@echo ""

bnd: FFLAGS = $(FFLAGS_CLIM)
bnd: LFLAGS = $(LFLAGS_CLIM)
bnd: CPPFLAGS = $(CPPFLAGS_CLIM)
bnd: LDFLAGS = $(FFLAGS) $(DFLAGS) $(CPPFLAGS)
bnd: $(obj_bnd)
	@echo ""
	@echo "bnd files compiled successfully."
	@echo ""

ch4: FFLAGS = $(FFLAGS_CLIM)
ch4: LFLAGS = $(LFLAGS_CLIM)
ch4: CPPFLAGS = $(CPPFLAGS_CLIM)
ch4: LDFLAGS = $(FFLAGS) $(DFLAGS) $(CPPFLAGS)
ch4: $(obj_ch4)
	@echo ""
	@echo "ch4 files compiled successfully."
	@echo ""

co2: FFLAGS = $(FFLAGS_CLIM)
co2: LFLAGS = $(LFLAGS_CLIM)
co2: CPPFLAGS = $(CPPFLAGS_CLIM)
co2: LDFLAGS = $(FFLAGS) $(DFLAGS) $(CPPFLAGS)
co2: $(obj_co2)
	@echo ""
	@echo "co2 files compiled successfully."
	@echo ""

geo: FFLAGS = $(FFLAGS_FULL)
geo: LFLAGS = $(LFLAGS_FULL)
geo: CPPFLAGS = $(CPPFLAGS_FULL)
geo: LDFLAGS = $(FFLAGS) $(DFLAGS) $(CPPFLAGS)
geo: $(obj_geo)
	@echo ""
	@echo "geo files compiled successfully."
	@echo ""

geo-clim: FFLAGS = $(FFLAGS_CLIM)
geo-clim: LFLAGS = $(LFLAGS_CLIM)
geo-clim: CPPFLAGS = $(CPPFLAGS_CLIM)
geo-clim: LDFLAGS = $(FFLAGS) $(DFLAGS) $(CPPFLAGS)
geo-clim: $(obj_geo)
	@echo ""
	@echo "geo files compiled successfully."
	@echo ""

ice: FFLAGS = $(FFLAGS_FULL)
ice: LFLAGS = $(LFLAGS_FULL)
ice: CPPFLAGS = $(CPPFLAGS_FULL)
ice: LDFLAGS = $(FFLAGS) $(DFLAGS) $(CPPFLAGS)
ice: ice-sico yelmo-static $(obj_ice)
	@echo ""
	@echo "ice files compiled successfully."
	@echo ""

ice-dummy: FFLAGS = $(FFLAGS_CLIM)
ice-dummy: LFLAGS = $(LFLAGS_CLIM)
ice-dummy: CPPFLAGS = $(CPPFLAGS_CLIM)
ice-dummy: LDFLAGS = $(FFLAGS) $(DFLAGS) $(CPPFLAGS)
ice-dummy: $(obj_ice_dummy)
	@echo ""
	@echo "ice-dummy files compiled successfully."
	@echo ""

ice-sico: FFLAGS = $(FFLAGS_FULL)
ice-sico: LFLAGS = $(LFLAGS_FULL)
ice-sico: CPPFLAGS = $(CPPFLAGS_FULL)
ice-sico: LDFLAGS = $(FFLAGS) $(DFLAGS) $(CPPFLAGS)
ice-sico: $(obj_ice_sico)
	@echo ""
	@echo "ice-sico files compiled successfully."
	@echo ""

bmb: FFLAGS = $(FFLAGS_CLIM)
bmb: LFLAGS = $(LFLAGS_CLIM)
bmb: CPPFLAGS = $(CPPFLAGS_CLIM)
bmb: LDFLAGS = $(FFLAGS) $(DFLAGS) $(CPPFLAGS)
bmb: $(obj_bmb)
	@echo ""
	@echo "bmb files compiled successfully."
	@echo ""

bmb-dummy: FFLAGS = $(FFLAGS_CLIM)
bmb-dummy: LFLAGS = $(LFLAGS_CLIM)
bmb-dummy: CPPFLAGS = $(CPPFLAGS_CLIM)
bmb-dummy: LDFLAGS = $(FFLAGS) $(DFLAGS) $(CPPFLAGS)
bmb-dummy: $(obj_bmb_dummy)
	@echo ""
	@echo "bmb-dummy files compiled successfully."
	@echo ""

lnd: FFLAGS = $(FFLAGS_CLIM)
lnd: LFLAGS = $(LFLAGS_CLIM)
lnd: CPPFLAGS = $(CPPFLAGS_CLIM)
lnd: LDFLAGS = $(FFLAGS) $(DFLAGS) $(CPPFLAGS)
lnd: $(obj_lnd)
	@echo ""
	@echo "lnd files compiled successfully."
	@echo ""

lndvc: FFLAGS = $(FFLAGS_CLIM)
lndvc: LFLAGS = $(LFLAGS_CLIM)
lndvc: CPPFLAGS = $(CPPFLAGS_CLIM)
lndvc: LDFLAGS = $(FFLAGS) $(DFLAGS) $(CPPFLAGS)
lndvc: $(objdir)/lndvc_model.o
	@echo ""
	@echo "lndvc_model module and dependencies compiled successfully."
	@echo ""

main: FFLAGS = $(FFLAGS_FULL)
main: LFLAGS = $(LFLAGS_FULL)
main: CPPFLAGS = $(CPPFLAGS_FULL)
main: LDFLAGS = $(FFLAGS) $(DFLAGS) $(CPPFLAGS)
main: $(obj_main)
	@echo ""
	@echo "main files compiled successfully."
	@echo ""

main-clim: FFLAGS = $(FFLAGS_CLIM)
main-clim: LFLAGS = $(LFLAGS_CLIM)
main-clim: CPPFLAGS = $(CPPFLAGS_CLIM)
main-clim: LDFLAGS = $(FFLAGS) $(DFLAGS) $(CPPFLAGS)
main-clim: $(obj_main_clim)
	@echo ""
	@echo "main files compiled successfully."
	@echo ""

main-clim-bgc: FFLAGS = $(FFLAGS_CLIM)
main-clim-bgc: LFLAGS = $(LFLAGS_CLIM)
main-clim-bgc: CPPFLAGS = $(CPPFLAGS_CLIM)
main-clim-bgc: LDFLAGS = $(FFLAGS) $(DFLAGS) $(CPPFLAGS)
main-clim-bgc: $(obj_main_clim_bgc)
	@echo ""
	@echo "main files compiled successfully."
	@echo ""

main-clim-ice: FFLAGS = $(FFLAGS_FULL)
main-clim-ice: LFLAGS = $(LFLAGS_FULL)
main-clim-ice: CPPFLAGS = $(CPPFLAGS_FULL)
main-clim-ice: LDFLAGS = $(FFLAGS) $(DFLAGS) $(CPPFLAGS)
main-clim-ice: $(obj_main_clim_ice)
	@echo ""
	@echo "main files compiled successfully."
	@echo ""

ocn: FFLAGS = $(FFLAGS_CLIM)
ocn: LFLAGS = $(LFLAGS_CLIM)
ocn: CPPFLAGS = $(CPPFLAGS_CLIM)
ocn: LDFLAGS = $(FFLAGS) $(DFLAGS) $(CPPFLAGS)
ocn: $(obj_ocn)
	@echo ""
	@echo "ocn files compiled successfully."
	@echo ""

sic: FFLAGS = $(FFLAGS_CLIM)
sic: LFLAGS = $(LFLAGS_CLIM)
sic: CPPFLAGS = $(CPPFLAGS_CLIM)
sic: LDFLAGS = $(FFLAGS) $(DFLAGS) $(CPPFLAGS)
sic: $(obj_sic)
	@echo ""
	@echo "sic files compiled successfully."
	@echo ""

smb: FFLAGS = $(FFLAGS_CLIM)
smb: LFLAGS = $(LFLAGS_CLIM)
smb: CPPFLAGS = $(CPPFLAGS_CLIM)
smb: LDFLAGS = $(FFLAGS) $(DFLAGS) $(CPPFLAGS)
smb: $(obj_smb)
	@echo ""
	@echo "smb files compiled successfully."
	@echo ""

smb-dummy: FFLAGS = $(FFLAGS_CLIM)
smb-dummy: LFLAGS = $(LFLAGS_CLIM)
smb-dummy: CPPFLAGS = $(CPPFLAGS_CLIM)
smb-dummy: LDFLAGS = $(FFLAGS) $(DFLAGS) $(CPPFLAGS)
smb-dummy: $(obj_smb_dummy)
	@echo ""
	@echo "smb-dummy files compiled successfully."
	@echo ""

utils: FFLAGS = $(FFLAGS_CLIM)
utils: LFLAGS = $(LFLAGS_CLIM)
utils: CPPFLAGS = $(CPPFLAGS_CLIM)
utils: LDFLAGS = $(FFLAGS) $(DFLAGS) $(CPPFLAGS)
utils: $(obj_utils)
	@echo ""
	@echo "utils files compiled successfully."
	@echo ""

yelmo_restart: tests/yelmo_from_restart.f90 $(ncio_obj) $(nml_obj)
	$(FC) $(LDFLAGS) -o yelmo_from_restart.x $^ $(LFLAGS)

check: $(climber_obj)
	@echo $(climber_obj)
	@echo ""

docs: 
	ford ford.md

.PHONY : files
files:
	@echo $(climber_obj)

.PHONY : usage
usage:
	@echo ""
	@echo "    * USAGE * "
	@echo ""
	@echo " make climber    : compiles the main program climber.x"
	@echo " make clean      : clean object files"
	@echo " make module     : show required modules to load on PIK cluster"
	@echo " make showconfig : indicate current makefile configuration"
	@echo " make docs       : create automatic documentation with FORD"
	@echo ""

# target name config interferes with the same-named subdirectory
showconfig:
	@echo "----------------------"
	@echo "Makefile configuration"
	@echo "----------------------"
	@echo "compiler 		: $(FC)"
	@echo "host 			: $(shell hostname)"
	@echo "machine  		: $(machine)"
	@echo "required modules : $(modules)"
	@echo "debug 			: $(debug)"
module:
	@echo module load $(modules)

clean-atm:
	rm -f  $(obj_atm)

clean-geo:
	rm -f  $(obj_geo) $(obj_geo_dummy)

clean-lnd:
	rm -f  $(obj_lnd)

clean-lndvc:
	rm -f  $(obj_lndvc)

clean-bgc:
	rm -f  $(obj_bgc) $(obj_bgc_dummy)

clean-ice:
	rm -f  $(obj_ice) $(obj_ice_sico) $(obj_ice_dummy)

clean-bmb:
	rm -f  $(obj_bmb) $(obj_bmb_dummy)

clean-sic:
	rm -f  $(obj_sic)

clean-smb:
	rm -f  $(obj_smb) $(obj_smb_dummy)

clean-ocn:
	rm -f  $(obj_ocn)

clean:
	rm -f $(bindir)/*.x
	rm -f  *.x gmon.out $(objdir)/*.o $(objdir)/.*.o $(objdir)/*.mod $(objdir)/*.a $(objdir)/*.so
	rm -rf *.x.dSYM

cleanall: clean
	$(MAKE) -C $(COORDROOT) clean
	if [ -d "$(YELMOROOT)" ]; then $(MAKE) -C $(YELMOROOT) clean; fi
