# Shared build configuration for climber-x (dependency wiring).
#
# Loaded *after* the compiler and machine fragments (configme assembles them in
# the order: compiler -> machine -> netCDF -> common). It references variables
# those provide:
#   FFLAGS_BASE / FFLAGS_OPENMP / CPPFLAGS_PP  (compiler fragment)
#   INC_NC / LIB_NC                            (machine or auto-detected netCDF)
#
# External dependency repos (coordinates, fesm-utils, yelmo) live at the
# climber-x root — the layout configme uses for every orchestrator. VILMA is
# optional and user-provided under src/vilma (not managed by configme); it is
# only referenced by the fully-coupled (FULL) build variant.
#
# NOTE: FFLAGS_CLIM / FFLAGS_FULL are composed from $(FFLAGS_BASE), NOT $(FFLAGS):
# the template reassigns `FFLAGS = $(FFLAGS_CLIM)` per build target, so referring
# to $(FFLAGS) here would be a recursive self-reference.

# --- coordinates
COORDROOT = coordinates
INC_COORD = -I${COORDROOT}/libcoordinates/include
LIB_COORD = -L${COORDROOT}/libcoordinates/include -lcoordinates

# --- fesm-utils (serial build by default; OpenMP variants swapped in below)
FESMUTILSROOT = fesm-utils/utils
INC_FESMUTILS = -I${FESMUTILSROOT}/include-serial
LIB_FESMUTILS = -L${FESMUTILSROOT}/include-serial -lfesmutils

FFTWROOT = fesm-utils/fftw-serial
INC_FFTW = -I${FFTWROOT}/include
LIB_FFTW = -L${FFTWROOT}/lib -lfftw3 -lm

LISROOT = fesm-utils/lis-serial
INC_LIS = -I${LISROOT}/include
LIB_LIS = -L${LISROOT}/lib/ -llis

# --- yelmo (ice sheet; used by the FULL build variant)
YELMOROOT = yelmo
INC_YELMO = -I${YELMOROOT}/libyelmo/include
LIB_YELMO = -L${YELMOROOT}/libyelmo/include -lyelmo

# --- VILMA (optional solid-earth lib; user-provided, not managed by configme)
VILMAROOT = src/vilma
INC_VILMA = -I${VILMAROOT}/include
LIB_VILMA = ${VILMAROOT}/lib/vega_pism.a

# OpenMP build (make openmp=1, the climber-x default): swap the serial
# dependency builds for their OpenMP variants. The compiler's $(FFLAGS_OPENMP)
# is appended to FFLAGS_CLIM / FFLAGS_FULL by the template itself.
ifeq ($(openmp), 1)
    INC_FESMUTILS = -I${FESMUTILSROOT}/include-omp
    LIB_FESMUTILS = -L${FESMUTILSROOT}/include-omp -lfesmutils

    FFTWROOT = fesm-utils/fftw-omp
    INC_FFTW = -I${FFTWROOT}/include
    LIB_FFTW = -L${FFTWROOT}/lib -lfftw3_omp -lfftw3 -lm

    LISROOT = fesm-utils/lis-omp
    INC_LIS = -I${LISROOT}/include
    LIB_LIS = -L${LISROOT}/lib/ -llis
endif

# --- compile-flag sets: climate-only (CLIM) vs fully-coupled (FULL, adds ice +
# VILMA). The template appends -DVERSION and, for openmp=1, $(FFLAGS_OPENMP).
# $(CPPFLAGS_PP) is the compiler's Fortran-preprocessor flag (-fpp Intel, -cpp GNU).
CPPFLAGS_CLIM = $(CPPFLAGS_PP)
CPPFLAGS_FULL = $(CPPFLAGS_PP) -DVILMA

FFLAGS_CLIM = $(FFLAGS_BASE) $(INC_NC) $(INC_FESMUTILS) $(INC_COORD) $(INC_FFTW)
FFLAGS_FULL = $(FFLAGS_BASE) $(INC_NC) $(INC_FESMUTILS) $(INC_LIS) $(INC_COORD) $(INC_YELMO) $(INC_VILMA) $(INC_FFTW)

# Extra link flags. -Wl,-zmuldefs works around duplicate symbols in the static
# deps (the default on Linux). A machine fragment disables it with
# `LFLAGS_EXTRA =` (macOS ld rejects -zmuldefs).
LFLAGS_EXTRA ?= -Wl,-zmuldefs
LFLAGS_CLIM = $(LIB_NC) $(LIB_COORD) $(LIB_FESMUTILS) $(LIB_FFTW) $(LFLAGS_EXTRA)
LFLAGS_FULL = $(LIB_NC) $(LIB_COORD) $(LIB_FFTW) $(LIB_LIS) $(LIB_YELMO) $(LIB_VILMA) $(LIB_FESMUTILS) $(LFLAGS_EXTRA)
