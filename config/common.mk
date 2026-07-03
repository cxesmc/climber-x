# Shared build configuration for climber-x (dependency wiring).
#
# Loaded *after* the compiler and machine fragments (configme assembles them in
# the order: compiler -> machine -> netCDF -> common). It references variables
# those provide:
#   FFLAGS_BASE / FFLAGS_OPENMP / CPPFLAGS_PP  (compiler fragment)
#   INC_NC / LIB_NC                            (machine or auto-detected netCDF)
#
# External dependency repos (fesm-utils, yelmo) live at the
# climber-x root — the layout configme uses for every orchestrator. VILMA is
# optional and user-provided under src/vilma (not managed by configme); it is
# only referenced by the fully-coupled (FULL) build variant.
#
# NOTE: FFLAGS_CLIM / FFLAGS_FULL are composed from $(FFLAGS_BASE), NOT $(FFLAGS):
# the template reassigns `FFLAGS = $(FFLAGS_CLIM)` per build target, so referring
# to $(FFLAGS) here would be a recursive self-reference.

# --- fesm-utils (serial build by default; OpenMP variants swapped in below)
# FESMUTILS_LIB is the on-disk libfesmutils.a (the artifact behind LIB_FESMUTILS).
# Climber objects that `use` a fesm-utils module (coords/ncio/nml) depend on it
# (see the blanket rule at the end of Makefile_climber.mk) so a rebuilt fesm-utils
# with a changed module interface forces a recompile against the new .mod, rather
# than silently relinking a stale object (ABI mismatch -> segfault).
FESMUTILSROOT = fesm-utils/utils
INC_FESMUTILS = -I${FESMUTILSROOT}/include-serial
LIB_FESMUTILS = -L${FESMUTILSROOT}/include-serial -lfesmutils
FESMUTILS_LIB = ${FESMUTILSROOT}/include-serial/libfesmutils.a

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

# --- FastHydrology (subglacial hydrology; required transitively by yelmo >= v2.2:
# yelmo_defs declares `type(hydro_class) :: hyd`, so every yelmo module references
# the fasthydro .mod files and libyelmo.a references libfasthydro symbols).
# Built into a single include/ dir (no serial/omp split). FastHydrology itself
# pulls in FFTW, which is already wired above.
FASTHYDROROOT = yelmo/FastHydrology
INC_FASTHYDRO = -I${FASTHYDROROOT}/include
LIB_FASTHYDRO = -L${FASTHYDROROOT}/include -lfasthydro

# --- VILMA (optional solid-earth lib; user-provided, not managed by configme)
VILMAROOT = src/vilma
INC_VILMA = -I${VILMAROOT}/include
LIB_VILMA = ${VILMAROOT}/lib/vega_pism.a

# --- FastEarth3D (alternative solid-earth model; cloned at the climber-x root by
# configme, like yelmo, and built from source against the SAME fesm-utils as
# climber-x: `make fastearth-static` -> obj/libfastearth.a, .mod files in obj/).
# It is a swappable alternative to VILMA, selected at runtime with i_geo=3.
FASTEARTHROOT = FastEarth3D
INC_FASTEARTH = -I${FASTEARTHROOT}/obj
LIB_FASTEARTH = ${FASTEARTHROOT}/obj/libfastearth.a

# --- SHTns (spherical-harmonic transforms; provided by fesm-utils, pulled in by
# FastEarth3D). FFTW (wired above) is linked after it (SHTns calls FFTW). Serial
# by default; the OpenMP variant (shtns-omp, -lshtns_omp) is swapped in below for
# openmp=1 to match FastEarth3D's own openmp= dependency swap.
SHTNSROOT = fesm-utils/shtns-serial
INC_SHTNS = -I${SHTNSROOT}/include
LIB_SHTNS = -L${SHTNSROOT}/lib -lshtns

# OpenMP build (make openmp=1, the climber-x default): swap the serial
# dependency builds for their OpenMP variants. The compiler's $(FFLAGS_OPENMP)
# is appended to FFLAGS_CLIM / FFLAGS_FULL by the template itself.
ifeq ($(openmp), 1)
    INC_FESMUTILS = -I${FESMUTILSROOT}/include-omp
    LIB_FESMUTILS = -L${FESMUTILSROOT}/include-omp -lfesmutils
    FESMUTILS_LIB = ${FESMUTILSROOT}/include-omp/libfesmutils.a

    FFTWROOT = fesm-utils/fftw-omp
    INC_FFTW = -I${FFTWROOT}/include
    LIB_FFTW = -L${FFTWROOT}/lib -lfftw3_omp -lfftw3 -lm

    LISROOT = fesm-utils/lis-omp
    INC_LIS = -I${LISROOT}/include
    LIB_LIS = -L${LISROOT}/lib/ -llis

    SHTNSROOT = fesm-utils/shtns-omp
    INC_SHTNS = -I${SHTNSROOT}/include
    LIB_SHTNS = -L${SHTNSROOT}/lib -lshtns_omp
endif

# --- solid-earth backend selection (FULL build). The vilma= / fastearth= toggles
# (defined in config/Makefile, default 1) gate each backend's -D flag plus its
# include/link flags, accumulated into the *_EARTH pieces spliced into the FULL
# sets below. -DVILMA and -DFASTEARTH are independent; both are on by default and
# chosen at runtime by i_geo (2=VILMA, 3=FastEarth3D). A backend toggled off is
# compiled as a stub that aborts cleanly if selected (see vilma.F90/fastearth.F90).
CPPFLAGS_EARTH =
INC_EARTH =
LIB_EARTH =
ifeq ($(vilma),1)
    CPPFLAGS_EARTH += -DVILMA
    INC_EARTH      += $(INC_VILMA)
    LIB_EARTH      += $(LIB_VILMA)
endif
ifeq ($(fastearth),1)
    CPPFLAGS_EARTH += -DFASTEARTH
    INC_EARTH      += $(INC_FASTEARTH) $(INC_SHTNS)
    # order: LIB_FASTEARTH precedes LIB_SHTNS (it calls SHTns), which precedes the
    # trailing LIB_FFTW in LFLAGS_FULL (SHTns calls FFTW).
    LIB_EARTH      += $(LIB_FASTEARTH) $(LIB_SHTNS)
endif

# --- compile-flag sets: climate-only (CLIM) vs fully-coupled (FULL, adds ice +
# the selected solid-earth backends). The template appends -DVERSION and, for
# openmp=1, $(FFLAGS_OPENMP). $(CPPFLAGS_PP) is the compiler's Fortran-preprocessor
# flag (-fpp Intel, -cpp GNU).
CPPFLAGS_CLIM = $(CPPFLAGS_PP)
CPPFLAGS_FULL = $(CPPFLAGS_PP) $(CPPFLAGS_EARTH)

FFLAGS_CLIM = $(FFLAGS_BASE) $(INC_NC) $(INC_FESMUTILS) $(INC_FFTW)
FFLAGS_FULL = $(FFLAGS_BASE) $(INC_NC) $(INC_FESMUTILS) $(INC_LIS) $(INC_YELMO) $(INC_FASTHYDRO) $(INC_EARTH) $(INC_FFTW)

# Extra link flags. -Wl,-zmuldefs works around duplicate symbols in the static
# deps (the default on Linux). A machine fragment disables it with
# `LFLAGS_EXTRA =` (macOS ld rejects -zmuldefs).
LFLAGS_EXTRA ?= -Wl,-zmuldefs
LFLAGS_CLIM = $(LIB_NC) $(LIB_FESMUTILS) $(LIB_FFTW) $(LFLAGS_EXTRA)
# LIB_FASTHYDRO follows LIB_YELMO (yelmo references its symbols); the trailing
# LIB_FFTW resolves the FFTW symbols pulled in by fasthydro and fesmutils
# (static archives resolve left-to-right, so deps must come after dependents).
# LIB_EARTH (the selected solid-earth libs, see above) precedes LIB_FESMUTILS,
# which resolves the coords/ncio symbols VILMA and FastEarth3D share.
LFLAGS_FULL = $(LIB_NC) $(LIB_FFTW) $(LIB_LIS) $(LIB_YELMO) $(LIB_FASTHYDRO) $(LIB_EARTH) $(LIB_FESMUTILS) $(LIB_FFTW) $(LFLAGS_EXTRA)
