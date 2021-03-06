#
# Makefile to build the aed2 water quality library
# with hydrodynamic driver wrapper
#

objdir=obj
srcdir=src
libdir=lib
VERS=2.0.1ptm

OSTYPE=$(shell uname -s)

ifeq ($(AED2DIR),)
  AED2DIR=../libaed2
endif
OUTLIB=libtuflowfv_external_wq


INCLUDES+=-I${AED2DIR}/include
INCLUDES+=-I${AED2DIR}/mod
LIBAED2=aed2
LIBFVAED2=fvaed2
moddir=mod

ifeq ($(OSTYPE),Darwin)
  SHARED=-dynamiclib -undefined dynamic_lookup
  so_ext=dylib
else
  SHARED=-shared
  so_ext=so
endif

ifeq ($(F90),ifort)
  INCLUDES+=-I/opt/intel/include
  DEBUG_FFLAGS=-g -traceback
  OPT_FFLAGS=-O3 -openmp
  FFLAGS=-fpp -warn all -module ${moddir} -i-static -mp1 -warn nounused $(DEFINES)
  ifeq ($(WITH_CHECKS),true)
    FFLAGS+=-check
  endif
  FFLAGS+=-real-size 64
  LIBS+=-lifcore -lsvml
  LIBS+=-limf -lintlc
  LIBS+=-L/opt/intel/lib -Wl,-rpath=/opt/intel/lib
endif

FFLAGS+=-fPIC
LIBS += -L${AED2DIR} -l${LIBAED2}

LIBS+=-lnetcdff -lnetcdf

FFLAGS+=$(OPT_FFLAGS)

ifeq ($(DEBUG),true)
  FFLAGS+=$(DEBUG_FFLAGS)
endif

ifeq ($(PRECISION),1)
  TFFLAGS += -D_PRECISION=1
else ifeq ($(PRECISION),2)
  TFFLAGS += -D_PRECISION=2
else
  TFFLAGS += -D_PRECISION=1
endif

TFFLAGS += -g -DAED2 -DEXTERNAL_WQ=2
INCLUDES += -I${moddir}


FVOBJECTS=${objdir}/fv_zones.o ${objdir}/fv_aed2.o
OBJECTS=${objdir}/tuflowfv_external_wq.o

SOFLAGS = ${libdir}/lib${LIBFVAED2}.a ${AED2DIR}/lib/lib${LIBAED2}.a

ifeq ($(EXTERNAL_LIBS),shared)
  TARGET = ${libdir}/$(OUTLIB).${so_ext}
else
  TARGET = ${libdir}/$(OUTLIB).a
endif

all: ${TARGET}

${libdir}/lib${LIBFVAED2}.a: ${objdir} ${moddir} ${libdir} ${FVOBJECTS}
	ar -rv $@ ${FVOBJECTS} ${LDFLAGS}
	ranlib $@

${libdir}/${OUTLIB}.a: ${libdir}/lib${LIBFVAED2}.a ${OBJECTS}
	ar -rv $@ ${OBJECTS} ${LDFLAGS}
	ranlib $@

${libdir}/${OUTLIB}.${so_ext}: ${libdir}/lib${LIBFVAED2}.a ${OBJECTS}
	$(CC) ${SHARED} -o $@.${VERS} ${OBJECTS} ${LDFLAGS} ${SOFLAGS}
	ln -sf ${OUTLIB}.${so_ext}.${VERS} $@

${objdir}/%.o: ${srcdir}/%.F90 ${AED2DIR}/include/aed2.h
	$(F90) ${FFLAGS} ${INCLUDES} -g -c $< -o $@

${objdir}/tuflowfv_external_wq.o: tuflowfv_external_wq/tuflowfv_external_wq.f90
	$(FC) ${FFLAGS} ${TFFLAGS} ${INCLUDES} -Ituflowfv_external_wq -c $< -o $@

${objdir}:
	@mkdir ${objdir}

${moddir}:
	@mkdir ${moddir}

${libdir}:
	@mkdir ${libdir}

clean:
	/bin/rm -f *.i90
	/bin/rm -f ${objdir}/*.o
	/bin/rm -f ${moddir}/*.mod
	/bin/rm -f ${libdir}/*.a
	/bin/rm -f ${libdir}/*.${so_ext}*

distclean: clean
	/bin/rm -rf ${libdir} ${moddir} ${objdir} mod_s
