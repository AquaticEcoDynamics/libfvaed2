#
# Makefile to build the aed2 water quality library
# with hydrodynamic driver wrapper
#

objdir=obj
moddir=mod
srcdir=src
libdir=lib
VERS=1.0.1

ifeq ($(EXTERNAL_LIBS),static)
  tuflowfv_external_wq=libtuflowfv_external_wq.a
else
  tuflowfv_external_wq=libtuflowfv_external_wq.so
endif

ifeq ($(AED2DIR),)
  AED2DIR=../libaed2
endif
OUTLIB=libtuflowfv_wq_aed


INCLUDES+=-I${AED2DIR}/include
ifeq ($(SINGLE),true)
  INCLUDES+=-I${AED2DIR}/mod_s
  LIBAED2=libaed2_s
else
  INCLUDES+=-I${AED2DIR}/mod
  LIBAED2=libaed2
endif

ifeq ($(F90),ifort)
  INCLUDES+=-I/opt/intel/include
  DEBUG_FFLAGS=-g -traceback
  OPT_FFLAGS=-O3 -openmp
  FFLAGS=-warn all -module ${moddir} -i-static -mp1 -warn nounused $(DEFINES)
  ifeq ($(WITH_CHECKS),true)
    FFLAGS+=-check
  endif
  ifeq ($(SINGLE),true)
    FFLAGS+=-real-size 32 -DSINGLE=1
  else
    FFLAGS+=-real-size 64
  endif
  LIBS+=-lifcore -lsvml
  LIBS+=-limf -lintlc
  LIBS+=-L/opt/intel/lib -Wl,-rpath=/opt/intel/lib
endif

FFLAGS+=-fPIC
LIBS += -L${AED2DIR} -laed2

LIBS+=-lnetcdff -lnetcdf

FFLAGS+=$(OPT_FFLAGS)

ifeq ($(PRECISION),1)
  TFFLAGS += -D_PRECISION=1
else ifeq ($(PRECISION),2)
  TFFLAGS += -D_PRECISION=2
else
  TFFLAGS += -D_PRECISION=1
endif

TFFLAGS += -g -DAED2 -DEXTERNAL_WQ=2
INCLUDES += -I${moddir}


OBJECTS=${objdir}/fv_zones.o ${objdir}/fv_aed2_csv_reader.o ${objdir}/fv_aed2.o ${objdir}/tuflowfv_external_wq.o


ifeq ($(EXTERNAL_LIBS),static)

all: ${objdir} ${moddir} ${libdir} ${libdir}/$(OUTLIB).a

else

  LDFLAGS += --start-group ${AED2DIR}/lib/${LIBAED2}.a --end-group

all: ${objdir} ${moddir} ${libdir} ${libdir}/$(OUTLIB).so

endif

${libdir}/$(OUTLIB).a: $(OBJECTS)
	ar -rv $@ $(OBJECTS) $(LDFLAGS)
	ranlib $@

${libdir}/$(OUTLIB).so: $(OBJECTS)
	ld -shared -o $@.$(VERS) $(OBJECTS) $(LDFLAGS)
	ln -sf $@.$(VERS) $@


${objdir}/%.o: ${srcdir}/%.F90 ${AED2DIR}/include/aed2.h
	$(F90) $(FFLAGS) $(INCLUDES) -g -c $< -o $@

${objdir}/tuflowfv_external_wq.o: tuflowfv_external_wq/tuflowfv_external_wq.f90
	$(FC) $(FFLAGS) $(TFFLAGS) $(INCLUDES) -Ituflowfv_external_wq -c $< -o $@


${objdir}:
	@mkdir ${objdir}

${moddir}:
	@mkdir ${moddir}

${libdir}:
	@mkdir ${libdir}


clean:
	/bin/rm -f ${objdir}/*.o
	/bin/rm -f ${moddir}/*.mod
	/bin/rm -f ${libdir}/*.a

distclean: clean
	/bin/rm -rf ${libdir} ${moddir} ${objdir}
