
OUTLIB=libtuflowfv_external_wq
VERS=1
FC=ifort
SRCDIR=.
INCLUDES=
LDFLAGS=

ifndef EXTERNAL_LIBS
  EXTERNAL_LIBS = static
endif

ifeq ($(COMPILATION_MODE),debug)
  FFLAGS=-g -traceback -fpp -check -fPIC
  buildt=_debug
else
  FFLAGS=-g -axSSE4.1 -fpp -fPIC
  buildt=_prod
endif

ifeq ($(PRECISION),1)
  FFLAGS += -D_PRECISION=1
else ifeq ($(PRECISION),2)
  FFLAGS += -D_PRECISION=2
else
  FFLAGS += -D_PRECISION=1
endif
MODDIR=mod
AEDLIB=libaed2.a
FVAEDLIB=fvaed2

ifdef AED2DIR
    OBJECTS=
    FFLAGS += -g -DAED2 -DEXTERNAL_WQ=2
    INCLUDES += -I${AED2DIR}/$(MODDIR) -I$(FVAED2DIR)/$(AED2MOD)
    OBJECTS+=tuflowfv_external_wq.o
    ifeq ($(EXTERNAL_LIBS),shared)
      LDFLAGS += --start-group ${AED2DIR}/lib/$(AEDLIB) ../tuflowfv_aed/$(FVAEDLIB) --end-group
    endif
else
    OBJECTS=tuflowfv_external_wq.o
endif

ifeq ($(EXTERNAL_LIBS),static)
all: $(OBJECTS)
	ar -rv $(OUTLIB).a $(OBJECTS) $(LDFLAGS)
else
all: $(OBJECTS)
	ld -shared -o $(OUTLIB).so.$(VERS) $(OBJECTS) $(LDFLAGS)
	ln -sf $(OUTLIB).so.$(VERS) $(OUTLIB).so
endif

%.o: $(SRCDIR)/%.f90
	$(FC) $(FFLAGS) $(INCLUDES) -c $<

clean:
	rm -f *.mod *.o *.a

distclean:
	rm -f *.mod *.o *.a *.so*

