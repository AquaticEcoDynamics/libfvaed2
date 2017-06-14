#!/bin/bash

if [ "$FV_CONFIGURED" != "true" ] ; then
  . ./FV_CONFIG
fi

if [ "$FORTRAN_COMPILER" = "IFORT" ] ; then
   if [ -d /opt/intel/bin ] ; then
      . /opt/intel/bin/compilervars.sh intel64
   fi
   which ifort >& /dev/null
   if [ $? != 0 ] ; then
      echo ifort compiler requested, but not found
      exit 1
   fi
fi

if [ "$FORTRAN_COMPILER" = "IFORT" ] ; then
   export PATH="/opt/intel/bin:$PATH"
   export FC=ifort
   export NETCDFHOME=/opt/intel
else
   export FC=gfortran
   export NETCDFHOME=/usr
   if [ "$OSTYPE" == "Darwin" ] ; then
     export NETCDFHOME=/opt/local
   fi
fi

export F77=$FC
export F90=$FC
export F95=$FC

export MPI=OPENMPI

if [ "$AED2DIR" = "" ] ; then
  export AED2DIR=../libaed2
fi

if [ "$DEBUG" = "" ] ; then
   export DEBUG=true
fi
if [ "$PRECISION" = "" ] ; then
   export PRECISION=1
fi
if [ "$EXTERNAL_LIBS" = "" ] ; then
  export EXTERNAL_LIBS=shared
fi
if [ "$FC" = "" ] ; then
  export FC=ifort
fi

export F77=$FC
export F90=$FC
export F95=$FC

cd ${AED2DIR}
make
cd ${CURDIR}

if [ -d ${AED2PLS} ] ; then
   cd ${AED2PLS}
   make
   cd ${CURDIR}
fi

#make distclean
make
