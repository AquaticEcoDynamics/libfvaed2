#!/bin/bash -e

if [ "$FV_CONFIGURED" != "true" ] ; then
  . ./FV_CONFIG
fi

if [ "$FC" = "" ] ; then
  export FC=ifort
fi

if [ "$FC" = "ifort" ] ; then
  if [ -d /opt/intel/bin ] ; then
    . /opt/intel/bin/compilervars.sh intel64
  fi
  which ifort >& /dev/null
  if [ $? != 0 ] ; then
    echo ifort compiler requested, but not found
    exit 1
  fi
fi

if [ "$AED2DIR" = "" ] ; then
  export AED2DIR=../libaed2
fi

if [ "$PRECISION" = "" ] ; then
  export PRECISION=1
fi
if [ "$EXTERNAL_LIBS" = "" ] ; then
  export EXTERNAL_LIBS=shared
fi
if [ "$DEBUG" = "" ] ; then
  export DEBUG=false
fi

export F77=$FC
export F90=$FC
export F95=$FC

cd ${AED2DIR}
make
cd ${CURDIR}

if [ "${AED2PLS}" != "" ] ; then
   if [ -d ${AED2PLS} ] ; then
      cd ${AED2PLS}
      make
      cd ${CURDIR}
   fi
fi

#make distclean
make
