#!/bin/bash

if [ "$FV_CONFIGURED" != "true" ] ; then
  . ./FV_CONFIG
fi

if [ "$AED2DIR" = "" ] ; then
  export AED2DIR=../libaed2
fi

if [ "$DEBUG" = "" ] ; then
   export DEBUG=true
fi
if [ "$SINGLE" = "" ] ; then
   export SINGLE=true
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

if [ "$SINGLE" = "true" ] ; then
  export PRECISION=1
else
  export PRECISION=2
fi

#make distclean
make
