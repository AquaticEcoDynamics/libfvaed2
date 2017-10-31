#!/bin/bash

if [ "$FV_CONFIGURED" != "true" ] ; then
  . ./FV_CONFIG
fi

if [ "$AED2DIR" = "" ] ; then
  export AED2DIR=../libaed2
fi

cd ${AED2DIR}
make distclean
cd ${CURDIR}

if [ "${AED2PLS}" != "" ] ; then
   if [ -d ${AED2PLS} ] ; then
      cd ${AED2PLS}
      make distclean
      cd ${CURDIR}
   fi
fi

make distclean
