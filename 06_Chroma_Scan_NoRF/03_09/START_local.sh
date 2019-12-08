#!/bin/bash

set -o errexit

if [ ! -n "$1" ]
  then
    echo "Usage: `basename $0` <name of the SC script> <N CPUs>"
    exit $E_BADARGS
fi

if [ ! -n "$2" ]
  then
    echo "Usage: `basename $0` <name of the SC script> <N CPUs>"
    exit $E_BADARGS
fi

export FI_PROVIDER=sockets
. setup_environment.sh

mpirun -np $2 ${ORBIT_ROOT}/bin/pyORBIT $1
