#!/bin/bash

# Change this as needed
CTAPATH=$HOME/cta
#CTAPATH=`pwd`

TMPDIR=/var/tmp

# Executables
CTA=$CTAPATH/calculate_torsion_angle
GETCHAIN=pdbgetchain

# Input PDB file
IN=$1

if [ "X$IN" == "X" ]; then
   echo "Usage: cta.sh in.pdb"
   exit 0;
fi

# Run program
LIGHT="$TMPDIR/L$$"
HEAVY="$TMPDIR/H$$"
$GETCHAIN L $IN >$LIGHT
$GETCHAIN H $IN >$HEAVY

RESULT=`$CTA -l $LIGHT -h $HEAVY | grep Torsion | awk '{print $3}'`

echo $RESULT

# Cleanup
\rm -f $LIGHT $HEAVY
