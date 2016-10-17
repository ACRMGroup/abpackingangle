#!/bin/bash

TMPDIR=/var/tmp

# Executables
CTA=../abpackingangle
GETCHAIN=pdbgetchain

# Input PDB file
IN=8FAB_1.pdb

# Run program
LIGHT="$TMPDIR/L$$"
HEAVY="$TMPDIR/H$$"
$GETCHAIN L $IN >$LIGHT
$GETCHAIN H $IN >$HEAVY

RESULT=`$CTA -l $LIGHT -h $HEAVY | grep Torsion | awk '{print $3}'`

echo $RESULT >8FAB_1_test.out

echo "you should see nothing after this line!"
diff 8FAB_1.out 8FAB_1_test.out



# Cleanup
\rm -f $LIGHT $HEAVY 8FAB_1_test.out

