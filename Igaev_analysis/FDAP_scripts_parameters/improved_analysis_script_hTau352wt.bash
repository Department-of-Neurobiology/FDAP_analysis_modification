#!/bin/bash

NCPU=6

if [ $# -eq 1 ]
then
  NCPU=$1
fi

echo "using ${NCPU} CPUs"

i=`expr 0`

for fn in ./*.dat; do
	echo Analysing the file ${fn}...
	mpiexec -n $NCPU ./cFDAP -m fullModel -d 12.6001 -kon0 0.5 -koff0 0.5 -i ${fn} -o result_$i >> output.log
	echo Done with the file ${fn}. Taking the next one...

	i=`expr $i + 1`
done
