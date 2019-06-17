#!/bin/bash

# --------------------------------------------------------------------
# System settings
# --------------------------------------------------------------------

TIMEFORMAT=%R
HISTORY_CLEAN_FILENAME="HISTORY_CLEAN"
HOSTFILE_FILENAME="host_file"
OUT_FILENAME="OUT"

echo "VACF AND COEFFICIENT OF DISSFUSION"

# --------------------------------------------------------------------
# User modified variables
# --------------------------------------------------------------------

FILE="HISTORY"
START=50000
STOP=86650
STEP=10
PARTICLES=864
ITERATIONS=84

CLEAN=false
PLOT=false
PARALLEL=true
LOCAL=false

THREADS=4

GRAPH_LOWER=0
GRAPH_UPPER=166

# --------------------------------------------------------------------
# Commands
# --------------------------------------------------------------------

if [ "$CLEAN" = true ] ; then
    echo "[INFO] Cleaning data"
    time python3 clean.py > $HISTORY_CLEAN_FILENAME
fi

if [ "$PARALLEL" = true ]; then
    if [ "$LOCAL" = true ]; then
        echo "[INFO] Generating VACF plot (local, parallel) for $ITERATIONS iterations."
        mpicc -o vacf par1.1.c -std=c99
        time mpirun -n $THREADS ./vacf -p $START,$STOP,$STEP -a $PARTICLES -i $ITERATIONS > $OUT_FILENAME
    else
        echo "[INFO] Generating VACF plot (parallel) for $ITERATIONS iterations."
        mpicc -o vacf par1.1.c -std=c99
        time mpirun -hostfile $HOSTFILE_FILENAME ./vacf -p $START,$STOP,$STEP -a $PARTICLES -i $ITERATIONS > $OUT_FILENAME
    fi
else
    echo "[INFO] Generating VACF plot (series) for $ITERATIONS iterations."
    gcc -o vacf seq.c -g -O3 -std=c99
    time ./vacf -p $START,$STOP,$STEP -a $PARTICLES -i $ITERATIONS > $OUT_FILENAME
    
fi

if [ "$PLOT" = true ] ; then
    python3 plot.py $GRAPH_LOWER $GRAPH_UPPER
fi

echo "[INFO] Calculating coefficient of diffusion"
time python3 diffusion.py
