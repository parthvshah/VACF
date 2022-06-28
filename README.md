# Velocity Auto-Correlation Function (VACF)

Serial and parallel implementation of Velocity Auto-Correlation Function in C with MPICH for the paper titled _"Calculation of the Velocity Auto-Correlation Function Using Parallel Algorithms"_ by Parth Vipul Shah, Shubhadeep Nag, Debnath Pal, Subramanian Yashonath. Submitted for review.

## Introduction

Abstract: WIP

TODO: Update

This repository contains all the resources for the reproduction of our results which are detailed in our paper. The data required can be generated using [LAMMPS](https://www.lammps.org/).

## Setup

To run these algorithms on a cluster, you will need [MPICH](https://www.mpich.org/), an implementation of the Message Passing Interface (MPI). To generate data, you will need LAMMPS.

#### MPICH

MPICH can be installed from [here](https://www.mpich.org/downloads/) or a package manager like apt. Documentation is available [here](https://www.mpich.org/documentation/guides/). Once installed, run `hello_world.c` to verify installation.

Compile and run with:

```
$mpicc -o test hello_world.c -std=c99

$mpiexec.hydra -n 4 ./test
```

More information can be found on the official MPICH documentation.

`new_script.sh` is a sample script provided for your reference. We used this to submit jobs to our PBS queue running on a 120 node cluster with Intel Xeon E5-2670 (Haswell) CPUs.

#### LAMMPS

LAMMPS can be downloaded from [here](https://www.lammps.org/download.html). Unpack the same and run. The Manual is available [here](https://docs.lammps.org/Manual.html).

Run a simulation using:
```
$(lammps executable) < (lammps input script)
```

## Usage

All programs take 3 arguments during execution. The massively parallel implementations take one additional parameter.

```
Flags
---------------------------
-p int,int,int details the start,stop, step timesteps
-a int details the number of particles in the system
-i int details the number of correlations that need to be calculated
---------------------------
[-bp int] (optional) details the number of particles in one batch for the massively parallel algorithms
```

#### Example

To compile and run `Mpar1.3.1.c` on 3 nodes with 24 processors each:

```
$mpicc -o vacf.Mpar1.3.1 Mpar1.3.1.c -std=c99

$time mpiexec.hydra -np 72 -genvall -ppn 24 ./vacf.Mpar1.3.1 -p 50000,9857960,10 -a 6912 -i 1000 -bp 144 > out.data 2>&1
```

## Organization

This project has multiple directories. Here is a brief description of each. Please navigate to the specific directories for more information.

- algos

    This directory contains all the algorithms, serial and parallel.

- lammps

    This directory contains the details for the data generation using LAMMPS, a widely used parallel software for Molecular Dynamic simulations.

- timing

    This directory contains CSV's with timing data. For more (and complete) timing information, visit [this](https://github.com/parthvshah/VACFTimingData) repository or refer to the tables in our paper.

- utils

    This directory contains the utility functions associated with this project. They are python scripts.

## Citing

Cite this work as: WIP

Please contact the authors for any additional information required. All rights reserved.
