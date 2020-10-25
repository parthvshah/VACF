# VACF
Serial and parallel implementation of Velocity Auto-Correlation Function

## Naming
- ```Mpar1.1.c``` - [Massively Parallel] Timestep decomposition (Algorithm 1.1)
- ```Mpar1.3.1.c``` - Double decomposition (Algorithm 2.1)
- ```par1.1.c``` - Correlation decomposition, load balancing
- ```par1.3.1.c``` - Particle decomposition, lesser communication (Algorithm 2)
- ```par1.3.c``` - Particle decomposition, barrier implementation
- ```par1.c``` - Correlation decomposition (Algorithm 1)
- ```par2.1.c``` - Double decomposition, batch processing (Algorithm 3)
- ```par2.c``` - Double decomposition, sub-root processes
- ```parOpenMP.c``` - OpenMP implementation, correlation decomposition
- ```seq.c``` - Sequential implementation 
- ```seq1.c``` - Sequential with Redis read optimization
