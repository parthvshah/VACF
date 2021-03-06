/*
VACF (Particle decomposition)
Command line arguments:
    -p START,STOP,STEP
    -a PARTICLES
    -i ITERATIONS
    -f FILENAME

Algorithm:
    size = no. of processes
    M = no. of timesteps
    N = no. of particles
    Vx, Vy, Vz = Velocity of a particle at timestep 
    (0)
    Parse command line arguments
    Read data
    (0-size)
    Receive broadcasted information from 0
    Calculate atoms in process
    for dt <- 1 to tmax
        count = 0
        accumalate = 0
        for t <- 0 to (M-dt)
            particle = 0.0
            for i <- localNStart to localNStop
                particle += Vx[i][t].Vx[i][t+dt] + Vy[i][t].Vy[i][t+dt] + Vz[i][t].Vz[i][t+dt]
            count += 1
            accumalate += particle
        localCorrelation[dt] = accumalate/((N-1)*count
    (0)
    Accumalate localCorrelation from each process and add it to globalCorrelation
    for i <- 1 to tmax
        print(i, globalCorrelation[i])

 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>

#define ROW 501
#define COL 20000001


void padding(int n, float **xData, float **yData, float **zData)
{
    for (int i = 0; i < n; i++)
    {
        xData[0][i] = yData[0][i] = zData[0][i] = 0.0;
    }
}

void readData(int start, int stop, int step, int N, float **xData, float **yData, float** zData)
{
    int timestep, particle, one, index, count, maxCount;
    float xVel, yVel, zVel;
    FILE *fp;

    fp = fopen("./HISTORY_atoms/HISTORY_CLEAN_500_l", "r");
    if (fp == NULL)
        return;
    
    count = 1;
    maxCount = (stop-start) / step * N;

    char line[256];
    while (fgets(line, sizeof line, fp) != NULL)
    {
        if(count == maxCount)
            return;
        sscanf(line, "%d %d %d %f %f %f", &timestep, &particle, &one, &xVel, &yVel, &zVel);
        index = (timestep - start) / step;

        xData[particle][index] = xVel;
        yData[particle][index] = yVel;
        zData[particle][index] = zVel;

        count++;
    }
    fclose(fp);
}

int main(int argc, char **argv)
{
    int rank, wSize, NStart, NEnd, chunk, remainder;
    MPI_Status status;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &wSize);

    int start;
    int stop;
    int step;
    int N, M;
    int tmax = 0;

    double accumalate, particle;
    int count;

    if (rank == 0)
    {
        char fileName[128];
        for (int c = 1; c < argc; c++)
        {
            if (!strcmp(argv[c], "-p"))
                sscanf(argv[++c], "%d,%d,%d", &start, &stop, &step);
            else if (!strcmp(argv[c], "-a"))
                N = atoi(argv[++c]);
            else if (!strcmp(argv[c], "-i"))
                tmax = atoi(argv[++c]);
            else if (!strcmp(argv[c], "-f"))
                sscanf(argv[++c], "%s", fileName);
            else
            {
                fprintf(stderr, "[Error] Command-line argument not recognized.\n");
                exit(-1);
            }
        }

        N++;                             // Number of particles
        M = ((stop - start) / step) + 1; // Number of timesteps

        tmax = (tmax == 0) ? (M / 3) : tmax;

        // padding(M);
        // readData(start, step, fileName);

    }

    float **xData = NULL;
    xData = (float**) malloc(sizeof(float*) * ROW);
    if(!xData)
    {
        free(xData);
        fprintf(stdout, "[Error - %d] xData not allocated.\n", rank);
        return EXIT_FAILURE;
    }
    for (int ii = 0; ii<ROW; ii++)
    {
        xData[ii] = NULL;
        xData[ii] = (float*) malloc(sizeof(float) * COL);
        if(!xData[ii])
        {
            fprintf(stdout, "[Error - %d] Internal xData not allocated. \n", rank);
            for(int jj = ii; jj>=0; jj--)
                free(xData[jj]);
            
            return EXIT_FAILURE;
        }
    }

    
    float **yData = NULL;
    yData = (float**) malloc(sizeof(float*) * ROW);
    if(!yData)
    {
        free(yData);
        fprintf(stdout, "[Error - %d] xData not allocated.\n", rank);
        return EXIT_FAILURE;
    }
    for (int ii = 0; ii<ROW; ii++)
    {
        yData[ii] = NULL;
        yData[ii] = (float*) malloc(sizeof(float) * COL);
        if(!yData[ii])
        {
            fprintf(stdout, "[Error - %d] Internal xData not allocated. \n", rank);
            for(int jj = ii; jj>=0; jj--)
                free(yData[jj]);
                
            return EXIT_FAILURE;
        }
    }

    
    float **zData = NULL;
    zData = (float**) malloc(sizeof(float*) * ROW);
    if(!zData)
    {
        free(zData);
        fprintf(stdout, "[Error - %d] xData not allocated.\n", rank);
        return EXIT_FAILURE;
    }
    for (int ii = 0; ii<ROW; ii++)
    {
        zData[ii] = NULL;
        zData[ii] = (float*) malloc(sizeof(float) * COL);
        if(!zData[ii])
        {
            fprintf(stdout, "[Error - %d] Internal xData not allocated. \n", rank);
            for(int jj = ii; jj>=0; jj--)
                free(zData[jj]);
            
            return EXIT_FAILURE;
        }
    }

    MPI_Bcast(&start, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&stop, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&step, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&M, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&N, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&tmax, 1, MPI_INT, 0, MPI_COMM_WORLD);

    float *lCorr = NULL;
    lCorr = (float *)malloc(sizeof(float) * (tmax + 1));
    if(!lCorr)
    {
        free(lCorr);
        fprintf(stdout, "[Error - %d] lCorr not allocated.\n", rank);
        return EXIT_FAILURE;
    }

    float *gCorr = NULL;
    gCorr = (float *)malloc(sizeof(float) * (tmax + 1));
    if(!gCorr)
    {
        free(gCorr);
        fprintf(stdout, "[Error - %d] gCorr not allocated.\n", rank);
        return EXIT_FAILURE;
    }

    padding(M, xData, yData, zData);
    readData(start, stop, step, N, xData, yData, zData);

    // Atom calculation
    chunk = (N - 1) / wSize;
    remainder = (N - 1) - (chunk * wSize);
    NStart = (rank * chunk) + 1;
    NEnd = NStart + chunk - 1;

    if ((remainder != 0) && (rank == (wSize - 1)))
    {
        NEnd += remainder;
    }

    for (int dt = 1; dt <= tmax; dt++)
    {
        count = 0;
        accumalate = 0.0;
        for (int t = 0; t < (M - dt); t++)
        {
            particle = 0.0;
            // for (int i = 1; i < N; i++)
            for (int i = NStart; i <= NEnd; i++)
            {
                particle += xData[i][t] * xData[i][t + dt] +
                            yData[i][t] * yData[i][t + dt] +
                            zData[i][t] * zData[i][t + dt];
            }
            count++;
            accumalate += particle;
        }
        accumalate /= ((N - 1) * count);
        lCorr[dt] = accumalate;
        // fprintf(stdout, "%d, %e\n", dt, accumalate);
    }

    MPI_Reduce(lCorr, gCorr, tmax + 1, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);

    if (rank == 0)
    {
        for (int k = 1; k <= tmax; k++)
        {
            fprintf(stdout, "%d, %e\n", k, gCorr[k]);
        }
    }

    MPI_Finalize();
    return 0;
}