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

#define ROW 1024
#define COL 5120

double xData[ROW][COL];
double yData[ROW][COL];
double zData[ROW][COL];

void padding(int n)
{
    for (int i = 0; i < n; i++)
    {
        xData[0][i] = yData[0][i] = zData[0][i] = 0.0;
    }
}

void readData(int start, int step, char *fileName)
{
    int timestep, particle, index;
    double xVel, yVel, zVel;
    FILE *fp;

    fp = fopen(fileName, "r");
    if (fp == NULL)
        return;

    char line[128];
    while (fgets(line, sizeof line, fp) != NULL)
    {
        sscanf(line, "%d %d %lf %lf %lf", &timestep, &particle, &xVel, &yVel, &zVel);
        index = (timestep - start) / step;

        xData[particle][index] = xVel;
        yData[particle][index] = yVel;
        zData[particle][index] = zVel;
    }
    fclose(fp);
}

int main(int argc, char **argv)
{
    int rank, wSize, NStart, NEnd, chunk, remainder;
    double output[2];
    MPI_Status status;
    MPI_File fp;

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
        char fileName[256];
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

        padding(M);
        readData(start, step, fileName);

    }

    MPI_Bcast(&tmax, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&M, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&N, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&xData, ROW * COL, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&yData, ROW * COL, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&zData, ROW * COL, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    double *lCorr = (double *)malloc(sizeof(double) * (tmax + 1));
    lCorr[0] = 0.0;

    double *gCorr = (double *)malloc(sizeof(double) * (tmax + 1));
    gCorr[0] = 0.0;

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

    MPI_Reduce(lCorr, gCorr, tmax + 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

    MPI_File_open(MPI_COMM_WORLD, "OUT", MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &fp);
    if (rank == 0)
    {
        for (int k = 1; k <= tmax; k++)
        {
            output[0] = (double)k;
            output[1] = gCorr[k];
            MPI_File_write(fp, &output, 2, MPI_DOUBLE, MPI_STATUS_IGNORE);
            // fprintf(stdout, "%d, %e\n", k, gCorr[k]);
        }
    }

    MPI_File_close(&fp);
    MPI_Finalize();
    return 0;
}