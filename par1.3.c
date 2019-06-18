/*
Calculate Velocity Auto-Correlation (Parallel)
This program calculates the VACF from a given HISTORY file. Scripts to clean the data,
plot the graph and calculate the coefficient of diffusion are external and in python.

Command line arguments:
    -p START,STOP,STEP
    -a PARTICLES
    -i ITERATIONS

Algorithm:
    size = no. of processes
    M = no. of timesteps
    N = no. of particles
    Vx, Vy, Vz = Velocity of a particle at timestep
    tmax = M/3 
    Parse command line arguments
    Read data
    for dt <- 1 to tmax
        count = 0
        accumalate = 0
        for t <- 0 to (M-dt)
            particle = 0.0
            for i <- 1 to N
                particle += Vx[i][t].Vx[i][t+dt] + Vy[i][t].Vy[i][t+dt] + Vz[i][t].Vz[i][t+dt]
            count += 1
            accumalate += particle
        print(dt, accumalate/((N-1)*count)

Date: 13-6-2019
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
    int rank, wSize, lStart, lEnd, chunk, remainder;
    MPI_Status status;

    MPI_Group orig_group, new_group;
    MPI_Comm MPI_COMM_SLAVES;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &wSize);

    // int slaveRanks[3] = {1, 2, 3};

    int *slaveRanks = (int *)malloc(sizeof(int) * (wSize - 1));
    for (int iS = 0; iS < (wSize - 1); iS++)
    {
        slaveRanks[iS] = iS + 1;
    }

    MPI_Comm_group(MPI_COMM_WORLD, &orig_group);
    MPI_Group_incl(orig_group, (wSize - 1), slaveRanks, &new_group);
    MPI_Comm_create(MPI_COMM_WORLD, new_group, &MPI_COMM_SLAVES);

    int start;
    int stop;
    int step;
    int N, M;
    int t, dt;
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

        fprintf(stdout, "timestep, vacf\n");

        chunk = (N - 1) / wSize;
        remainder = (N - 1) - (wSize * chunk);

        int sendStart, sendEnd;

        lStart = 1;
        lEnd = lStart - 1 + chunk;

        for (int process = 1; process < wSize - 1; process++)
        {
            sendStart = (process * chunk) + 1;
            sendEnd = sendStart - 1 + chunk;
            MPI_Send(&sendStart, 1, MPI_INT, process, 101, MPI_COMM_WORLD);
            MPI_Send(&sendEnd, 1, MPI_INT, process, 102, MPI_COMM_WORLD);
        }
        sendStart = ((wSize - 1) * chunk) + 1;
        sendEnd = sendStart - 1 + chunk + remainder;
        MPI_Send(&sendStart, 1, MPI_INT, wSize - 1, 101, MPI_COMM_WORLD);
        MPI_Send(&sendEnd, 1, MPI_INT, wSize - 1, 102, MPI_COMM_WORLD);
    }

    if (rank != 0)
    {
        MPI_Recv(&lStart, 1, MPI_INT, 0, 101, MPI_COMM_WORLD, &status);
        MPI_Recv(&lEnd, 1, MPI_INT, 0, 102, MPI_COMM_WORLD, &status);
    }

    MPI_Bcast(&xData, ROW * COL, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&yData, ROW * COL, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&zData, ROW * COL, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    if (rank == 0)
    {
        for (int dt = 1; dt <= tmax; dt++)
        {
            count = 0;
            accumalate = 0.0;
            for (int t = 0; t < (M - dt); t++)
            {
                particle = 0.0;
                for (int sendIndex = 1; sendIndex < wSize; sendIndex++)
                {
                    MPI_Send(&t, 1, MPI_INT, sendIndex, 201, MPI_COMM_WORLD);
                    MPI_Send(&dt, 1, MPI_INT, sendIndex, 202, MPI_COMM_WORLD);
                }

                for (int i = lStart; i <= lEnd; i++)
                {
                    accumalate += xData[i][t] * xData[i][t + dt] +
                                yData[i][t] * yData[i][t + dt] +
                                zData[i][t] * zData[i][t + dt];
                }

                for (int recIndex = 1; recIndex < wSize; recIndex++)
                {
                    MPI_Recv(&particle, 1, MPI_DOUBLE, recIndex, 203, MPI_COMM_WORLD, &status);
                    accumalate += particle;
                }

                count++;
                // accumalate += particle;
            }
            accumalate /= ((N - 1) * count);
            fprintf(stdout, "%d, %e\n", dt, accumalate);
        }
    }

    if (rank != 0)
    {
        while (1)
        {
            MPI_Recv(&t, 1, MPI_INT, 0, 201, MPI_COMM_WORLD, &status);
            MPI_Recv(&dt, 1, MPI_INT, 0, 202, MPI_COMM_WORLD, &status);

            particle = 0.0;
            for (int i = lStart; i <= lEnd; i++)
            {
                particle += xData[i][t] * xData[i][t + dt] +
                            yData[i][t] * yData[i][t + dt] +
                            zData[i][t] * zData[i][t + dt];
            }

            MPI_Send(&particle, 1, MPI_DOUBLE, 0, 203, MPI_COMM_WORLD);

            MPI_Barrier(MPI_COMM_SLAVES);
        }
    }

    exit(1);
    MPI_Finalize();

    return 0;
}