/*
VACF (Correlation Decompositon with Load Balancing)
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
    Calculate localStart1, localStart2, localStop1, localStop2
    for dt <- localStart1 to localStop1
        count = 0
        accumalate = 0
        for t <- 0 to (M-dt)
            particle = 0.0
            for i <- 1 to N
                particle += Vx[i][t].Vx[i][t+dt] + Vy[i][t].Vy[i][t+dt] + Vz[i][t].Vz[i][t+dt]
            count += 1
            accumalate += particle
        print(dt, accumalate/((N-1)*count)
    for dt <- localStart2 to localStop2
        count = 0
        accumalate = 0
        for t <- 0 to (M-dt)
            particle = 0.0
            for i <- 1 to N
                particle += Vx[i][t].Vx[i][t+dt] + Vy[i][t].Vy[i][t+dt] + Vz[i][t].Vz[i][t+dt]
            count += 1
            accumalate += particle
        print(dt, accumalate/((N-1)*count)
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>

#define ROW 2064
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

    // "./HISTORY_atoms/HISTORY_CLEAN_864"
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
    int rank, wSize, lStart1, lEnd1, lStart2, lEnd2, chunk, remainder;
    double batch1Start, batch1End, batch2Start, batch2End;
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

        if (tmax % wSize != 0)
        {
            tmax = tmax - (tmax % wSize);
        }

    }

    MPI_Bcast(&M, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&N, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&tmax, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&xData, ROW * COL, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&yData, ROW * COL, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&zData, ROW * COL, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    chunk = tmax / (wSize * 2);
    remainder = tmax - (chunk * wSize * 2);
    lStart1 = (rank * chunk) + 1;
    lEnd1 = lStart1 + chunk - 1;

    MPI_File_open(MPI_COMM_WORLD, "OUT", MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &fp);
    int viewLength = chunk;
    MPI_File_set_view(fp, (rank * sizeof(double) * 2 * viewLength), MPI_INT, MPI_INT, "native", MPI_INFO_NULL);

    for (int dt = lStart1; dt <= lEnd1; dt++)
    {
        count = 0;
        accumalate = 0.0;
        for (int t = 0; t < (M - dt); t++)
        {
            particle = 0.0;
            for (int i = 1; i < N; i++)
            {
                particle += xData[i][t] * xData[i][t + dt] +
                            yData[i][t] * yData[i][t + dt] +
                            zData[i][t] * zData[i][t + dt];
            }
            count++;
            accumalate += particle;
        }
        accumalate /= ((N - 1) * count);
        output[0] = (double)dt;
        output[1] = accumalate;
        MPI_File_write_all(fp, &output, 2, MPI_DOUBLE, MPI_STATUS_IGNORE);
        // fprintf(stdout, "%d, %e\n", dt, accumalate);
    }

    lEnd2 = (tmax - remainder) - (rank * chunk);
    lStart2 = lEnd2 - chunk + 1;

    if (rank == 0)
    {
        lEnd2 += remainder;
    }
    // Batch2
    for (int dt = lStart2; dt <= lEnd2; dt++)
    {
        count = 0;
        accumalate = 0.0;
        for (int t = 0; t < (M - dt); t++)
        {
            particle = 0.0;
            for (int i = 1; i < N; i++)
            {
                particle += xData[i][t] * xData[i][t + dt] +
                            yData[i][t] * yData[i][t + dt] +
                            zData[i][t] * zData[i][t + dt];
            }
            count++;
            accumalate += particle;
        }
        accumalate /= ((N - 1) * count);
        output[0] = (double)dt;
        output[1] = accumalate;
        MPI_File_write_all(fp, &output, 2, MPI_DOUBLE, MPI_STATUS_IGNORE);
        // fprintf(stdout, "%d, %e\n", dt, accumalate);
    }

    MPI_File_close(&fp);
    MPI_Finalize();

    return 0;
}