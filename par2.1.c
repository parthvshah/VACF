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

Date: 10-6-2019
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
    int wRank, wSize, lChunk, NChunk, lRemainder, NRemainder, NStart, NEnd, lStart, lEnd;
    double output[2];
    MPI_Status status;
    MPI_File fp;

    MPI_Init(&argc, &argv);

    MPI_Comm_rank(MPI_COMM_WORLD, &wRank);
    MPI_Comm_size(MPI_COMM_WORLD, &wSize);

    int colour = wRank / 12;
    MPI_Comm MPI_COMM_ROW, MPI_COMM_ROWROOT;
    MPI_Comm_split(MPI_COMM_WORLD, colour, wRank, &MPI_COMM_ROW);

    int rRank, rSize;
    MPI_Comm_rank(MPI_COMM_ROW, &rRank);
    MPI_Comm_size(MPI_COMM_ROW, &rSize);

    int *ranks = (int *)malloc(sizeof(int)*wSize/12);
    int iCount = 0;
    for(int iRank = 0; iRank<wSize; iRank++)
        if(iRank%12==0)
            ranks[iCount++] = iRank;
    MPI_Group MPI_GROUP_WORLD, MPI_GROUP_ROWROOT;
    MPI_Comm_group(MPI_COMM_WORLD, &MPI_GROUP_WORLD);
    MPI_Group_incl(MPI_GROUP_WORLD, (wSize / 12), ranks, &MPI_GROUP_ROWROOT);
    MPI_Comm_create(MPI_COMM_WORLD, MPI_GROUP_ROWROOT, &MPI_COMM_ROWROOT);

    int start;
    int stop;
    int step;
    int N, M;
    int tmax = 0;

    double accumalate, particle;
    int count;

    if (wRank == 0)
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

    MPI_Bcast(&M, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&N, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&tmax, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&xData, ROW * COL, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&yData, ROW * COL, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&zData, ROW * COL, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    double *lCorr = (double *)malloc(sizeof(double) * (tmax + 1));
    lCorr[0] = 0.0;

    double *gCorr = (double *)malloc(sizeof(double) * (tmax + 1));
    gCorr[0] = 0.0;

    // Correlation calculation
    lChunk = tmax / (wSize / rSize);
    lRemainder = tmax - (lChunk * (wSize / rSize));
    lStart = ((wRank - (wRank % 12)) / 12) * lChunk + 1;
    lEnd = lStart + lChunk - 1;

    if ((lRemainder != 0) && (((wRank - (wRank % 12)) / 12) == 12))
    {
        lEnd += lRemainder;
    }

    // Atom calculation
    NChunk = (N - 1) / rSize;
    NRemainder = (N - 1) - (NChunk * rSize);
    NStart = (rRank * NChunk) + 1;
    NEnd = NStart + NChunk - 1;

    if ((NRemainder != 0) && (rRank == (rSize - 1)))
    {
        NEnd += NRemainder;
    }

    MPI_File_open(MPI_COMM_ROWROOT, "OUT", MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &fp);

    if(rRank == 0)
    {
        int viewLength = lChunk;
        MPI_File_set_view(fp, (((wRank - (wRank % 12)) / 12) * sizeof(double) * 2 * viewLength), MPI_INT, MPI_INT, "native", MPI_INFO_NULL);
    }

    for (int dt = lStart; dt <= lEnd; dt++)
    {
        count = 0;
        accumalate = 0.0;
        for (int t = 0; t < (M - dt); t++)
        {
            particle = 0.0;
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

    MPI_Reduce(lCorr, gCorr, tmax + 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_ROW);
    
    if (rRank == 0)
    {
        for (int k = lStart; k <= lEnd; k++)
        {
            output[0] = (double)k;
            output[1] = gCorr[k];
            // fprintf(stdout, "%d, %e\n", k, gCorr[k]);
            MPI_File_write_all(fp, &output, 2, MPI_DOUBLE, MPI_STATUS_IGNORE);
        }
    }

    MPI_File_close(&fp);
    MPI_Finalize();
    
    return 0;
}