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
    int rank, wSize, lStart1, lEnd1, lStart2, lEnd2, chunk, remainder;
    double readStart, readEnd;
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

        if (tmax % wSize != 0)
        {
            tmax = tmax - (tmax % wSize);
        }

    }

    float **xData = (float **)malloc(ROW * sizeof(float *));
    if(xData==NULL)
    {
        fprintf(stderr, "[Error] xData allocation failed.\n");
        exit(1);
    } 
    for (int i=0; i<ROW; i++)
    { 
        xData[i] = (float *)malloc(COL * sizeof(float));
        if(xData[i]==NULL)
        {
            fprintf(stderr, "[Error] xData internal allocation failed.\n");
            exit(1);
        }
    }
    
    float **yData = (float **)malloc(ROW * sizeof(float *));
    if(yData==NULL)
    {
        fprintf(stderr, "[Error] yData allocation failed.\n");
        exit(1);        
    }  
    for (int i=0; i<ROW; i++)
    { 
        yData[i] = (float *)malloc(COL * sizeof(float));
        if(yData[i]==NULL)
        {
            fprintf(stderr, "[Error] yData internal allocation failed.\n");
            exit(1);
        }  
    }
    
    float **zData = (float **)malloc(ROW * sizeof(float *));
    if(zData==NULL)
    {
        fprintf(stderr, "[Error] zData allocation failed.\n");
        exit(1);        
    }  
    for (int i=0; i<ROW; i++)
    { 
        zData[i] = (float *)malloc(COL * sizeof(float));
        if(zData[i]==NULL)
        {
            fprintf(stderr, "[Error] zData internal allocation failed.\n");
            exit(1);
        }    
    }

    MPI_Bcast(&start, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&stop, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&step, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&M, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&N, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&tmax, 1, MPI_INT, 0, MPI_COMM_WORLD);

    readStart = MPI_Wtime();

    padding(M, xData, yData, zData);
    readData(start, stop, step, N, xData, yData, zData);

    readEnd = MPI_Wtime();

    chunk = tmax / (wSize * 2);
    remainder = tmax - (chunk * wSize * 2);
    lStart1 = (rank * chunk) + 1;
    lEnd1 = lStart1 + chunk - 1;

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
        fprintf(stdout, "%d, %e\n", dt, accumalate);
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
        fprintf(stdout, "%d, %e\n", dt, accumalate);
    }
    if(rank == 0)
        fprintf(stdout, "Read Time: %lf \n", (readEnd-readStart));

    MPI_Finalize();

    return 0;
}