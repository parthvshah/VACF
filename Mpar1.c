/*
VACF (Massively Parallel Correlation Decompositon)
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

void padding(int n, double **xData, double **yData, double **zData)
{
    for (int i = 0; i < n; i++)
    {
        xData[0][i] = yData[0][i] = zData[0][i] = 0.0;
    }
}

void readData(int lTimesteps, int start, int stop, int step, int N, int rank, double **xData, double **yData, double** zData)
{
    int timestep, particle, one, index, count;
    double xVel, yVel, zVel;
    FILE *fp;

    fp = fopen("./HISTORY_atoms/HISTORY_CLEAN_20", "r");
    if (fp == NULL)
        return;

    // Seek
    char line[256];

    short int skipLines = rank*(N-1);
    count = 0;
    while (count != skipLines)
    {
        fgets(line, sizeof line, fp);
        count++;
    }
    
    short int readLines = lTimesteps*(N-1);
    count = 0;
    while (fgets(line, sizeof line, fp) != NULL)
    {
        if(count == readLines)
            return;
        sscanf(line, "%d %d %d %lf %lf %lf", &timestep, &particle, &one, &xVel, &yVel, &zVel);
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
    int rank, wSize, lStart, lEnd, chunk, remainder;
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

    }

    MPI_Bcast(&start, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&stop, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&step, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&M, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&N, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&tmax, 1, MPI_INT, 0, MPI_COMM_WORLD);

    int tStart, tStop, lTimesteps;
    // Calculate using memory
    lTimesteps = (M-1)/wSize;
    tStart = (rank*step) + start;
    tStop = tStart + (lTimesteps*step);

    double **xData = (double**) malloc(sizeof(double*) * N);
    if(!xData)
    {
        free(xData);
        fprintf(stdout, "[Error - %d] xData not allocated.\n", EXIT_FAILURE);
        return EXIT_FAILURE;
    }
    for (int ii = 0; ii<N; ii++)
    {
        xData[ii] = (double*) malloc(sizeof(double) * lTimesteps);
        if(!xData[ii])
        {
            fprintf(stdout, "[Error - %d] Internal xData not allocated. \n", EXIT_FAILURE);
            for(int jj = ii; jj>=0; jj--)
                free(xData[jj]);
            
            return EXIT_FAILURE;
        }
    }

    
    double **yData = (double**) malloc(sizeof(double*) * N);
    if(!yData)
    {
        free(yData);
        fprintf(stdout, "[Error - %d] xData not allocated.\n", EXIT_FAILURE);
        return EXIT_FAILURE;
    }
    for (int ii = 0; ii<N; ii++)
    {
        yData[ii] = (double*) malloc(sizeof(double) * lTimesteps);
        if(!yData[ii])
        {
            fprintf(stdout, "[Error - %d] Internal xData not allocated. \n", EXIT_FAILURE);
            for(int jj = ii; jj>=0; jj--)
                free(yData[jj]);
                
            return EXIT_FAILURE;
        }
    }

    
    double **zData = (double**) malloc(sizeof(double*) * N);
    if(!zData)
    {
        free(zData);
        fprintf(stdout, "[Error - %d] xData not allocated.\n", EXIT_FAILURE);
        return EXIT_FAILURE;
    }
    for (int ii = 0; ii<N; ii++)
    {
        zData[ii] = (double*) malloc(sizeof(double) * lTimesteps);
        if(!zData[ii])
        {
            fprintf(stdout, "[Error - %d] Internal xData not allocated. \n", EXIT_FAILURE);
            for(int jj = ii; jj>=0; jj--)
                free(zData[jj]);
            
            return EXIT_FAILURE;
        }
    }
    
    readData(lTimesteps, tStart, tStop, step, N, rank, xData, yData, zData);

    if(rank==0)
    {
        printf("%d)\n", rank);

        for(int iii = 0; iii<N; iii++)
        {
            for(int jjj = 0; jjj<lTimesteps; jjj++)
            {
                printf("%e  ", xData[iii][jjj]);
            }
            printf("\n");
        }
    }

    // chunk = tmax / wSize;
    // remainder = tmax - (wSize * chunk);

    // lStart = (rank * chunk) + 1;
    // lEnd = lStart + chunk - 1;

    // // for (int dt = 1; dt <= tmax; dt++)
    // for (int dt = lStart; dt <= lEnd; dt++)
    // {
    //     count = 0;
    //     accumalate = 0.0;
    //     for (int t = 0; t < (M - dt); t++)
    //     {
    //         particle = 0.0;
    //         for (int i = 1; i < N; i++)
    //         {
    //             particle += xData[i][t] * xData[i][t + dt] +
    //                         yData[i][t] * yData[i][t + dt] +
    //                         zData[i][t] * zData[i][t + dt];
    //         }
    //         count++;
    //         accumalate += particle;
    //     }
    //     accumalate /= ((N - 1) * count);
    //     fprintf(stdout, "%d, %e\n", dt, accumalate);
    // }

    MPI_Finalize();

    return 0;
}