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

double clockTime(double startT)
{
    double endT;
    endT = MPI_Wtime() - startT;
    return endT;
}

int main(int argc, char **argv)
{

    int wRank, wSize, lChunk, NChunk, lRemainder, NRemainder, NStart, NEnd, lStart, lEnd;
    MPI_Status status;

    MPI_Init(&argc, &argv);

    double mpi_setup = MPI_Wtime();

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

    double *lTime = NULL;
    lTime = (double *)malloc(sizeof(double) * 5);
    if(!lTime)
    {
        free(lTime);
        fprintf(stdout, "[Error - %d] lTime not allocated.\n", wRank);
        return EXIT_FAILURE;
    }

    double *gTime = NULL;
    gTime = (double *)malloc(sizeof(double) * 5);
    if(!gTime)
    {
        free(gTime);
        fprintf(stdout, "[Error - %d] gTime not allocated.\n", wRank);
        return EXIT_FAILURE;
    }

    lTime[0] = clockTime(mpi_setup);

    int start;
    int stop;
    int step;
    int N, M;
    int tmax = 0;

    double accumalate, particle;
    int count;

    if (wRank == 0)
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

    double mem_alloc = MPI_Wtime();

    float **xData = NULL;
    xData = (float**) malloc(sizeof(float*) * ROW);
    if(!xData)
    {
        free(xData);
        fprintf(stdout, "[Error - %d] xData not allocated.\n", wRank);
        return EXIT_FAILURE;
    }
    for (int ii = 0; ii<ROW; ii++)
    {
        xData[ii] = NULL;
        xData[ii] = (float*) malloc(sizeof(float) * COL);
        if(!xData[ii])
        {
            fprintf(stdout, "[Error - %d] Internal xData not allocated. \n", wRank);
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
        fprintf(stdout, "[Error - %d] xData not allocated.\n", wRank);
        return EXIT_FAILURE;
    }
    for (int ii = 0; ii<ROW; ii++)
    {
        yData[ii] = NULL;
        yData[ii] = (float*) malloc(sizeof(float) * COL);
        if(!yData[ii])
        {
            fprintf(stdout, "[Error - %d] Internal xData not allocated. \n", wRank);
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
        fprintf(stdout, "[Error - %d] xData not allocated.\n", wRank);
        return EXIT_FAILURE;
    }
    for (int ii = 0; ii<ROW; ii++)
    {
        zData[ii] = NULL;
        zData[ii] = (float*) malloc(sizeof(float) * COL);
        if(!zData[ii])
        {
            fprintf(stdout, "[Error - %d] Internal xData not allocated. \n", wRank);
            for(int jj = ii; jj>=0; jj--)
                free(zData[jj]);
            
            return EXIT_FAILURE;
        }
    }

    lTime[1] = clockTime(mem_alloc);

    double bcast = MPI_Wtime();

    MPI_Bcast(&start, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&stop, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&step, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&M, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&N, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&tmax, 1, MPI_INT, 0, MPI_COMM_WORLD);

    lTime[2] = clockTime(bcast);

    float *lCorr = NULL;
    lCorr = (float *)malloc(sizeof(float) * (tmax + 1));
    if(!lCorr)
    {
        free(lCorr);
        fprintf(stdout, "[Error - %d] lCorr not allocated.\n", wRank);
        return EXIT_FAILURE;
    }

    float *gCorr = NULL;
    gCorr = (float *)malloc(sizeof(float) * (tmax + 1));
    if(!gCorr)
    {
        free(gCorr);
        fprintf(stdout, "[Error - %d] gCorr not allocated.\n", wRank);
        return EXIT_FAILURE;
    }

    padding(M, xData, yData, zData);
    double read_d = MPI_Wtime();
    readData(start, stop, step, N, xData, yData, zData);
    lTime[3] = clockTime(read_d);

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

    // for (int dt = lStart; dt <= lEnd; dt++)
    // {
    //     count = 0;
    //     accumalate = 0.0;
    //     for (int t = 0; t < (M - dt); t++)
    //     {
    //         particle = 0.0;
    //         for (int i = NStart; i <= NEnd; i++)
    //         {
    //             particle += xData[i][t] * xData[i][t + dt] +
    //                         yData[i][t] * yData[i][t + dt] +
    //                         zData[i][t] * zData[i][t + dt];
    //         }
    //         count++;
    //         accumalate += particle;
    //     }
    //     accumalate /= ((N - 1) * count);
    //     lCorr[dt] = accumalate;
    //     // fprintf(stdout, "%d, %e\n", dt, accumalate);
    // }

    double collect = MPI_Wtime();

    MPI_Reduce(lCorr, gCorr, tmax + 1, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_ROW);

    lTime[4] = clockTime(collect);

    MPI_Reduce(lTime, gTime, 5, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

    if (wRank == 0)
    {
        for (int k = 0; k < 5; k++)
        {
            fprintf(stdout, "%d. %e\n", k+1, gTime[k]/wSize);
        }
    }
    
    // if (rRank == 0)
    // {
    //     for (int k = lStart; k <= lEnd; k++)
    //     {
    //         fprintf(stdout, "%d, %e\n", k, gCorr[k]);
    //     }
    // }

    MPI_Finalize();
    return 0;
}