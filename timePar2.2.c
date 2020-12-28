#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>

int readData(int wRank, int batch, int corr, int N, int start, int stop, int step, float **xData, float **yData, float **zData, int batchTimesteps)
{
    int timestep, particle, one, count, skip, maxCount, index, firstTimestep;
    float xVel, yVel, zVel;
    FILE *fp;
    char line[256];

    fp = fopen("./HISTORY_atoms/HISTORY_CLEAN_500_l", "r");
    if (fp == NULL)
        return 0;

    count = 1;
    skip = batch * batchTimesteps * (N - 1);

    if (skip != 0)
    {
        while (fgets(line, sizeof line, fp) != NULL)
        {
            if (count == skip)
                break;
            count++;
        }
    }

    count = 1;
    maxCount = (batchTimesteps + corr) * (N - 1);

    while (fgets(line, sizeof line, fp) != NULL)
    {
        if (count == maxCount)
            break;
        sscanf(line, "%d %d %d %f %f %f", &timestep, &particle, &one, &xVel, &yVel, &zVel);

        if (count == 1)
            firstTimestep = timestep;

        index = (timestep - firstTimestep) / step;
        xData[particle][index] = xVel;
        yData[particle][index] = yVel;
        zData[particle][index] = zVel;

        count++;
    }
    fclose(fp);
    return index;
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

    double total = MPI_Wtime();
    double mpi_setup = MPI_Wtime();

    MPI_Comm_rank(MPI_COMM_WORLD, &wRank);
    MPI_Comm_size(MPI_COMM_WORLD, &wSize);

    int colour = wRank / 12;
    MPI_Comm MPI_COMM_ROW, MPI_COMM_ROWROOT;
    MPI_Comm_split(MPI_COMM_WORLD, colour, wRank, &MPI_COMM_ROW);

    int rRank, rSize;
    MPI_Comm_rank(MPI_COMM_ROW, &rRank);
    MPI_Comm_size(MPI_COMM_ROW, &rSize);

    int *ranks = (int *)malloc(sizeof(int) * wSize / 12);
    int iCount = 0;
    for (int iRank = 0; iRank < wSize; iRank++)
        if (iRank % 12 == 0)
            ranks[iCount++] = iRank;
    MPI_Group MPI_GROUP_WORLD, MPI_GROUP_ROWROOT;
    MPI_Comm_group(MPI_COMM_WORLD, &MPI_GROUP_WORLD);
    MPI_Group_incl(MPI_GROUP_WORLD, (wSize / 12), ranks, &MPI_GROUP_ROWROOT);
    MPI_Comm_create(MPI_COMM_WORLD, MPI_GROUP_ROWROOT, &MPI_COMM_ROWROOT);

    double *lTime = NULL;
    lTime = (double *)malloc(sizeof(double) * 6);
    if (!lTime)
    {
        free(lTime);
        fprintf(stdout, "[Error - %d] lTime not allocated.\n", wRank);
        return EXIT_FAILURE;
    }

    double *gTime = NULL;
    gTime = (double *)malloc(sizeof(double) * 6);
    if (!gTime)
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
    int batchTimesteps, batches;
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
            else if (!strcmp(argv[c], "-bt"))
                batchTimesteps = atoi(argv[++c]);
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

        batches = (M / batchTimesteps) + 1;

        // padding(M);
        // readData(start, step, fileName);
    }

    double bcast = MPI_Wtime();

    MPI_Bcast(&batches, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&batchTimesteps, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&start, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&stop, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&step, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&M, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&N, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&tmax, 1, MPI_INT, 0, MPI_COMM_WORLD);

    lTime[2] = clockTime(bcast);

    double mem_alloc = MPI_Wtime();

    float **xData = NULL;
    xData = (float **)malloc(sizeof(float *) * N);
    if (!xData)
    {
        free(xData);
        fprintf(stdout, "[Error - %d] xData not allocated.\n", wRank);
        return EXIT_FAILURE;
    }
    for (int ii = 0; ii < N; ii++)
    {
        xData[ii] = NULL;
        xData[ii] = (float *)malloc(sizeof(float) * (batchTimesteps + tmax + 1));
        if (!xData[ii])
        {
            fprintf(stdout, "[Error - %d] Internal xData not allocated. \n", wRank);
            for (int jj = ii; jj >= 0; jj--)
                free(xData[jj]);

            return EXIT_FAILURE;
        }
    }

    float **yData = NULL;
    yData = (float **)malloc(sizeof(float *) * N);
    if (!yData)
    {
        free(yData);
        fprintf(stdout, "[Error - %d] xData not allocated.\n", wRank);
        return EXIT_FAILURE;
    }
    for (int ii = 0; ii < N; ii++)
    {
        yData[ii] = NULL;
        yData[ii] = (float *)malloc(sizeof(float) * (batchTimesteps + tmax + 1));
        if (!yData[ii])
        {
            fprintf(stdout, "[Error - %d] Internal xData not allocated. \n", wRank);
            for (int jj = ii; jj >= 0; jj--)
                free(yData[jj]);

            return EXIT_FAILURE;
        }
    }

    float **zData = NULL;
    zData = (float **)malloc(sizeof(float *) * N);
    if (!zData)
    {
        free(zData);
        fprintf(stdout, "[Error - %d] xData not allocated.\n", wRank);
        return EXIT_FAILURE;
    }
    for (int ii = 0; ii < N; ii++)
    {
        zData[ii] = NULL;
        zData[ii] = (float *)malloc(sizeof(float) * (batchTimesteps + tmax + 1));
        if (!zData[ii])
        {
            fprintf(stdout, "[Error - %d] Internal xData not allocated. \n", wRank);
            for (int jj = ii; jj >= 0; jj--)
                free(zData[jj]);

            return EXIT_FAILURE;
        }
    }

    float *lCorr = NULL;
    lCorr = (float *)malloc(sizeof(float) * (tmax + 1));
    if (!lCorr)
    {
        free(lCorr);
        fprintf(stdout, "[Error - %d] lCorr not allocated.\n", wRank);
        return EXIT_FAILURE;
    }

    float *gCorr = NULL;
    gCorr = (float *)malloc(sizeof(float) * (tmax + 1));
    if (!gCorr)
    {
        free(gCorr);
        fprintf(stdout, "[Error - %d] gCorr not allocated.\n", wRank);
        return EXIT_FAILURE;
    }

    lTime[1] = clockTime(mem_alloc);

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

    int timestepsRead;
    for (int batch = 0; batch < batches; batch++)
    {
        double read_d = MPI_Wtime();
        if(batch == 0)
            lTime[3] = 0.0;
        timestepsRead = readData(wRank, batch, tmax, N, start, stop, step, xData, yData, zData, batchTimesteps);
        lTime[3] += clockTime(read_d);

        if (timestepsRead == 0)
        {
            fprintf(stdout, "[Error - %d] No timesteps read.\n", wRank);
            return EXIT_FAILURE;
        }

        for (int dt = lStart; dt <= lEnd; dt++)
        {
            count = 0;
            accumalate = 0.0;

            // for (int t = 0; t < (M - dt); t++)
            for (int t = 0; t < (timestepsRead - dt); t++)
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
            // MPI_Barrier(MPI_COMM_WORLD);

            accumalate /= ((N - 1) * count);
            lCorr[dt] = accumalate;
            // fprintf(stdout, "%d, %e\n", dt, accumalate);
        }
    }

    double collect = MPI_Wtime();

    MPI_Reduce(lCorr, gCorr, tmax + 1, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_ROW);

    lTime[4] = clockTime(collect);
    lTime[5] = clockTime(total);

    MPI_Reduce(lTime, gTime, 5, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

    if (rRank == 0)
    {
        for (int k = lStart; k <= lEnd; k++)
        {
            fprintf(stdout, "%d, %e\n", k, gCorr[k]);
        }
    }

    if (wRank == 0)
    {
        for (int k = 0; k < 6; k++)
        {
            fprintf(stdout, "%d. %lf\n", k + 1, gTime[k] / wSize);
        }
    }

    MPI_Finalize();
    return 0;
}