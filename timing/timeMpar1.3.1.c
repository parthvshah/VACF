#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>

void readData(int row, int col, int batch, int rank, int wSize, double **xData, double **yData, double **zData, int start, int stop, int step, int batchParticles)
{
    int startParticle, endParticle;
    startParticle = (batch * (batchParticles * wSize)) + (batchParticles * rank) + 1;
    endParticle = startParticle + batchParticles - 1;

    FILE *fp;

    char **fileNames = NULL;
    fileNames = (char**) malloc(sizeof(char*) * batchParticles);
    if(!fileNames)
    {
        fprintf(stderr, "[Error - %d] fileNames not allocated. \n", rank);
        free(fileNames);
        return;
    }
    for(int z = 0; z<batchParticles; z++)
    {
        fileNames[z] = NULL;
        fileNames[z] = (char*) malloc(sizeof(char) * 4);
        if(!fileNames[z])
        {
            fprintf(stderr, "[Error - %d] Internal fileNames not allocated. \n", rank);
            for(int zz = z; zz>=0; zz--)
                free(fileNames[zz]);

            return;
        }
    }

    char buffer[7];
    for(int i = 0; i < batchParticles; i++)
    {
        snprintf(buffer, 7, "./%04d", (startParticle+i));
        strcpy(fileNames[i], buffer);
    }

    int timestep, rParticle, one, index;
    double xVel, yVel, zVel;

    for (int p = 0; p < batchParticles; p++)
    {
        fp = fopen(fileNames[p], "r");
        if (fp == NULL)
        {
            fprintf(stderr, "[Error - %d] %s: File open error. \n", rank, fileNames[p]);
            return;
        }


        int count = 0;
        int readLines = (stop - start) / step;

        char line[256];
        while (fgets(line, sizeof line, fp) != NULL)
        {
            if(count == readLines)
                break;
            sscanf(line, "%d %d %d %lf %lf %lf", &timestep, &rParticle, &one, &xVel, &yVel, &zVel);
            index = (timestep - start) / step;

            xData[p][index] = xVel;
            yData[p][index] = yVel;
            zData[p][index] = zVel;
            count++;
        }
        fclose(fp);
    }

    for(int j = 0; j<batchParticles; j++)
    {
        free(fileNames[j]);
    }
    free(fileNames);
}

double clockTime(double startT)
{
    double endT;
    endT = MPI_Wtime() - startT;
    return endT;
}

int main(int argc, char **argv)
{
    int rank, wSize;
    MPI_Status status;

    MPI_Init(&argc, &argv);

    double mpi_setup = MPI_Wtime();

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &wSize);

    double *lTime = NULL;
    lTime = (double *)malloc(sizeof(double) * 5);
    if(!lTime)
    {
        free(lTime);
        fprintf(stdout, "[Error - %d] lTime not allocated.\n", rank);
        return EXIT_FAILURE;
    }

    double *gTime = NULL;
    gTime = (double *)malloc(sizeof(double) * 5);
    if(!gTime)
    {
        free(gTime);
        fprintf(stdout, "[Error - %d] gTime not allocated.\n", rank);
        return EXIT_FAILURE;
    }

    lTime[0] = clockTime(mpi_setup);

    int start, stop, step, N, M, tmax, batchParticles;
   

    if (rank == 0)
    {
        for (int c = 1; c < argc; c++)
        {
            if (!strcmp(argv[c], "-p"))
                sscanf(argv[++c], "%d,%d,%d", &start, &stop, &step);
            else if (!strcmp(argv[c], "-a"))
                N = atoi(argv[++c]);
            else if (!strcmp(argv[c], "-i"))
                tmax = atoi(argv[++c]);
            else if (!strcmp(argv[c], "-bp"))
                batchParticles = atoi(argv[++c]);
            else
            {
                fprintf(stderr, "[Error] Command-line argument not recognized. \n");
                exit(-1);
            }
        }

        N++;                             // Number of particles
        M = ((stop - start) / step) + 1; // Number of timesteps

        tmax = (tmax == 0) ? (M / 3) : tmax;
    }

    double bcast = MPI_Wtime();

    MPI_Bcast(&batchParticles, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&start, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&stop, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&step, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&tmax, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&M, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&N, 1, MPI_INT, 0, MPI_COMM_WORLD);

    lTime[2] = clockTime(bcast);

    int batches = (N-1) / (wSize * batchParticles);    

    double *lCorr = NULL;
    lCorr = (double*) malloc(sizeof(double) * (tmax + 1));
    if(!lCorr)
    {
        fprintf(stderr, "[Error - %d] lCorr not allocated.\n", rank);
        free(lCorr);
        return EXIT_FAILURE;
    }

    double *gCorr = NULL;
    gCorr = (double*) malloc(sizeof(double) * (tmax + 1));
    if(!gCorr)
    {
        fprintf(stderr, "[Error - %d] gCorr not allocated.\n", rank);
        free(gCorr);
        return EXIT_FAILURE;
    }

    double mem_alloc = MPI_Wtime();

    double **xData = NULL;
    xData = (double**) malloc(sizeof(double*) * batchParticles);
    if(!xData)
    {
        fprintf(stderr, "[Error - %d] xData not allocated.\n", rank);
        free(xData);
        return EXIT_FAILURE;
    }
    for (int ii = 0; ii<batchParticles; ii++)
    {
        xData[ii] = NULL;
        xData[ii] = (double*) malloc(sizeof(double) * M);
        if(!xData[ii])
        {
            fprintf(stderr, "[Error - %d] Internal xData not allocated. \n", rank);
            for(int jj = ii; jj>=0; jj--)
                free(xData[jj]);

            return EXIT_FAILURE;
        }
    }

    double **yData = NULL;
    yData = (double**) malloc(sizeof(double*) * batchParticles);
    if(!yData)
    {
        fprintf(stderr, "[Error - %d] yData not allocated.\n", rank);
        free(yData);
        return EXIT_FAILURE;
    }
    for (int ii = 0; ii<batchParticles; ii++)
    {
        yData[ii] = NULL;
        yData[ii] = (double*) malloc(sizeof(double) * M);
        if(!yData[ii])
        {
            fprintf(stderr, "[Error - %d] Internal yData not allocated. \n", rank);
            for(int jj = ii; jj>=0; jj--)
                free(yData[jj]);

            return EXIT_FAILURE;
        }
    }

    double **zData = NULL;
    zData = (double**) malloc(sizeof(double*) * batchParticles);
    if(!zData)
    {
        fprintf(stderr, "[Error - %d] zData not allocated.\n", rank);
        free(zData);
        return EXIT_FAILURE;
    }
    for (int ii = 0; ii<batchParticles; ii++)
    {
        zData[ii] = NULL;
        zData[ii] = (double*) malloc(sizeof(double) * M);
        if(!zData[ii])
        {
            fprintf(stderr, "[Error - %d] Internal zData not allocated. \n", rank);
            for(int jj = ii; jj>=0; jj--)
                free(zData[jj]);

            return EXIT_FAILURE;
        }
    }

    lTime[1] = clockTime(mem_alloc);

    unsigned long int count;
    double accumalate, particle;
    
    // Read time is total, for all batches. To calculate per batch, roughly divide by no. of batches
    double read_d = MPI_Wtime();

    for (int batch = 0; batch < batches; batch++)
    {
        // Read into 3 * arrays based on batch and rank
        readData(batchParticles, M, batch, rank, wSize, xData, yData, zData, start, stop, step, batchParticles);

        // Compute lCorr

        // for (int dt = 1; dt <= tmax; dt++)
        // {
        //     count = 0;
        //     accumalate = 0.0;
        //     for (int t = 0; t < (M - dt); t++)
        //     {
        //         particle = 0.0;
        //         // for (int i = 1; i < N; i++)
        //         for (int i = 0; i < batchParticles; i++)
        //         {
        //             particle += xData[i][t] * xData[i][t + dt] +
        //                         yData[i][t] * yData[i][t + dt] +
        //                         zData[i][t] * zData[i][t + dt];
        //         }
        //         count++;
        //         accumalate += particle;
        //     }
        //     accumalate /= ((N - 1) * count);
        //     lCorr[dt] += accumalate;
        //     // fprintf(stdout, "%d, %e\n", dt, accumalate);
        // }
        MPI_Barrier(MPI_COMM_WORLD);
    }

    lTime[3] = clockTime(read_d);


    double collect = MPI_Wtime();
    MPI_Reduce(lCorr, gCorr, tmax + 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    lTime[4] = clockTime(collect);

    MPI_Reduce(lTime, gTime, 5, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);


    if (rank == 0)
    {
        for (int k = 0; k < 5; k++)
        {
            fprintf(stdout, "%d. %lf\n", k+1, gTime[k]/wSize);
        }

        // for(int k = (M-1); k>=0; k--)
        // {
        //     free(xData[k]);
        //     free(yData[k]);
        //     free(zData[k]);
        // }

        // free(xData);
        // free(yData);
        // free(zData);
        // free(lCorr);
        // free(gCorr);  
    }

    MPI_Finalize();
    return 0;
}