#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>

void readData(int batch, int rank, int wSize, int N, int batchTimesteps, int tmax, double **xData, double **yData, double **zData, int start, int stop, int step)
{
    int startTimestep, endTimestep;
    startTimestep = start + (batch * wSize * step) + (batchTimesteps * rank * step);
    endTimestep = startTimestep + ((batchTimesteps + tmax - 1) * step);

    int skipLines, readLines;
    skipLines = (startTimestep - start) / step * (N-1);
    readLines = (batchTimesteps + tmax) * (N-1);

    int timestep, particle, one, index, count, maxCount;
    double xVel, yVel, zVel;
    FILE *fp;
    // printf("%d)%d) %d %d\n", batch, rank, startTimestep, endTimestep);


    fp = fopen("./HISTORY_atoms/HISTORY_CLEAN_25", "r");
    if (fp == NULL)
        return;
    
    count = 0;
    char line[256];
    while (count != skipLines)
    {
        fgets(line, sizeof line, fp);
        count++;
    }

    count = 0;
    while (fgets(line, sizeof line, fp) != NULL)
    {
        if(count == readLines)
            break;
        sscanf(line, "%d %d %d %lf %lf %lf", &timestep, &particle, &one, &xVel, &yVel, &zVel);
        index = (timestep - startTimestep) / step;
        xData[particle][index] = xVel;
        yData[particle][index] = yVel;
        zData[particle][index] = zVel;

        count++;
    }
    fclose(fp);

} 

int main(int argc, char **argv)
{
    int rank, wSize;
    MPI_Status status;

    MPI_Init(&argc, &argv);

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &wSize);

    int start, stop, step, N, M, tmax;
    int batchTimesteps = 1;    

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
            else if (!strcmp(argv[c], "-bt"))
                batchTimesteps = atoi(argv[++c]);
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

    MPI_Bcast(&batchTimesteps, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&start, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&stop, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&step, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&tmax, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&M, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&N, 1, MPI_INT, 0, MPI_COMM_WORLD);

    // int batchTimesteps = 1;
    int batches = (M-1) / (wSize * batchTimesteps);
    if(batches == 0 && rank == 0)
    {
        fprintf(stderr, "[Error] Incorrect number of batches. \n");
        return 1;
    }

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

    double **xData = NULL;
    xData = (double**) malloc(sizeof(double*) * N);
    if(!xData)
    {
        fprintf(stderr, "[Error - %d] xData not allocated.\n", rank);
        free(xData);
        return EXIT_FAILURE;
    }
    for (int ii = 0; ii<N; ii++)
    {
        xData[ii] = NULL;
        xData[ii] = (double*) malloc(sizeof(double) * (batchTimesteps+tmax));
        if(!xData[ii])
        {
            fprintf(stderr, "[Error - %d] Internal xData not allocated. \n", rank);
            for(int jj = ii; jj>=0; jj--)
                free(xData[jj]);

            return EXIT_FAILURE;
        }
    }

    double **yData = NULL;
    yData = (double**) malloc(sizeof(double*) * N);
    if(!yData)
    {
        fprintf(stderr, "[Error - %d] yData not allocated.\n", rank);
        free(yData);
        return EXIT_FAILURE;
    }
    for (int ii = 0; ii<N; ii++)
    {
        yData[ii] = NULL;
        yData[ii] = (double*) malloc(sizeof(double) * (batchTimesteps+tmax));
        if(!yData[ii])
        {
            fprintf(stderr, "[Error - %d] Internal yData not allocated. \n", rank);
            for(int jj = ii; jj>=0; jj--)
                free(yData[jj]);

            return EXIT_FAILURE;
        }
    }

    double **zData = NULL;
    zData = (double**) malloc(sizeof(double*) * N);
    if(!zData)
    {
        fprintf(stderr, "[Error - %d] zData not allocated.\n", rank);
        free(zData);
        return EXIT_FAILURE;
    }
    for (int ii = 0; ii<N; ii++)
    {
        zData[ii] = NULL;
        zData[ii] = (double*) malloc(sizeof(double) * (batchTimesteps+tmax));
        if(!zData[ii])
        {
            fprintf(stderr, "[Error - %d] Internal zData not allocated. \n", rank);
            for(int jj = ii; jj>=0; jj--)
                free(zData[jj]);

            return EXIT_FAILURE;
        }
    }

    int count, skipTimes = -1;
    double accumalate, particle;

    
    for (int batch = 0; batch < batches; batch++)
    {
        // Read into 3 * arrays based on batch and rank
        readData(batch, rank, wSize, N, batchTimesteps, tmax, xData, yData, zData, start, stop, step);
    
        if(rank != 0)
            skipTimes = tmax - 1;
        // Compute lCorr
        for (int dt = 1; dt <= tmax; dt++)
        {
            count = 0;
            accumalate = 0.0;
            for (int t = 0; t < (batchTimesteps + tmax - dt); t++)
            {
                particle = 0.0;
                for (int i = 1; i < N; i++)
                {
                    particle += xData[i][t] * xData[i][t + dt] +
                                yData[i][t] * yData[i][t + dt] +
                                zData[i][t] * zData[i][t + dt];
                }
                if(count >= skipTimes)
                {
                    accumalate += particle;
                }
                
                count++;
                
            }
            // accumalate /= ((N - 1) * count);
            lCorr[dt] += accumalate;
            skipTimes--;
            // fprintf(stdout, "%d, %e\n", dt, accumalate);
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }

    MPI_Reduce(lCorr, gCorr, tmax + 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

    if (rank == 0)
    {
        for (int l = 1; l <= tmax; l++)
        {
            fprintf(stdout, "%d, %e\n", l, gCorr[l] / ((N-1) * (M-l)));
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