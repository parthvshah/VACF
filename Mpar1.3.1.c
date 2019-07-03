#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>

void readData()
{
    
}

int main(int argc, char **argv)
{
    int rank, wSize;
    MPI_Status status;

    MPI_Init(&argc, &argv);

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &wSize);

    int start, stop, step, N, M, tmax;
   

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
                fprintf(stderr, "[Error] Command-line argument not recognized. \n");
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
    MPI_Bcast(&tmax, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&M, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&N, 1, MPI_INT, 0, MPI_COMM_WORLD);

    // Based on memory
    int batchParticles = 5;
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

    for (int batch = 0; batch < batches; batch++)
    {
        // Read into 3 * arrays based on batch and rank
        // Compute lCorr
    }

    for(int k = (M-1); k>=0; k--)
    {
        free(xData[k]);
        free(yData[k]);
        free(zData[k]);
    }

    free(xData);
    free(yData);
    free(zData);
    free(lCorr);
    free(gCorr);

    MPI_Reduce(lCorr, gCorr, tmax + 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

    if (rank == 0)
    {
        for (int k = 1; k <= tmax; k++)
            fprintf(stdout, "%d, %e\n", k, gCorr[k]);
        
    }

    MPI_Finalize();
    return 0;
}