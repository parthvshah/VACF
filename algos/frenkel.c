#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#define ROW 6913
#define COL 200001

void readData(int start, int stop, int step, int N, float **xData, float **yData, float **zData, char *fileName)
{
    int timestep, particle, one, index, count, maxCount;
    float xVel, yVel, zVel;
    FILE *fp;

    fp = fopen(fileName, "r");
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
/*
    For particle = 1, n
        index <- 0
        For timestep = 1, total timesteps
            accumulator = x, y, z velocities at timestep
            If timestep % blockSize == 0
                xData, yData, zData at index <- accumulator / blockSize
                reset accumulator
                index++
*/
int blockAverage(float **xData, float **yData, float **zData, int N, int M, int blockSize)
{
    int index = 0;

    for(int i = 1; i < N; i++)
    {
        index = 0;
        float xAccumulator = 0.0, yAccumulator = 0.0, zAccumulator = 0.0;
        for(int j = 0; j < M; j++)
        {
            xAccumulator += xData[i][j];
            yAccumulator += yData[i][j];
            zAccumulator += zData[i][j];

            if((j+1) % blockSize == 0)
            {
                xData[i][index] = xAccumulator / blockSize;
                yData[i][index] = yAccumulator / blockSize;
                zData[i][index] = zAccumulator / blockSize;

                xAccumulator = 0.0; yAccumulator = 0.0; zAccumulator = 0.0;
                index++;
            }
        }
    }

    return index+1;
}

int main(int argc, char **argv)
{
    int start;
    int stop;
    int step;
    int tmax = 0;

    int blockSize;
    int level = 1;

    int N, M;
    
    char fileName[128];
    for (int c = 1; c < argc; c++)
    {
        if (!strcmp(argv[c], "-p"))
            sscanf(argv[++c], "%d,%d,%d", &start, &stop, &step);
        else if (!strcmp(argv[c], "-a"))
            N = atoi(argv[++c]);
        else if (!strcmp(argv[c], "-i"))
            tmax = atoi(argv[++c]);
        else if (!strcmp(argv[c], "-b"))
            blockSize = atoi(argv[++c]);
        else if (!strcmp(argv[c], "-l"))
            level = atoi(argv[++c]);
        else if (!strcmp(argv[c], "-f"))
            sscanf(argv[++c], "%s", fileName);
        else
        {
            fprintf(stderr, "[Error] Command-line argument not recognized.\n");
            exit(-1);
        }
    }

    float **xData = NULL;
    xData = (float**) malloc(sizeof(float*) * ROW);
    if(!xData)
    {
        free(xData);
        fprintf(stdout, "[Error] xData not allocated.\n");
        return EXIT_FAILURE;
    }
    for (int ii = 0; ii<ROW; ii++)
    {
        xData[ii] = NULL;
        xData[ii] = (float*) malloc(sizeof(float) * COL);
        if(!xData[ii])
        {
            fprintf(stdout, "[Error] Internal xData not allocated. \n");
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
        fprintf(stdout, "[Error] xData not allocated.\n");
        return EXIT_FAILURE;
    }
    for (int ii = 0; ii<ROW; ii++)
    {
        yData[ii] = NULL;
        yData[ii] = (float*) malloc(sizeof(float) * COL);
        if(!yData[ii])
        {
            fprintf(stdout, "[Error] Internal xData not allocated. \n");
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
        fprintf(stdout, "[Error] xData not allocated.\n");
        return EXIT_FAILURE;
    }
    for (int ii = 0; ii<ROW; ii++)
    {
        zData[ii] = NULL;
        zData[ii] = (float*) malloc(sizeof(float) * COL);
        if(!zData[ii])
        {
            fprintf(stdout, "[Error] Internal xData not allocated. \n");
            for(int jj = ii; jj>=0; jj--)
                free(zData[jj]);
            
            return EXIT_FAILURE;
        }
    }

    N++;                             // Number of particles
    M = ((stop - start) / step) + 1; // Number of timesteps

    tmax = (tmax == 0) ? (M / 3) : tmax;

    readData(start, stop, step, N, xData, yData, zData, fileName);

    while(level != 0)
    {
        M = blockAverage(xData, yData, zData, N, M, blockSize);
        level--;
    }

    double accumalate, particle;
    int count;
    for (int dt = 1; dt <= tmax; dt++)
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
        fprintf(stdout, "%d, %e\n", dt*blockSize/2, accumalate);
    }

    return 0;
}