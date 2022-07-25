#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#define ROW 6913
#define COL 200001

int readData(int start, int stop, int step, int N, float **xData, float **yData, float **zData, int *indexData, int m, char *fileName)
{
    int timestep, particle, one, index, count, maxCount;
    float xVel, yVel, zVel;
    FILE *fp;
    
    fp = fopen(fileName, "r");
    if (fp == NULL)
        return 0;

    char line[256];

    N--;
    count = 0;
    maxCount = 100 * m * N;

    while (count <= maxCount)
    {
        fgets(line, sizeof line, fp);
        sscanf(line, "%d %d %d %f %f %f", &timestep, &particle, &one, &xVel, &yVel, &zVel);
        index = (timestep - start) / step;

        xData[particle][index] = xVel;
        yData[particle][index] = yVel;
        zData[particle][index] = zVel;

        count++;
    }

    int skip, i, j;
    skip = m;

    do
    {
        count = 0;

        while(count <= 100 * N)
        {
            for(int n = 0; n < N; n++)
            {
                fgets(line, sizeof line, fp);
                count++;
                sscanf(line, "%d %d %d %f %f %f", &timestep, &particle, &one, &xVel, &yVel, &zVel);
                xData[particle][index] = xVel;
                yData[particle][index] = yVel;
                zData[particle][index] = zVel;
            }
            indexData[index] = timestep;
            index++;

            j = 0;
            while(j < (skip - 1) * N)
            {
                fgets(line, sizeof line, fp);
                j++;
            }
        }

        skip *= m;

    } while(!feof(fp));

    fclose(fp);
    return index;
}

int main(int argc, char **argv)
{
    int start;
    int stop;
    int step;
    int tmax = 0;

    int blockSize;
    int level = 1;

    int N, M, m;
    
    char fileName[128];

    for (int c = 1; c < argc; c++)
    {
        if (!strcmp(argv[c], "-p"))
            sscanf(argv[++c], "%d,%d,%d", &start, &stop, &step);
        else if (!strcmp(argv[c], "-a"))
            N = atoi(argv[++c]);
        else if (!strcmp(argv[c], "-i"))
            tmax = atoi(argv[++c]);
        else if (!strcmp(argv[c], "-m"))
            m = atoi(argv[++c]);
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

    int *indexData = (int*) malloc(sizeof(int) * COL);

    N++;                             // Number of particles
    M = ((stop - start) / step) + 1; // Number of timesteps

    tmax = (tmax == 0) ? (M / 3) : tmax;

    int size = readData(start, stop, step, N, xData, yData, zData, indexData, m, fileName);

    double accumalate, particle;
    int count;
    for (int dt = 1; dt <= tmax; dt++)
    {
        count = 0;
        accumalate = 0.0;
        for (int t = 0; t < (size - dt); t++)
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

    return 0;
}