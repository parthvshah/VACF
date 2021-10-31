/*
Calculate Velocity Auto-Correlation (Parallel - OpenMPI)
This program calculates the VACF from a given HISTORY file. Scripts to clean the data,
plot the graph and calculate the coefficient of diffusion are external and in python.

Command line arguments:
    -p START,STOP,STEP
    -a PARTICLES
    -i ITERATIONS

Date: 10-6-2019
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <omp.h>

#define ROW 1024
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

void readData(int start, int step)
{
    int timestep, particle, index;
    double xVel, yVel, zVel;
    FILE *fp;

    fp = fopen("HISTORY_CLEAN", "r");
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
    int start;
    int stop;
    int step;
    int tmax = 0;
    int rank, wSize, chunk, remainder, lStart, lEnd;

    int N, M;

    for (int c = 1; c < argc; c++)
    {
        if (!strcmp(argv[c], "-p"))
            sscanf(argv[++c], "%d,%d,%d", &start, &stop, &step);
        else if (!strcmp(argv[c], "-a"))
            N = atoi(argv[++c]);
        else if (!strcmp(argv[c], "-i"))
            tmax = atoi(argv[++c]);
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
    readData(start, step);

    fprintf(stdout, "timestep, vacf\n");

    omp_set_num_threads(4);
    wSize = omp_get_max_threads();

    chunk = tmax / wSize;
    remainder = tmax - (wSize * chunk);

    double accumalate, particle;
    int count;
#pragma omp parallel
    {
        rank = omp_get_thread_num();
        lStart = (rank * chunk) + 1;
        lEnd = lStart - 1 + chunk;
        if (rank == (wSize - 1))
        {
            lEnd += remainder;
        }
        // for (int dt = 1; dt <= tmax; dt++)
        for (int dt = lStart; dt <= lEnd; dt++)
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
    }

    return 0;
}
