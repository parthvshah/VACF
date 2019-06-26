/*
Calculate Velocity Auto-Correlation (Serial)
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

Date: 6-6-2019
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define ROW 2064
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

void readData(int start, int step, char *fileName)
{
    int timestep, particle, index;
    double xVel, yVel, zVel;
    FILE *fp;

    fp = fopen(fileName, "r");
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

    int N, M;
    
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

    padding(M);
    readData(start, step, fileName);

    fprintf(stdout, "timestep, vacf\n");

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
        fprintf(stdout, "%d, %e\n", dt, accumalate);
    }

    return 0;
}