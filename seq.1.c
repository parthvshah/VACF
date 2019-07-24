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
#include "hiredis/hiredis.h"

void padding(int n, float **xData, float **yData, float **zData)
{
    for (int i = 0; i < n; i++)
    {
        xData[0][i] = yData[0][i] = zData[0][i] = 0.0;
    }
}

float atofc(char *in)
{
    float out;
    if(in != NULL)
    {
        out = atof(in);
    }
    else
    {
        out = 0.0;
    }
    return out;
}

float retreive(int i, int t, int dt, redisContext *c)
{
    // particle += retreive(0, i, t, c) * retreive(0, i, (t+dt), c) + 
    //             retreive(1, i, t, c) * retreive(1, i, (t+dt), c) + 
    //             retreive(2, i, t, c) * retreive(2, i, (t+dt), c);
    redisReply *reply;
    float value = 0.0;
    char buffer1[256];
    char buffer2[256];
    sprintf(buffer1, "%d:%d", i, t);        
    sprintf(buffer2, "%d:%d", i, (t+dt));        

    reply = redisCommand(c, "HMGET %s %s %s", "0", buffer1, buffer2);
    value += atofc(reply->element[0]->str) * atofc(reply->element[1]->str);

    reply = redisCommand(c, "HMGET %s %s %s", "1", buffer1, buffer2);
    value += atofc(reply->element[0]->str) * atofc(reply->element[1]->str);

    reply = redisCommand(c, "HMGET %s %s %s", "2", buffer1, buffer2);
    value += atofc(reply->element[0]->str) * atofc(reply->element[1]->str);

    freeReplyObject(reply);
    return(value);
}

void readData(int start, int stop, int step, int N, redisContext *c)
{
    int timestep, particle, one, index, count, maxCount;
    float xVel, yVel, zVel;
    char xVelStr[256], yVelStr[256], zVelStr[256];
    FILE *fp;

    fp = fopen("./HISTORY_atoms/HISTORY_CLEAN_108", "r");
    if (fp == NULL)
        return;
    
    count = 1;
    maxCount = (stop-start) / step * N;

    redisReply *reply;

    char line[256];
    char buffer[256];
    while (fgets(line, sizeof line, fp) != NULL)
    {
        if(count == maxCount)
            return;
        sscanf(line, "%d %d %d %f %f %f", &timestep, &particle, &one, &xVel, &yVel, &zVel);
        index = (timestep - start) / step;

        sprintf(buffer, "%d:%d", particle, index);
        sprintf(xVelStr, "%f", xVel);
        reply = redisCommand(c, "HSET %s %s %s", "0", buffer, xVelStr);
        sprintf(buffer, "%d:%d", particle, index); 
        sprintf(yVelStr, "%f", yVel);        
        reply = redisCommand(c, "HSET %s %s %s", "1", buffer, yVelStr);
        sprintf(buffer, "%d:%d", particle, index);
        sprintf(zVelStr, "%f", zVel);         
        reply = redisCommand(c, "HSET %s %s %s", "2", buffer, zVelStr);
        freeReplyObject(reply);

        count++;
    }
    fclose(fp);
}

int main(int argc, char **argv)
{

    redisContext *c = redisConnect("127.0.0.1", 6379);
    if(c != NULL && c->err)
    {
        fprintf(stderr, "[Error] %s \n", c->errstr);
    }

    int start;
    int stop;
    int step;
    int tmax = 0;

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
        else
        {
            fprintf(stderr, "[Error] Command-line argument not recognized.\n");
            exit(-1);
        }
    }

    float xData, yData, zData;

    N++;                             // Number of particles
    M = ((stop - start) / step) + 1; // Number of timesteps

    tmax = (tmax == 0) ? (M / 3) : tmax;

    // padding(M, xData, yData, zData);
    // readData(start, stop, step, N, c);    
    float accumalate, particle;
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
                // particle += xData[i][t] * xData[i][t + dt] +
                //             yData[i][t] * yData[i][t + dt] +
                //             zData[i][t] * zData[i][t + dt];
                particle += retreive(i, t, dt, c);
            }
            count++;
            accumalate += particle;
        }
        accumalate /= ((N - 1) * count);
        fprintf(stdout, "%d, %e\n", dt, accumalate);
    }
    redisFree(c);
    return 0;
}