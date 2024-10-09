#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include <time.h>

//hypergraph
#define MAX_VERTICES 99 // number of vertices
#define MAX_EDGES 33 // number of hyperedges
#define Max_Gen 1 // number of subhypergraphs
#define Max_k 500 // number of largest degree, it should be larger than the max hyperdegree in the hypergraph

#define s 2 // selection intensity

char filename[] = "./UniformRing/uniform_ring_6_99_33.txt";
char datafile[] = "./UniformRing/data_uniform_ring_6_99_33_test.txt";
char name_program[] = "Uniform Ring";

int g = 100; // generation number on each subhypergraphs

#define TIME 1000 // number of simulation times
#define RUN 4e6  // number of rounds in each simulation
#define CRUN 1e3 // averaging the last 1e3 rounds to get fC

// Synergy factor
float r0 = 0.6, dr = 1, r_bound = 3;

// Adj Matrix of hypergraph
int Pic0[MAX_VERTICES][MAX_EDGES];

//hyperdegrees
int hd[Max_Gen][MAX_VERTICES];
int hd0[MAX_VERTICES];

//individual's hyperedges
int hl[Max_Gen][MAX_VERTICES][Max_k];
int hl0[MAX_VERTICES][MAX_EDGES];

//number of hyperedge's individuals
int n[MAX_EDGES];
int nc[MAX_EDGES];

// payoff of each individual
double payi[MAX_VERTICES];

int randnum[MAX_VERTICES];

int randset[MAX_EDGES];


int strategy[MAX_VERTICES];
int strategy0[MAX_VERTICES]; //initial condition for evolve same time


void matrix_read()
{
    // read the Adj matrix of temporal hypergraph
    FILE* fp;
    int i;
    int lines = 0;
    printf("Matrix reading\n");
    fp = fopen(filename, "r");
    if (fp == NULL)
    {
        return;
        printf("fail to open");
    }

    while (lines < MAX_VERTICES)
    {
        for (i = 0; i < MAX_EDGES; i++)
            fscanf(fp, "%d, ", &Pic0[lines][i]);
        if (feof(fp)) break;
        lines++;
    }

    fclose(fp);
}


void buildGraph()
{
    // constructing networks

    int numberOfVertices = MAX_VERTICES;

    // initializing
    int i, j, k;
    for (k = 0; k < Max_Gen; k++)
    {
        for (i = 0; i < MAX_VERTICES; i++)
        {
            hd0[i] = 0;
            hd[k][i] = 0;
            for (j = 0; j < MAX_EDGES; j++)
            {
                //Pic[k][i][j] = 0;
                n[j] = 0;
                nc[j] = 0;

            }
        }
    }

    printf("\n");
    //fprintf(fp,"\n");

    // read the list of hyperedges and hyperdegrees from Adj matrix   
    for (i = 0; i < MAX_VERTICES; i++)
    {
        for (j = 0; j < MAX_EDGES; j++)
        {
            k = Pic0[i][j];
            if (k > 0)//if they are linked
            {
                k--;
                hl0[i][hd0[i]] = j;
                hd0[i]++; //totoal hyperdegrees
                hl[k][i][hd[k][i]] = j;
                hd[k][i]++; //k-th subhyperdegree

                n[j] ++;
            }
        }

    }

    //print the hyperdegree and adj-list
    // int vertexCounter = 0;
    // int edgeCounter = 0;
    // printf("Hyperdegrees\n");
    // for (vertexCounter = 0; vertexCounter < numberOfVertices; vertexCounter++) {
    //     printf("%d\n", hd0[vertexCounter]);
    // }
    // printf("\n\n\n\n\n\nHypernetwork:\n");
    // for (vertexCounter = 0; vertexCounter < numberOfVertices; vertexCounter++) {
    //     printf("%d:\t", vertexCounter);
    //     for (edgeCounter = 0; edgeCounter < hd0[vertexCounter]; edgeCounter++) {

    //         printf("%d, ", hl0[vertexCounter][edgeCounter]);

    //     }
    //     printf("\n");
    // }
}


double paycal(int die, double r, int gen)
{
    // payoff calculation
    double pay = 0;
    int i, k, cond = 0;
    if (strategy[die] == 1) cond = 1;
    for (i = 0; i < hd[gen][die]; i++)
    {
        k = hl[gen][die][i];
        pay += nc[k] * r - cond;
    }

    return pay;
}


void swap(int* tempa, int* tempb)
{
    int temp;
    temp = *tempa;
    *tempa = *tempb;
    *tempb = temp;
}


void randomize(int l1)
{
    // randomize series
    int a[MAX_VERTICES], i = 0;
    while (i < l1)
    {
        a[i] = i;
        i++;
    }
    int j = 0;
    for (i = 0; i < (int)l1; i++)
    {
        j = rand() % (int)l1;
        swap(&a[i], &a[j]);
    }
    i = 0;
    while (i < l1)
    {
        randnum[i] = a[i];
        i++;
    }

}


int main()
{
    // Read the adjacency matrix
    matrix_read();

    //Build the hypernetwork
    printf("Simulation on tempotal hypernetwork\n\n");

    printf("%s\n", name_program);
    printf("Max_gen: %d\n", Max_Gen);

    buildGraph();

    srand(time(NULL));

    // number of cooperators
    int Ncoo = 0;

    // Iterators for steps and runs
    int i = 0, j = 0;
    long double j0 = 0;

    long double aa;
    float freC = 0, freCs = 0;

    int strategychange;
    int numberofSubHG;
    int role, nodeChosen;
    int hyperedgechosen;
    double payrole, payNodeChosen;
    int numberofCoplayer;


    FILE* fp;


    printf("started\n");

    double r = r0;
    while (r < r_bound) 
    {
        int vertexCounter = 0;
        int edgeCounter = 0;

        printf("%f", r);
        i = 0;
        freC = 0; // fraction of cooperators
        while (i < TIME)
        {
            // Initiation for strategies
            randomize(MAX_VERTICES);
            for (vertexCounter = 0; vertexCounter < MAX_VERTICES; vertexCounter++)
            {
                if (vertexCounter < MAX_VERTICES / 2) strategy[randnum[vertexCounter]] = 0; //0 stands for defectors
                else strategy[randnum[vertexCounter]] = 1;
            }

            for (edgeCounter = 0; edgeCounter < MAX_EDGES; edgeCounter++)
            {
                nc[edgeCounter] = 0; // number of cooperators in the hyperedge
                for (vertexCounter = 0; vertexCounter < MAX_VERTICES; vertexCounter++)
                    if (strategy[vertexCounter] == 1 && Pic0[vertexCounter][edgeCounter] > 0) 
                        nc[edgeCounter]++;
            }

            Ncoo = 0; // number of cooperators of the population
            for (vertexCounter = 0; vertexCounter < MAX_VERTICES; vertexCounter++)
                if (strategy[vertexCounter] == 1)  
                    Ncoo++;

            j = 0;
            freCs = 0; // fraction of cooperators in a single run
            while (j < RUN)
            {
                /*Simulation each round*/

                numberofSubHG = j / g % Max_Gen; // defining the number of subhypernetwork

                /*Calculate payoff of each node*/
                for (vertexCounter = 0; vertexCounter < MAX_VERTICES; vertexCounter++)
                {
                    payi[vertexCounter] = paycal(vertexCounter, r, numberofSubHG);
                    strategy0[vertexCounter] = strategy[vertexCounter];
                }

                /*synchronous update*/
                for (int v = 0; v < MAX_VERTICES; v++)
                {

                    nodeChosen = v; //the node chosen
                    if (hd[numberofSubHG][nodeChosen] == 0) continue; // ignore isolated node

                    // select a random hyperedge
                    hyperedgechosen = rand() % hd[numberofSubHG][nodeChosen];
                    //printf("\nl1: %d, hd0: %d\n", l1, hd0[die]);
                    hyperedgechosen = hl[numberofSubHG][nodeChosen][hyperedgechosen];
                    //printf("\nl1: %d, gc: %d\n", l1, gc[l1]);

                    /*Find the coplayer in hyperedge*/
                    int coplayer[MAX_VERTICES]; // list of the neighbours in the hyperedge
                    numberofCoplayer = 0;
                    for (vertexCounter = 0; vertexCounter < MAX_VERTICES; vertexCounter++)
                    {
                        if (vertexCounter == nodeChosen) continue;
                        if (Pic0[vertexCounter][hyperedgechosen] > 0)
                        {
                            coplayer[numberofCoplayer] = vertexCounter;
                            numberofCoplayer++;
                        }
                    }

                    // select a random neighbor as role model
                    role = rand() % numberofCoplayer;
                    role = coplayer[role];

                    // payoff of the two players
                    payrole = payi[role];
                    payNodeChosen = payi[nodeChosen];

                    aa = rand() / (RAND_MAX + 1.0);
                    if (aa < 1 / (1 + exp(-s * (payrole - payNodeChosen))))
                    {
                        strategychange = strategy0[role] - strategy0[nodeChosen];
                        if (strategychange != 0)
                            for (edgeCounter = 0; edgeCounter < hd0[nodeChosen]; edgeCounter++)
                                nc[hl0[nodeChosen][edgeCounter]] += strategychange;// update nc
                        Ncoo += strategychange;
                        strategy[nodeChosen] = strategy0[role]; // change strategy
                    }
                }

                if (Ncoo == 0)
                {
                    freCs = 0;
                    break;
                }
                else if (Ncoo == MAX_VERTICES)
                {
                    freCs = 1;
                    break;
                }

                // To take average for the last CRUN runs
                j0 = j - CRUN;
                if (j0 >= 0)
                    freCs = freCs * (long double)j0 / (j0 + 1) + (double)Ncoo / MAX_VERTICES * (long double)1 / (j0 + 1);
                j++;
            }

            // printf("%d, %f; ", i, freCs);
            freC = freC * (double)i / (i + 1) + freCs * (double)1 / (i + 1); // fraction of cooperators
            i++;
        }

        printf(", %f", freC);
        fp = fopen(datafile, "a");
        fprintf(fp, "%.2f, %f\n", r, freC);
        fclose(fp);

        printf("\n");

        // the gap of r
        dr = 0.02;
        if (freC > 0 && freC < 1)
            dr = 0.04;
        else if (freC == 1.0)
            dr = 1;
        r += dr;
    }

    return 0;
}