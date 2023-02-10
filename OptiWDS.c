/*OptiWDS Thesis by Felipe R. Goleta III for CE 199*/
/*This has been created with the help of Engr. Guimba for his OptiPipe and Louie Gallos for debugging and consultation and being my friend*/
/*This was last edited on May 1, 2022*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <conio.h>
#include <omp.h>
#include <mem.h>
#include <string.h>
#include "EPANET_2_2_WIN_32_64/include/epanet2_2.h"         //Contains the header file for the execution of EPANET 2.2 dll
#include "EPANET_2_2_WIN_32_64/include/epanet2_enums.h"

typedef struct PipeData {           //Input data of commercially available pipes
    int pipeID;                     //Decimal ID to convert pipe size to binary coding
    double pipeDiam;                //Diameter of commercially available pipes in mm
    double pipeUcost;               //Unit cost of commercially available pipes per linear meter in Pesos
    double pipeRough;               //Roughness Coefficient of commercially available pipes
} PipeData;

typedef struct PumpData {           //Input data of commercially available pumps
    int pumpID;                     //Decimal ID to convert pipe size to binary coding
    double pumpPower;               //Power rating of commercially available pumps in kiloWatt
    double pumpUcost;               //Unit cost of commercially available pumps
} PumpData;

typedef struct Node {               //Node parameters where water goes out
    int nodeIndex;                  //Index Name of Node
    char epanetID[100];             //Node ID in EPANET
    double pressure;                //Node pressure
    double ppenalty;                //Pressure Penalty for the node
} Node;

typedef struct Pipe {               //Pipe parameters for water medium
    int pipeID;                     //Pipe ID based on Pump Data
    int pipeIndex;                  //Index Number of the pipe based on EPANET
    char epanetID[100];             //EPANET ID name of the Pipe
    double diam;                    //Pipe diameter
    double ucost;                   //Unit Cost per Pipe based on material only; Real Cost
    double rough;                   //Hazen-Williams coefficient for pipe roughness
    double length;                  //Pipe length
    Node* ends[2];                 //Start and End nodes of the pipe
    double velocity;                //Ave. Velocity of Water in pipe
    double totcost;                 //Total Cost of all Pipes
    double realcost;                //Unit Cost of the Pipe multiplied to the length of the pipe
    double vpenalty;                //Penalty Cost for violating Velocity Constraints
} Pipe;

typedef struct Pump {               //Pump parameters for water medium
    int pumpID;                     //Pump ID based on Pump Data
    int pumpIndex;                  //Index Number of the pump based on EPANET
    char epanetID[100];             //EPANET ID name of the Pump
    double power;                   //Constant Power Rating of Pump
    double ucost;                   //Unit Cost of the Pump based on power rating
    double realcost;                //Unit Cost of the Pump based on power rating
} Pump;

typedef struct chrom {              //chromosome array of one solution
    Pipe** pipes;                  //Array of pointers to pipes
    Pump** pumps;                  //Array of pointers to pumps
    double fitness;                 //Fitness point of the chromosome
    double pipeRealCost;            //Total real cost of all pipes
    double pumpRealCost;            //Total real cost of all pumps
    double pipeTotalCost;           //Total real cost + penalty cost of all pipes
    double totcost;                 //(Total real cost + penalty cost of all pipe) + (Total real cost of all pumps)
} chrom;

/*GLOBAL VARIABLES*/
FILE* pipeFiles;                   //File pointer to open the text file containing the commercially available pipes
FILE* pumpFiles;                   //File pointer to open the text file containing the commercially available pumps
FILE* resultFiles;                 //File pointer to open the text file where the results of GA will be presented
FILE* ptr1;                        //Dummy file pointer used to measure the number of pipes for the dynamic memory allocation of pipedata
FILE* ptr2;                        //Dummy file pointer used to measure the number of pumps for the dynamic memory allocation of pumpdata
PipeData* pipedata;                //Storage of the data of commercially available pipes
PumpData* pumpdata;                //Storage of the data of commercially available pumps
int pipeCount;                      //Variable used to count the number of pipes for the dynamic memory allocation of pipedata
int pumpCount;                      //Variable used to count the number of pumps for the dynamic memory allocation of pumpdata
int pipeStock;                      //Variable used to count the number of pipes for storing the attributes of the pipes in pipedata
int pumpStock;                      //Variable used to count the number of pumps for storing the attributes of the pumps in pumpdata

/*Genetic Algorithm Parameters*/
int trials;                         //Number of GA trials
int gensize;                        //Generation Size or Number of Solutions to be generated for each generation
int iterations;                     //Maximum number of iterations to end the program for each trial
double mrate;                       //Mutation Rate (%) or possibility of mutation per crossover
int nElite;                         //Number of Elite Individuals that will not undergo crossover but will be kept in the next Generation
int nDiscarded;                     //Number of Discarded Individuals that will be discarded in the current Generation so that next Generation will not inherit their weak genes
double nElitePercentage;            //Number of Elites Percentage with respect to the generation size
double nDiscardedPercentage;        //Number of Discarded Percentage with respect to the generation size
double pmin;                        //Min allowable Pressure in nodes (m)
double pmax;                        //Max allowable Pressure in nodes (m)
double vmin;                        //Min allowable Velocity in pipes (m/s)
double vmax;                        //Max allowable Velocity in pipes (m/s)
double PP1;                         //Value of PP1 in Penalty Function Equations
double PP2;                         //Value of PP2 in Penalty Function Equations
double VP1;                         //Value of VP1 in Penalty Function Equations
double VP2;                         //Value of VP2 in Penalty Function Equations
int errcode;                        //Numerical representation of Errors in EPANET
char inpFile[100];                  //Name of the Input Pipe File
char rptFile[100];                  //Name of the Output Report File
char outFile[100];                  //Name of the Optional binary output file. Storage of hydraulic results
int nNodes;                         //Number of nodes (junctions + tanks +reservoirs) in the system
int nLinks;                         //Number of links (pipes + valves + pumps) in the system
int nPipes;                         //Number of pipes in the system
int nPumps;                         //Number of pumps in the system
int nTanks;                         //Number of tanks + reservoirs in the system
int linkType;                       //Returns the link type of the link: Pipe = 1, Pump = 2,
int nodeType;                       //Returns the node type of the node: Junction = 0, Reservoir = 1, Tank = 2
EN_Project ph;                      //EPANET Project File as specified in OWA-EPANET 2.2 Toolkit

/*USER-DEFINED FUNCTIONS*/
void generateChromosomes(chrom** popGen);                          //Generate permutations of pipes and pumps for all the chromosomes in the beginning of each trial of GA
double penaltycost(chrom* pop);                                    //Computes penalty cost of the chromosome or solution when real pressure/velocity exceeds the boundary condition as specified in pmin,pmax, vmin, and vmas
double fit(chrom* pop, int iter);                                  //Computes the fitness of a chromosome
void fitnessRanking(chrom** pop, chrom** dummy);                  //Ranks the chromosomes per iteration based on their fitness level
int roulette(chrom** popGen);                                      //Weighted randomizer for parent selection
void crossover(chrom** currentGen, chrom** nextGen, int iter);    //Crosses the chosen parents to get the offspring for the next Generation
void finalRanking(chrom** pop, chrom** dummy, int limit);              //Ranks the best solutions in each GA trial

int main() {
    char pipeTxt[100];                          //Name of the text file containing the commercially available pipes
    char pumpTxt[100];                          //Name of the text file containing the commercially available pumps
    char reportTxt[100];                        //Name of the text file containing the results and reports of GA
    char dummy1[100];                           //Dummy variable required by fgets for counting the number of lines in pipeFiles and Memory Allocation for pipedata
    char dummy2[100];                           //Dummy variable required by fgets for counting the number of lines in pumpFiles and Memory Allocation for pumpdata
    char a1[20], a2[20], a3[20], a4[20];        //Parameters in the beginning of the text file containing the commercially available pipes  (Pipe Number, Diameter, Unit Cost, and Roughness)
    char b1[20], b2[20], b3[20];                //Parameters in the beginning of the text file containing the commercially available pumps  (Pump Number, Power Rating, and Unit Cost)
    srand(time(NULL));                          //Initialization of number randomizer
    clock_t t;                                  //Store time to get the time spent using GA
    double runtime, sec;                        //Variables for getting the time spent in conducting GA
    int hr, mins;                               //Variables for getting the time spent in conducting GA
    int printNumber;                            //Variable used to store the number of solutions to be reported in the Report File
    int n, m, linktype;                         //Dummy variables to store the count of pipes and pumps, and type of link (1=pipe & 2=pump), respectively
    int terminationTimes;                       //Variable used to store the number of times iteration will still continue even if GA converges prematurely. Once this is exceeded, termination criterion will be triggered and iteration will be stopped and a new trial will begin
    int counter;                                //Count variable for iterations happening before terminationTimes is triggered
    double nodePressure;                        //Variable for storing the pressure of a node for the best solution
    double pipeVelocity;                        //Variable for storing the velocity of a pipe for the best solution
    char nodeEPANETID[20];                      //String name for the EPANET ID of the node for the best solution
    char pipeEPANETID[20];                      //String name for the EPANET ID of the pipe for the best solution
    pipeCount = 0; pumpCount = 0; pipeStock = 0; pumpStock = 0;

    printf("Welcome to OptiWDS:\n\n");
    /*Passing Data Input of commercially available pipes*/
    do {

        printf("Please enter filename for commercially available pipes (include.txt): ");
        scanf("%s", pipeTxt);
        pipeFiles = fopen(pipeTxt, "r");

        if (pipeFiles == NULL) {
            printf("\nError Reading File!\n\n");
        }
    } while (pipeFiles == NULL);

    /*Counting the number of lines in pipeFiles and Memory Allocation for pipedata*/
    ptr1 = fopen(pipeTxt, "r");
    do {
        fgets(dummy1, 100, ptr1);
        pipeCount++;
    } while (!feof(ptr1));
    pipeCount--;
    fclose(ptr1);
    pipedata = (PipeData*)calloc(pipeCount, sizeof(PipeData));
    rewind(pipeFiles);

    /*Reading and Writing Data Input of commercially available pipes*/
    fscanf(pipeFiles, "%s\t\t%s\t%s\t%s", a1, a2, a3, a4);
    printf("\nPipeID\t\tDiameter(m)\tUCost(Php/m)\tRoughnes(C)\n");
    while (!feof(pipeFiles)) {
        fscanf(pipeFiles, "%d\t\t%lf\t\t%lf\t\t%lf", &pipedata[pipeStock].pipeID, &pipedata[pipeStock].pipeDiam, &
            pipedata[pipeStock].pipeUcost, &pipedata[pipeStock].pipeRough);
        printf("%d\t\t%.1lf\t\t%.2lf\t\t%.2lf\n", pipedata[pipeStock].pipeID, pipedata[pipeStock].pipeDiam,
            pipedata[pipeStock].pipeUcost, pipedata[pipeStock].pipeRough);
        pipeStock++;
    }
    fclose(pipeFiles);

    /*Passing Data Input of commercially available pumps*/
    do {
        printf("\nPlease enter filename for commercially available pumps (include.txt): ");
        scanf("%s", pumpTxt);
        pumpFiles = fopen(pumpTxt, "r");

        if (pumpFiles == NULL) {
            printf("\nError Reading File!\n\n");
        }
    } while (pumpFiles == NULL);

    /*Counting the number of lines in pumpFiles and Memory Allocation for pumpdata*/
    ptr2 = fopen(pumpTxt, "r");
    do {
        fgets(dummy2, 100, ptr2);
        pumpCount++;
    } while (!feof(ptr2));
    pumpCount--;
    fclose(ptr2);
    pumpdata = (PumpData*)calloc(pumpCount, sizeof(PumpData));
    rewind(pumpFiles);

    /*Reading and Writing Data Input of commercially available pump*/
    fscanf(pumpFiles, "%s\t\t%s\t\t%s", b1, b2, b3);
    printf("\nPumpID\t\tPower(kW)\t\tUCost(Php) \n");
    while (!feof(pumpFiles)) {
        fscanf(pumpFiles, "%d\t\t%lf\t\t\t%lf", &pumpdata[pumpStock].pumpID, &
            pumpdata[pumpStock].pumpPower, &pumpdata[pumpStock].pumpUcost);
        printf("%d\t\t%.3lf\t\t\t%.2lf\n", pumpdata[pumpStock].pumpID,
            pumpdata[pumpStock].pumpPower, pumpdata[pumpStock].pumpUcost);
        pumpStock++;
    }
    fclose(pumpFiles);

    /*Request for the User-defined Variables*/
    printf("\nGENETIC ALGORITHM PARAMETERS:\n\n");
    printf("Number of Trials (Recommended Value by Guimba:20): ");
    scanf("%d", &trials);
    printf("Generation Size (Recommended Value by Guimba:200): ");
    scanf("%d", &gensize);
    printf("Number of GA Iterations (Recommended Value by Guimba:1500): ");
    scanf("%d", &iterations);
    printf("Mutation Rate in Percentage (Recommended Value by Guimba:5): ");
    scanf("%lf", &mrate);
    printf("Elite Individual Percentage with respect to Generation Size(Recommended Value by Guimba:4): ");
    scanf("%lf", &nElitePercentage);
    printf("Discarded Individual (Percentage with respect to Generation Size(Recommended Value by Guimba:10): ");
    scanf("%lf", &nDiscardedPercentage);
    printf("\nPROBLEM CONDITIONS:\n\n");
    printf("Minimum Required Pressure in meters (Recommended Value by RWSDM:3) ");
    scanf("%lf", &pmin);
    printf("Maximum Allowable Pressure in meters (Recommended Value by RWSDM:70): ");
    scanf("%lf", &pmax);
    printf("Minimum Required Velocity in meters per second (Recommended Value by RWSDM:0.45): ");
    scanf("%lf", &vmin);
    printf("Maximum Allowable Velocity in meters per second (Recommended Value by RWSDM:3): ");
    scanf("%lf", &vmax);
    printf("Pressure Penalty 1 (Recommended Value by Guimba:0.1): ");
    scanf("%lf", &PP1);
    printf("Pressure Penalty 2 (Recommended Value by Guimba:10): ");
    scanf("%lf", &PP2);
    printf("Velocity Penalty 1 (Recommended Value by Guimba:1): ");
    scanf("%lf", &VP1);
    printf("Velocity Penalty 2 (Recommended Value by Guimba:5): ");
    scanf("%lf", &VP2);
    printf("Up to how many times do you wish to proceed before terminating the iterations when the solution converges prematurely?: ");
    scanf("%d", &terminationTimes);

    nElite = gensize * (nElitePercentage / 100);                                //Computes the rounded number of elite chromosomes based on user-defined Elite Percentage
    nDiscarded = gensize * (nDiscardedPercentage / 100);                        //Computes the rounded number of discared chromosomes based on user-defined Discarded Percentage

    /*Enter EPANET Files*/
    EN_createproject(&ph);                                                    //OWA-EPANET API Function first used before any API Functions can be used
    printf("Please enter the EPANET File names.\nMake sure that the default in your EPANET File is in Metric System and has a Hazen-William Equation for the Headloss Computation.\n");
    printf("Also make sure that the EPANET input file is in single period analysis\n\n");
    printf("\nPlease enter EPANET MODEL FILES:\n\n");
    do {
        printf("Filename of Input File (include .inp): ");
        scanf("%s", inpFile);
        printf("Desired Filename for Report File (include .rpt): ");
        scanf("%s", rptFile);
        printf("Desired Filename for Binary File (include .out): ");
        scanf("%s", outFile);
        errcode = EN_open(ph, inpFile, rptFile, outFile);                       //Opens the EPANET Input File and saves results in Report and Output File
        if (errcode > 0) {                                                      //If there are no errors, then the program will continue, otherwise it will loop and ask again
            printf("\nFile cannot be opened! [Error %d]\n", errcode);
        }
    } while (errcode > 0);
    printf("\n\nFile successfully opened!\n\n");
    /*Enter the Name the Report Text File you want the report of GA to be written*/
    printf("Please enter text file name where the results of GA will be stored(include.txt): ");
    scanf("%s", reportTxt);
    resultFiles = fopen(reportTxt, "w");

    /*Getting the number of pipes,pumps, and nodes*/
    nPipes = 0, nPumps = 0, linkType = 0, nodeType = 0;
    EN_getcount(ph, EN_NODECOUNT, &nNodes);                                    //EPANET Function to get the count of nodes(junctions + tanks + reservoirs) in the EPANET input file
    EN_getcount(ph, EN_TANKCOUNT, &nTanks);                                    //EPANET Function to get the count of tanks + reservoirs in the EPANET input file
    EN_getcount(ph, EN_LINKCOUNT, &nLinks);                                    //EPANET Function to get the count of links(pipes + pumps + valves...) in the EPANET input file
    for (int i = 1; i <= nLinks; i++) {                                         //EPANET Function to classify the link type of the ith link (1=Pipes, 2=Pumps)
        EN_getlinktype(ph, i, &linkType);
        if (linkType == 1) {
            nPipes++;
        }
        if (linkType == 2) {
            nPumps++;
        }
    }

    /*Shows the number of Junctions, Nodes, Tanks, Links, Pipes, and Pumps in the Input File*/
    printf("\n\nWe have the following number of WDS components in your inputted EPANET FIle:\n\n");
    printf("Number of Nodes: %d\n", nNodes);
    printf("Number of Tanks and Reservoirs: %d\n", nTanks);
    printf("Number of Links: %d\n", nLinks);
    printf("Number of Pipes: %d\n", nPipes);
    printf("Number of Pumps: %d\n", nPumps);

    /*Chromosome data struct variable initialization*/
    chrom** currentGen;                                                        //Parent Generation of Chromosomes
    chrom** nextGen;                                                           //Offspring Generation of Chromosomes
    chrom** topGen;                                                            //Chromosomes containing the best Chromosomes for each trial
    chrom** dummy;                                                             //Dummy struct variable for fitnessRanking and finalRanking functions

    /*Dynamic Memory Allocation for Struct Variables and their respective attributes*/
    topGen = (chrom**)malloc(gensize * sizeof(chrom*));
    currentGen = (chrom**)malloc(gensize * sizeof(chrom*));
    nextGen = (chrom**)malloc(gensize * sizeof(chrom*));
    dummy = (chrom**)malloc(gensize * sizeof(chrom*));
    for (int i = 0; i < gensize; i++) {
        topGen[i] = (chrom*)malloc(sizeof(chrom));
        currentGen[i] = (chrom*)malloc(sizeof(chrom));
        nextGen[i] = (chrom*)malloc(sizeof(chrom));
        dummy[i] = (chrom*)malloc(sizeof(chrom));
    }
    for (int i = 0; i < trials; i++) {
        topGen[i]->pipes = (Pipe**)malloc(nPipes * sizeof(Pipe*));
        topGen[i]->pumps = (Pump**)malloc(nPumps * sizeof(Pump*));
        for (int j = 0; j < nPipes; j++) {
            topGen[i]->pipes[j] = (Pipe*)malloc(sizeof(Pipe));
            topGen[i]->pipes[j]->ends[0] = (Node*)malloc(sizeof(Node));
            topGen[i]->pipes[j]->ends[1] = (Node*)malloc(sizeof(Node));
        }
        for (int j = 0; j < nPumps; j++) {
            topGen[i]->pumps[j] = (Pump*)malloc(sizeof(Pump));
        }
    }
    for (int i = 0; i < gensize; i++) {
        currentGen[i]->pipes = (Pipe**)malloc(nPipes * sizeof(Pipe*));
        currentGen[i]->pumps = (Pump**)malloc(nPumps * sizeof(Pump*));
        nextGen[i]->pipes = (Pipe**)malloc(nPipes * sizeof(Pipe*));
        nextGen[i]->pumps = (Pump**)malloc(nPumps * sizeof(Pump*));
        dummy[i]->pipes = (Pipe**)malloc(nPipes * sizeof(Pipe*));
        dummy[i]->pumps = (Pump**)malloc(nPumps * sizeof(Pump*));
        for (int j = 0; j < nPipes; j++) {
            currentGen[i]->pipes[j] = (Pipe*)malloc(sizeof(Pipe));
            nextGen[i]->pipes[j] = (Pipe*)malloc(sizeof(Pipe));
            dummy[i]->pipes[j] = (Pipe*)malloc(sizeof(Pipe));
            currentGen[i]->pipes[j]->ends[0] = (Node*)malloc(sizeof(Node));
            currentGen[i]->pipes[j]->ends[1] = (Node*)malloc(sizeof(Node));
            nextGen[i]->pipes[j]->ends[0] = (Node*)malloc(sizeof(Node));
            nextGen[i]->pipes[j]->ends[1] = (Node*)malloc(sizeof(Node));
            dummy[i]->pipes[j]->ends[0] = (Node*)malloc(sizeof(Node));
            dummy[i]->pipes[j]->ends[1] = (Node*)malloc(sizeof(Node));
        }
        for (int j = 0; j < nPumps; j++) {
            currentGen[i]->pumps[j] = (Pump*)malloc(sizeof(Pump));
            nextGen[i]->pumps[j] = (Pump*)malloc(sizeof(Pump));
            dummy[i]->pumps[j] = (Pump*)malloc(sizeof(Pump));
        }
    }

    /*Running Genetic Algorithm*/
    printf("\nRunning Genetic Algorithm\n");
    t = clock();                                                                                                          //Gets the starting time of GA
    for (int i = 0; i < trials; i++) {
        printf("\nFor Trial %d:\n", i + 1);
        generateChromosomes(currentGen);                                                                                //Generating Starting Chromosomes for each trial using the commercially available pipes and pumps
        fitnessRanking(currentGen, dummy);                                                                              //Ranks the fitness of the Generated Chromosomes
        counter = 0;
        for (int j = 0; j < iterations; j++) {
            for (int k = 0; k < gensize; k++) {
                nextGen[k] = currentGen[k];                                                                             //The next generated chromosomes get the permanend data from their parents (e.g, node and link index, pipe length, etc.)
            }
            crossover(currentGen, nextGen, j + 1);                                                                      //Selects parents the Parent Pool and crosses their genes until a new offspring is created. Mutation is probable. This will happen until all the pool of offspring generation has been filled up
            fitnessRanking(nextGen, dummy);                                                                             //Ranks the fitness of the Offspring Chromosomes
            if (currentGen[0]->totcost == nextGen[0]->totcost) {                                                    //If the CurrentGen Total cost is the same for an iteration, there will be an additional value to the counter
                counter++;
            }
            else {
                counter = 0;
            }
            for (int k = 0; k < gensize; k++) {
                currentGen[k] = nextGen[k];                                                                             //The offspring chromosomes mature here that they become eligible as parents
            }
            printf("Iteration:%d \tBest Total Cost: Php %lf\n", j + 1, nextGen[0]->totcost);
            if (counter == terminationTimes) {                                                                         //Once termination Times is triggered, termination criterion for the iteration will commenced and the iteration will end and the new trial will start
                break;
            }
        }
        topGen[i] = nextGen[0];                                                                                         //The best offspring in a trial of many iterations is saved on the top generated chromosomes
        printf("\nThe Best solution for this iteration:\n");
        printf("Trial [%d] Real Cost: Php %.2lf\n\n", i + 1, topGen[i]->totcost);
    }
    finalRanking(topGen, dummy, trials);                                                                                //Ranks all of the best chromosomes in every trial in increasing amount of total cost
    printf("\nRanking\t\tPipe Real Cost(Php)\t\tPump Real Cost(Php)\t\tTotal Real Cost(Php)\t\tTotal Cost of the System(Php)\n");
    double totalRealCost[trials];                                                                                               //Place holder for Total Real Cost of each Top Gen chromosomes
    for (int i = 0; i < trials; i++) {                                                                                  //Solving for the Total Real Cost of each of the best chromosomes in Top Gen
        totalRealCost[i] = topGen[i]->pipeRealCost + topGen[i]->pumpRealCost;
        printf("\n%d\t\t%.2lf\t\t\t%.2lf\t\t\t%.2lf\t\t\t%.2lf", i + 1, topGen[i]->pipeRealCost, topGen[i]->pumpRealCost, totalRealCost[i], topGen[i]->totcost);
    }

    /*Gets the Time at the end of GA and solves for the time spent running GA*/
    t = clock() - t;
    runtime = (double)t / CLOCKS_PER_SEC;
    hr = (int)runtime / 3600;
    mins = (int)(runtime - (3600 * hr)) / 60;
    sec = (double)(runtime - (3600 * hr) - (60 * mins));
    fprintf(resultFiles, "\n\nProgram took %d hrs %d mins %.3lf seconds to execute.", hr, mins, sec);

    printf("\n\nHow many top solutions do you want to be written in the Report File? ");
    scanf("%d", &printNumber);

    /*Printing the results to the Report File*/
    fprintf(resultFiles, "These are the Genetic Algorithm and Problem Conditions you have used in this setup: \n\n");
    fprintf(resultFiles, "Number of Trials: %d\n", trials);
    fprintf(resultFiles, "Number of Generation Size: %d\n", gensize);
    fprintf(resultFiles, "Number of GA Iterations: %d\n", iterations);
    fprintf(resultFiles, "Mutation Rate Percentage: %.2lf%%\n", mrate);
    fprintf(resultFiles, "Elite Individual Percentage: %.2lf%%\n", nElitePercentage);
    fprintf(resultFiles, "Discarded Individual Percentage: %.2lf%%\n", nDiscardedPercentage);
    fprintf(resultFiles, "Min Pressure: %.2lf\n", pmin);
    fprintf(resultFiles, "Max Pressure: %.2lf\n", pmax);
    fprintf(resultFiles, "Min Velocity: %.2lf\n", vmin);
    fprintf(resultFiles, "Max Velocity: %.2lf\n", vmax);
    fprintf(resultFiles, "Pressure Penalty 1: %.2lf\n", PP1);
    fprintf(resultFiles, "Pressure Penalty 2: %.2lf\n", PP2);
    fprintf(resultFiles, "Velocity Penalty 1: %.2lf\n", VP1);
    fprintf(resultFiles, "Velocity Penalty 2: %.2lf\n", VP2);
    fprintf(resultFiles, "The number of iterations before preliminary termination criteria is reached: %d\n\n", terminationTimes);
    fprintf(resultFiles, "The Top %d Best Trials are as follows:\n\n", printNumber);
    fprintf(resultFiles, "\nRanking\t\tPipe Real Cost(Php)\t\tPump Real Cost(Php)\t\tTotal Real Cost(Php)\t\tTotal Cost of the System(Php)\n");
    for (int i = 0; i < printNumber; i++) {
        fprintf(resultFiles, "\n%d\t\t%.2lf\t\t\t%.2lf\t\t\t%.2lf\t\t\t%.2lf\n", i + 1, topGen[i]->pipeRealCost, topGen[i]->pumpRealCost, totalRealCost[i], topGen[i]->totcost);
    }
    fprintf(resultFiles, "\n\nPipesThe Pipe Sizing and Power Rating of each EPANET ID and Index are as follows:\n\n\n\n");
    for (int i = 0; i < printNumber; i++) {
        fprintf(resultFiles, "For Top %d Chromosome:\n\n", i + 1);
        fprintf(resultFiles, "Pipe ID\t\tEPANETID\t\tPipe Diameter (mm)\t\tRoughness(C)\n\n");
        for (int j = 0; j < nPipes; j++) {
            fprintf(resultFiles, "%d\t\t%s\t\t\t%.2lf\t\t\t\t%.0lf\n", topGen[i]->pipes[j]->pipeIndex, topGen[i]->pipes[j]->epanetID, topGen[i]->pipes[j]->diam, topGen[i]->pipes[j]->rough);
        }
        fprintf(resultFiles, "\nPump ID\t\tEPANETID\tPower Rating (kW)\n\n");
        for (int j = 0; j < nPumps; j++) {
            fprintf(resultFiles, "%d\t\t%s\t\t%.2lf\n", topGen[i]->pumps[j]->pumpIndex, topGen[i]->pumps[j]->epanetID, topGen[i]->pumps[j]->power);
        }
        fprintf(resultFiles, "\n\n");
    }
    fprintf(resultFiles, "\n\nPressure and Velocity of each node and pipe for the best solution are as follows:\n\n");
    fprintf(resultFiles, "Node Index\t\tEPANET ID\t\tPressure(m)\n\n");
    for (int i = 1; i <= nNodes; i++) {
        EN_getnodevalue(ph, i, EN_PRESSURE, &nodePressure);
        EN_getnodeid(ph, i, nodeEPANETID);
        fprintf(resultFiles, "%d\t\t\t%s\t\t\t%.2lf\n", i, nodeEPANETID, nodePressure);
    }
    fprintf(resultFiles, "\n\nPipe Index\tEPANET ID\tVelocity(m/s)\n\n");
    for (int i = 1; i <= nPipes; i++) {
        EN_getlinkvalue(ph, i, EN_VELOCITY, &pipeVelocity);
        EN_getlinkid(ph, i, pipeEPANETID);
        fprintf(resultFiles, "%d\t\t%s\t\t%.2lf\n", i, pipeEPANETID, pipeVelocity);
    }

    /*Saves the best solution to an EPANET Input File*/
    errcode = EN_open(ph, inpFile, rptFile, outFile);
    n = 0; m = 0;
    for (int i = 1; i <= nLinks; i++) {
        EN_getlinktype(ph, i, &linktype);
        if (linktype == 1) {
            EN_setlinkvalue(ph, i, EN_DIAMETER, topGen[0]->pipes[n]->diam);
            EN_setlinkvalue(ph, i, EN_ROUGHNESS, topGen[0]->pipes[n]->rough);
            n++;
        }
        if (linktype == 2) {
            EN_setlinkvalue(ph, i, EN_PUMP_POWER, topGen[0]->pumps[m]->power);
            m++;
        }
    }
    EN_saveinpfile(ph, "NewInputFile.inp");
    EN_close(ph);                                                                                                       //EPANET Function to close the project
    EN_deleteproject(ph);                                                                                               //EPANET Function to delete and free up memory for the project
    printf("\n\nGA has been finished. You can now see the Report File entitled \"%s\" and the New EPANET Input File for the best solution entitled \"NewInputFile.inp\" in the folder of OptiWDS. \nThank you for using our program!\n\n", reportTxt);
    return 0;
}

/*Generation of Solutions and Solves for the Cost and Fitness of the Generated Solutions*/
void generateChromosomes(chrom** popGen) {
    int randomPipeStock;                    //Index Randomizer in selecting Commercially Available Pipes for each chromosome
    int randomPumpStock;                    //Index Randomizer in selecting Commercially Available Pumps for each chromosome
    int nodes[nPipes][2];                   //Variable for storing the start and end node index of each pipe
    int pipeIndex[nPipes];                  //Variable storing the EPANET index of the pipes of the input file
    int pumpIndex[nPumps];                  //Variable storing the EPANET index of the pumps of the input file
    char epanetPipeID[nPipes][20];          //Variable storing the EPANET ID character of the pipes of the input file
    char epanetPumpID[nPumps][20];          //Variable storing the EPANET ID character of the pumps of the input file
    int n = 0, m = 0;                       //Dummy Variables for counting the pipes and pumps in the input file, respectively
    double length[nPipes];                  //Array for initial storage of lengths of each pipe index in the input file
    double pipeTotalRealCost;               //Variable storing the total real cost of pipes
    double pumpTotalRealCost;               //Variable storing the total real cost of pumps

    /*Categorizing link types for each link into pipes and pumps and using API Functions to get relevant data from input file*/
    for (int i = 1; i <= nLinks; i++) {
        EN_getlinktype(ph, i, &linkType);
        if (linkType == 1) {
            pipeIndex[n] = i;
            EN_getlinknodes(ph, i, &(nodes[n][0]), &(nodes[n][1]));
            EN_getlinkvalue(ph, i, EN_LENGTH, &(length[n]));
            EN_getlinkid(ph, i, epanetPipeID[n]);
            n++;
        }
        if (linkType == 2) {
            pumpIndex[m] = i;
            EN_getlinkid(ph, i, epanetPumpID[m]);
            m++;
        }
    }

    /*Generation of Chromosomes from 1 to Generation Size based on random selection from commercially available pipes and pumps*/
    for (int i = 0; i < gensize; i++) {
        pipeTotalRealCost = 0; pumpTotalRealCost = 0;
        for (int j = 0; j < nPipes; j++) {
            randomPipeStock = (rand() % (pipeStock)) + 1;                                                                   //Randomizes number from 1 to the number of commercially available Pipes
            popGen[i]->pipes[j]->pipeIndex = pipeIndex[j];                                                              //Equating the Link Index Number from EPANET to the chromosome
            popGen[i]->pipes[j]->ends[0]->nodeIndex = nodes[j][0];                                                    //Equating the Start Node of Link Index from EPANET to the chromosome
            popGen[i]->pipes[j]->ends[1]->nodeIndex = nodes[j][1];                                                    //Equating the End Node of Link Index from EPANET  to the chromosome
            EN_getnodeid(ph, popGen[i]->pipes[j]->ends[0]->nodeIndex, popGen[i]->pipes[j]->ends[0]->epanetID);  //EPANET Function to get the EPANET ID of the start nodes
            EN_getnodeid(ph, popGen[i]->pipes[j]->ends[1]->nodeIndex, popGen[i]->pipes[j]->ends[1]->epanetID);  //EPANET Function to get the EPANET ID of the end nodes
            popGen[i]->pipes[j]->length = length[j];                                                                    //Equating the Length of PLink Index from EPANET  to the chromosome
            popGen[i]->pipes[j]->pipeID = pipedata[randomPipeStock - 1].pipeID;                                         //Equating the inputs of pipe i index to the pipe ID of commercially available pipes
            popGen[i]->pipes[j]->diam = pipedata[randomPipeStock - 1].pipeDiam;                                         //Equating the diameter of pipe i to the diameter of commercially available pipes
            popGen[i]->pipes[j]->ucost = pipedata[randomPipeStock - 1].pipeUcost;                                       //Equating the unit cost of pipe i to the unit cost of commercially available pipes
            popGen[i]->pipes[j]->rough = pipedata[randomPipeStock - 1].pipeRough;                                       //Equating the roughness of pipe i to the roughness of commercially available pipes
            popGen[i]->pipes[j]->realcost = popGen[i]->pipes[j]->ucost * popGen[i]->pipes[j]->length;           //Real Cost of this pipe based on Ucost * pipe length
            strcpy(popGen[i]->pipes[j]->epanetID, epanetPipeID[j]);                                                     //Copies the EPANET Pipe ID to the pipe attribute of pipe structure
            pipeTotalRealCost += popGen[i]->pipes[j]->realcost;                                                         //Adds all the total real cost of the pipes for all the pipes in the input file
        }
        for (int j = 0; j < nPumps; j++) {
            randomPumpStock = (rand() % (pumpStock)) + 1;                                                                   //Randomizes number from 1 to the number of commercially available Pumps
            popGen[i]->pumps[j]->pumpIndex = pumpIndex[j];                                                              //Equating the Link Index Number from EPANET to the chromosome
            popGen[i]->pumps[j]->pumpID = pumpdata[randomPumpStock - 1].pumpID;                                         //Equating the inputs of pump i index to the pump ID of commercially available pumps
            popGen[i]->pumps[j]->power = pumpdata[randomPumpStock - 1].pumpPower;                                       //Equating the power rating of pump i to the power rating of commercially available pumps
            popGen[i]->pumps[j]->ucost = pumpdata[randomPumpStock - 1].pumpUcost;                                       //Equating the unit cost of pipe i to the unit cost of commercially available pumps
            popGen[i]->pumps[j]->realcost = popGen[i]->pumps[j]->ucost;                                             //Equating the unit cost of pipe i to the unit cost of commercially available pumps to the real cost of pump
            strcpy(popGen[i]->pumps[j]->epanetID, epanetPumpID[j]);                                                     //Copies the EPANET Pump ID to the pump attribute of pump structure
            pumpTotalRealCost += popGen[i]->pumps[j]->realcost;                                                         //Adds all the total real cost of the pumps for all the pumps in the input file
        }
        popGen[i]->pipeRealCost = pipeTotalRealCost;                                                                      //Equation the Pipe Total Real Cost to the chromosome attribute
        popGen[i]->pumpRealCost = pumpTotalRealCost;                                                                      //Equation the Pump Total Real Cost to the chromosome attribute
        popGen[i]->pipeTotalCost = penaltycost(popGen[i]);                                                                //Using the Penalty Function to get the the Pipe Total (Real+Penalty) Cost and equating it to the chromosome attribute
        popGen[i]->totcost = popGen[i]->pipeTotalCost + popGen[i]->pumpRealCost;                                      //Solving for the total cost of pipes and pumps and equating it to the chromosome attribute
        popGen[i]->fitness = fit(popGen[i], 0);                                                                           //Using the Fitness function to get the fitness of the chromosome and equating it to the chromosome attribute
    }
}

/*Evaluates the relative cost of a generated solution, subject to penalty costs*/
double penaltycost(chrom* pop) {
    double Cost = 0;                                    //Real + Penalty Cost of the Pipes in a chromosome
    double vpenalty, ppenalty1, ppenalty2;              //Penalty Computations for the Pipe Total Cost
    double pressure1, pressure2, velocity;              //Pressure at the start and end nodes of the pipe and the velocity of water in the pipe
    int n, m, linktype;                                 //Dummy variables to store the count of pipes and pumps, and type of link (1=pipe & 2=pump), respectively
    int nodetype[nNodes];                               //Array to store the Node Type of Node i (0=junction, 1=reservoir, 2=tank)
    double nodePressure[nNodes];                        //Array to store the Pressure for each node
    long time;                                          //Current Simulation Time in seconds

    /*Determine the link type from pipes and pumps and use EPANET API Functions to get information from the EPANET input file*/
    n = 0; m = 0;
    for (int i = 1; i <= nLinks; i++) {
        EN_getlinktype(ph, i, &linktype);
        if (linktype == 1) {
            EN_setlinkvalue(ph, i, EN_DIAMETER, pop->pipes[n]->diam);
            EN_setlinkvalue(ph, i, EN_ROUGHNESS, pop->pipes[n]->rough);
            n++;
        }
        if (linktype == 2) {
            EN_setlinkvalue(ph, i, EN_PUMP_POWER, pop->pumps[m]->power);
            m++;
        }
    }

    EN_openH(ph);                                       //EPANET Function to open the project hydraulic solver
    EN_initH(ph, EN_NOSAVE);                            //EPANET Function to initialize hydraulic analysis based on chromosome information but not saving its result to the input file
    EN_runH(ph, &time);                                //EPANET Function to compute the hydraulic solution from the given time

    /*Determine the node type of each node and use EPANET API Functions to get pressure of each node obtained from hydraulic analysis EPANET input file*/
    for (int i = 1; i <= nNodes; i++) {
        EN_getnodetype(ph, i, &nodetype[i - 1]);
        EN_getnodevalue(ph, i, EN_PRESSURE, &nodePressure[i - 1]);
    }

    /*Determine the link type if it is a pipe and use EPANET API Functions to get the velocity of each pipe obtained from hydraulic analysis EPANET input file*/
    n = 0; m = 0;
    for (int i = 1; i <= nLinks; i++) {
        EN_getlinktype(ph, i, &linktype);
        if (linktype == 1) {
            EN_getlinkvalue(ph, i, EN_VELOCITY, &(pop->pipes[n]->velocity));
            n++;
        }
    }
    /*Equating the pressure and velocity from each node and pipes respectively to the pipe and pump attributes and solves for the Penalty Points to get the Total Cost of the each pipe*/
    for (int i = 0; i < nPipes; i++) {
        pop->pipes[i]->ends[0]->pressure = nodePressure[pop->pipes[i]->ends[0]->nodeIndex - 1];     //Note EPANET Index Number is based 1 thus "-1" in array
        pop->pipes[i]->ends[1]->pressure = nodePressure[pop->pipes[i]->ends[1]->nodeIndex - 1];
        velocity = pop->pipes[i]->velocity;
        pressure1 = pop->pipes[i]->ends[0]->pressure;
        pressure2 = pop->pipes[i]->ends[1]->pressure;
        if (nodetype[pop->pipes[i]->ends[0]->nodeIndex - 1] != 0) {                                       //If nodetype is not Junction (!=0), then there would be no penalty
            ppenalty1 = 0;
        }
        else if (nodetype[pop->pipes[i]->ends[0]->nodeIndex - 1] == 0) {                                //Solving for Pressure Penalties of the Left Node
            ppenalty1 = pow(((2 * pressure1) - pmax - pmin), 2) / pow((pmax - pmin), 2);
            if (pressure1 >= pmin && pressure1 <= pmax) {
                ppenalty1 = PP1 * ppenalty1;
            }
            else {
                ppenalty1 = ppenalty1 * PP2;
            }
        }
        pop->pipes[i]->ends[0]->ppenalty = ppenalty1;
        if (nodetype[pop->pipes[i]->ends[1]->nodeIndex - 1] != 0) {                                       //Solving for Pressure Penalties of the Right Node
            ppenalty2 = 0;
        }
        else if (nodetype[pop->pipes[i]->ends[1]->nodeIndex - 1] == 0) {
            ppenalty2 = pow(((2 * pressure2) - pmax - pmin), 2) / pow((pmax - pmin), 2);
            if (pressure2 >= pmin && pressure2 <= pmax) {
                ppenalty2 = PP1 * ppenalty2;
            }
            else {
                ppenalty2 = ppenalty2 * PP2;
            }
        }
        pop->pipes[i]->ends[1]->ppenalty = ppenalty2;
        if (velocity >= vmin && velocity <= vmax) {                                                             //Solving for Velocity Penalties of the Pipe
            vpenalty = 1 + (VP1 * pow(((2 * velocity) - vmax - vmin), 2) / pow((vmax - vmin), 2));
        }
        else {
            vpenalty = 1 + (VP2 * pow(((2 * velocity) - vmax - vmin), 2) / pow((vmax - vmin), 2));
        }
        pop->pipes[i]->vpenalty = vpenalty;
        pop->pipes[i]->totcost = pop->pipes[i]->realcost * (1 + pop->pipes[i]->ends[0]->ppenalty + pop->pipes[i]->ends[1]->ppenalty) * pop->pipes[i]->vpenalty;
        Cost += pop->pipes[i]->totcost;
    }
    EN_saveH(ph);                   //EPANET Function to save hydraulic file to its binary output file
    EN_closeH(ph);                  //EPANET Function that closes the hydraulic solver
    return Cost;
}

/*Evaluates the fitness of a generated solution using the proposed scaled fitness*/
double fit(chrom* pop, int iter) {
    double fitness, n;
    if (iter == 1) {
        n = 0.5;
    }
    else {
        n = 0.5 + (2.5 * iter / iterations);
    }
    fitness = pow((1 / pop->totcost), n);
    return fitness;
}

/*Ranks the fitness of each chromosome*/
void fitnessRanking(chrom** pop, chrom** dummy) {
    for (int i = 1; i < gensize; i++) {
        for (int j = 0; j < (gensize - 1); j++) {
            if (pop[j]->fitness < pop[j + 1]->fitness) {
                dummy[j] = pop[j];
                pop[j] = pop[j + 1];
                pop[j + 1] = dummy[j];
            }
        }
    }
}

/*Ranks the GA Trials by decreasing cost*/
void finalRanking(chrom** pop, chrom** dummy, int limit) {
    printf("\nRanking the Best Solutions...\n");
    for (int i = 1; i < limit; i++) {
        for (int j = 0; j < (limit - 1); j++) {
            if (pop[j]->totcost > pop[j + 1]->totcost) {
                dummy[j] = pop[j];
                pop[j] = pop[j + 1];
                pop[j + 1] = dummy[j];
            }
        }
    }
}

/*Roulette wheel selection using the proportionate fitness method. Selects a candidate parent for crossover. Gives chromosome number.*/
int roulette(chrom** popGen) {
    int location;                                                                   //Location of the parent chosen in crossover function
    double prob[gensize - nDiscarded], cum_prob[gensize - nDiscarded];              //Array of the the Probability and Cumulative of each chromosome parent in the parent Pool in a generation
    double tot_fit, random;                                                         //Variable to get the total fitness and random generator for roulette spin
    tot_fit = 0;

    for (int i = 0; i < gensize - nDiscarded; i++) {
        tot_fit += popGen[i]->fitness;
    }
    for (int i = 0; i < gensize - nDiscarded; i++) {
        prob[i] = popGen[i]->fitness / tot_fit;
    }
    cum_prob[0] = prob[0];
    for (int i = 1; i < gensize - nDiscarded; i++) {
        cum_prob[i] = 0;
        for (int j = 0; j <= i; j++) {
            cum_prob[i] += prob[j];
        }
    }
    random = (double)rand() / RAND_MAX;
    for (int i = 0; i < gensize - nDiscarded; i++) {
        if (i != gensize - nDiscarded - 1) {
            if (random < cum_prob[i]) {
                location = i;
                break;
            }
        }
        else {
            location = i;
            break;
        }
    }
    return location;
}

/*Generates the new population by Uniform Crossover + Mutation Operator*/
void crossover(chrom** currentGen, chrom** nextGen, int iter) {
    int location1, location2;                                                       //Parent 1 and Parent 2 Chromosomes
    int random, dummy;                                                              //Coinflip randomizer and dummy variable for pipe and pump offsprings
    int pipeParent0[nPipes], pipeParent1[nPipes];                                   //Temporary storage of pipe chromosome selected as parents
    int pumpParent0[nPumps], pumpParent1[nPumps];                                   //Temporary storage of pump chromosome selected as parents
    int pumpOffspring[nPumps], pipeOffspring[nPipes];                               //Temporary storage of pipe and pump chromosomes selected as offsprings, respectively
    double mutationRandom;                                                          //Mutation randomizer variable
    int randomPipeStock, randomPumpStock;                                           //Index variable for the chosen parents on pipe and pump chromosomes
    int pipeCost = 0, pumpCost = 0;                                                 //Pipe and Pump Cost placeholder, respectively

    for (int i = 0; i < gensize; i++) {
        if (i < nElite) {
            nextGen[i]->fitness = fit(nextGen[i], iter);
        }
        if (i >= nElite) {
            do {                                                                    //Parent Selection does not allow parent 1 = parent 2 (Asexual Reproduction) because it will converge prematurely
                location1 = roulette(currentGen);
                location2 = roulette(currentGen);
            } while (location1 == location2);
            for (int j = 0; j < nPipes; j++) {
                pipeParent0[j] = currentGen[location1]->pipes[j]->pipeID;
                pipeParent1[j] = currentGen[location2]->pipes[j]->pipeID;
            }
            for (int j = 0; j < nPumps; j++) {
                pumpParent0[j] = currentGen[location1]->pumps[j]->pumpID;
                pumpParent1[j] = currentGen[location2]->pumps[j]->pumpID;
            }
            for (int j = 0; j < nPipes; j++) {
                random = (int)rand() % 2;
                if (random == 0) {
                    pipeOffspring[j] = pipeParent0[j];
                }
                else if (random == 1) {
                    pipeOffspring[j] = pipeParent1[j];
                }
                /*Pipe Mutation*/
                mutationRandom = rand() % 101;
                if (mutationRandom <= mrate) {
                    dummy = pipeOffspring[j];
                    do {
                        randomPipeStock = (rand() % (pipeStock)) + 1;
                        pipeOffspring[j] = pipedata[randomPipeStock - 1].pipeID;
                    } while (dummy == pipeOffspring[j]);
                }
                nextGen[i]->pipes[j]->pipeID = pipedata[pipeOffspring[j] - 1].pipeID;
                nextGen[i]->pipes[j]->diam = pipedata[pipeOffspring[j] - 1].pipeDiam;
                nextGen[i]->pipes[j]->ucost = pipedata[pipeOffspring[j] - 1].pipeUcost;
                nextGen[i]->pipes[j]->rough = pipedata[pipeOffspring[j] - 1].pipeRough;
                nextGen[i]->pipes[j]->realcost = nextGen[i]->pipes[j]->ucost * nextGen[i]->pipes[j]->length;
                pipeCost += nextGen[i]->pipes[j]->realcost;
            }
            for (int j = 0; j < nPumps; j++) {
                random = (int)rand() % 2;
                if (random == 0) {
                    pumpOffspring[j] = pumpParent0[j];
                }
                else if (random == 1) {
                    pumpOffspring[j] = pumpParent1[j];
                }
                /*Pump Mutation*/
                mutationRandom = rand() % 101;
                if (mutationRandom <= mrate) {
                    do {
                        dummy = pumpOffspring[j];
                        randomPumpStock = (rand() % (pumpStock)) + 1;
                        pumpOffspring[j] = pumpdata[randomPumpStock - 1].pumpID;
                    } while (dummy == pumpOffspring[j]);
                }
                nextGen[i]->pumps[j]->pumpID = pumpdata[pumpOffspring[j] - 1].pumpID;
                nextGen[i]->pumps[j]->power = pumpdata[pumpOffspring[j] - 1].pumpPower;
                nextGen[i]->pumps[j]->ucost = pumpdata[pumpOffspring[j] - 1].pumpUcost;
                nextGen[i]->pumps[j]->realcost = nextGen[i]->pumps[j]->ucost;
                pumpCost += nextGen[i]->pumps[j]->realcost;
            }
            nextGen[i]->pipeRealCost = pipeCost;
            nextGen[i]->pumpRealCost = pumpCost;
            nextGen[i]->pipeTotalCost = penaltycost(nextGen[i]);
            nextGen[i]->totcost = nextGen[i]->pipeTotalCost + nextGen[i]->pumpRealCost;
            nextGen[i]->fitness = fit(nextGen[i], iter);
            pipeCost = 0;
            pumpCost = 0;
        }
    }
}