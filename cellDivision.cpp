#define _GLIBCXX_USE_CXX11_ABI 0
#include "std_include.h"
#include "cuda_runtime.h"
#include "cuda_profiler_api.h"

#define ENABLE_CUDA

#include "vertexQuadraticEnergy.h"
#include "noiseSource.h"
#include "voronoiQuadraticEnergyWithTension.h"
#include "selfPropelledCellVertexDynamics.h"
#include "brownianParticleDynamics.h"
#include "DatabaseNetCDFAVM.h"
#include "DatabaseNetCDFSPV.h"
#include "DatabaseTextVoronoi.h"
#include "gpubox.h"

/*!
 This file performs the simulation of a finite number of cells in a square box
 with periodic boundary conditions.
 There is a growth rate whose form form can be chosen by indexrate.
 Birth and death events are stochastic and modeled using a Gillespie algorithm.
*/
int main(int argc, char*argv[])
{
    //There are a bunch of default parameters that
    //can be changed from the command line
    int numpts = 400;
    int USE_GPU = 0;
    int USE_TENSION = 0;
    int c;
    int tSteps = 1000000;
    int initSteps = 100;

    Dscalar dt = 0.01;
    Dscalar p0 = 3.64;
    Dscalar p0_0 = 3.65; // Resident
    Dscalar p0_1 = 3.90; // Mutant
    Dscalar a0 = 1.0;
    Dscalar v0 = 0.05;
    Dscalar Dv0 = 0;
    Dscalar Dr = 1.0;
    Dscalar DrM = 1.0; //Rotational diffusion of the mutant
    Dscalar gamma = 0.0;
    Dscalar init_fraction = 0.5;

    int random_param = 0; // Determines if parameters are drawn randomly
    Dscalar KAA = 1;
    Dscalar KPP = 1;
    Dscalar KAA_0 = 1;
    Dscalar KPP_0 = 1;
    Dscalar KAA_1 = 1;
    Dscalar KPP_1 = 1;
    Dscalar alpha = 1; //area preference of type 2
    Dscalar lambda = 1; // scaling factor for the mechanics coupling of resident
    Dscalar lambdaM = lambda; // scaling factor for the mechanics coupling of mutant
    Dscalar lambda_1 = 1; //Coupling between division rate and area
    Dscalar lambda_2 = KPP_0*p0_0/2/a0; //Coupling between division rate and perimeter

    //program_switch plays an awkward role in this example of both selecting vertex vs Voronoi model,
    //and also determining whether to save output files... read below for details
    int program_switch = +2;
    //division sets if there is division and death if division!=0.
    int division = 1;
    //pressure_dep determines wether division is complitely random (0 in
    //that case) or not (any non zero value)
    int pressure_dep = 1;
    //init_type declares the type of initial mixture: 
    //0 for random, 1 for half plane, 2 for circular colony
    int init_type = 0;
    //Id of the job, useful to name the file
    int id = 0;
    int gillespie = 1;
    //Choice of rate used in gillespie
    int indexrate = 1;
    //Timescale for division rate
    Dscalar divrate = 0.05/numpts;
    //Turnover rate for cells;
    Dscalar r0=1;
    Dscalar r0M=1;//Value for the mutants
    //Just count total number of divisions
    int divCount = 0;
    //Controling a fixed amount of cells after initialization. Hence non
    //normalized number of cells
    int Ntarget = numpts; // By default it is the same as initial condition
    // Introduce selective advantage
    Dscalar Ph_adv = 0; // Selective advantage on the homeostatic pressure
    Dscalar Ph = 0;
    while((c=getopt(argc,argv,"n:g:m:s:r:a:i:v:b:x:z:o:p:u:y:t:e:d:k:l:f:q:h:w:j:")) != -1)
        switch(c)
        {
            case 'n': numpts = atoi(optarg); break;
            case 't': tSteps = atoi(optarg); break;
            case 'g': KAA_1 = atof(optarg); break;
            case 'x': KPP_1 = atof(optarg); break;
            case 'i': indexrate = atoi(optarg); break;
            case 'o': pressure_dep = atoi(optarg); break;
            case 'z': r0M = atof(optarg); break;
            case 'e': dt = atof(optarg); break;
            case 'b': lambdaM = atof(optarg); break;
            case 's': Ph_adv = atof(optarg); break;
            case 'p': a0 = atof(optarg); break;
            case 'u': p0_0 = atof(optarg); break;
            case 'y': p0_1 = atof(optarg); break;
            case 'a': alpha = atof(optarg); break; 
            case 'v': v0 = atof(optarg); break;
            case 'd': DrM = atof(optarg); break;
            case 'k': division = atof(optarg); break;
            case 'l': init_type = atof(optarg); break;
            case 'f': init_fraction = atof(optarg); break;
            case 'q': Ph = atof(optarg); break;
            case 'h': Dv0 = atof(optarg); break;
            case 'w': id = atoi(optarg); break;
            case 'm': Ntarget = atoi(optarg); break;
            case 'j': random_param = atoi(optarg); break;
            case '?':
                    if(optopt=='c')
                        std::cerr<<"Option -" << optopt << "requires an argument.\n";
                    else if(isprint(optopt))
                        std::cerr<<"Unknown option '-" << optopt << "'.\n";
                    else
                        std::cerr << "Unknown option character.\n";
                    return 1;
            default:
                       abort();
        };
    clock_t t1,t2;

    bool reproducible = false;
    bool initializeGPU = true;
    //shear variable to decide if we do a shear or not
    bool shear = false;
    //if (USE_GPU >= 0)
    //    {
    //    bool gpu = chooseGPU(USE_GPU);
    //    if (!gpu) return 0;
    //    cudaSetDevice(USE_GPU);
    //    }
    //else
    //    initializeGPU = false;

    char dataname[256];
    sprintf(dataname,"../results/test1.nc");

    //Deciding which cell types gets the fitness advantage or
    //disadvantage
    Dscalar adv = 0;
    Dscalar adv0 = 0;
    Dscalar adv1 = Ph_adv;

    //Might be useful sometimes
    //r0=r0M

    Dscalar Ph_1 = Ph; 

    Dscalar gamma_1 = 1; //reference value
    Dscalar gamma_2 = gamma; // value inputed at the prompt


    //noiseSources are random number generators used within the code... here we'll randomly select
    //cells to divide
    noiseSource noise;
    noise.Reproducible = reproducible;

    if (random_param == 1)
    {
        //Set up to compute the homeostatic pressure from the index
        //of the file to go through value of parameters
        //
        int numsample = 10;
        int exp = (int) (id-1)/numsample;
        int expInd = (id-1)%numsample;
        switch(exp)
        {
            case 0:
                Ph_adv = -0.1;
                //Ph_1 = noise.getRealUniform(-1.8,-1.0);
                //KPP_1 = 0.5;
                //Ph_1 = noise.getRealUniform(-0.8,-0.2);
                //Ph_2 = Ph_1;
                //Ph_adv = Ph_1; // This is for bookkeeping purpuses, it doesn't actually do anything, just for the output file to know the pressure
                //p0_1 = noise.getRealUniform(3.5,5.66);
                //alpha = noise.getRealUniform(0.5,1.0);
                //v0 = noise.getRealUniform(0,2);
                //KPP_1 = 0.1+expInd*10/numsample;
                break;
            case 1:
                alpha = 1.1;
                //Ph_adv = -0.29;
                //a0 = 1.1;
                //Ph_1 = noise.getRealUniform(-1.8,-1.0);
                //Ph_2 = Ph_1;
                //Ph_adv = Ph_1; // This is for bookkeeping purpuses, it doesn't actually do anything, just for the output file to know the pressure
                //KPP_1 = noise.getRealUniform(0.1,10);
                //Dv0 = expInd*1/numsample;
                break;
            case 2:
                p0_1 = 3.86;
                //Ph_adv = -0.57;
                //KAA_1 = noise.getRealUniform(0.1,10);
                //DrM = 0.1+expInd*10/numsample;
                break;
            case 3:
                KPP_1 = 0.5;
                //Ph_adv = 0.5;
                //v0 = noise.getRealUniform(0,2);
                //p0_1 = 3+expInd*1.4/numsample;
                break;
            case 4:
                KAA_1 = 0.5;
                //DrM = noise.getRealUniform(0,2);
                //alpha = 0.5+expInd*1.5/numsample;
                break;
            case 5:
                lambdaM = 0.5;
                //KAA_1 = 0.1+expInd*10/numsample;
                break;
            case 6:
                Dv0 = -0.02;
            case 7:
                DrM = 0.5; 

        }

        //Randomly determine value of parameters

        //p0_1 = noise.getRealUniform(3,4.4);
        //alpha = noise.getRealUniform(0.5,2);
        //Ph_1 = noise.getRealUniform(-1.8,-1.2);
        //Ph_2 = Ph_1;
        //Ph_adv = Ph_1; // This is for bookkeeping purpuses, it doesn't actually do anything, just for the output file to know the pressure
        //KAA_1 = noise.getRealUniform(0.1,4);
        //KPP_1 = noise.getRealUniform(0.1,4);
        //lambda_1 = noise.getRealUniform(0.2,4);
        //lambda_2 = noise.getRealUniform(0.2*KPP_0*p0_0/2/a0,4*KPP_0*p0_0/2/a0);
    }

    // Defines the homeostatic pressure of type 2 after specifying
    // the change in pressure
    Dscalar Ph_2 = Ph+Ph_adv;

    int numpts_normalized = numpts;
    //Dscalar boxScale = sqrt(alpha);
    Dscalar boxScale = 1;

    //program_switch >= 0 --> self-propelled voronoi model
    if(program_switch >=0)
        {
        EOMPtr spp = make_shared<selfPropelledParticleDynamics>(numpts_normalized);
        //shared_ptr<VoronoiQuadraticEnergyWithTension> spv = make_shared<VoronoiQuadraticEnergyWithTension>(numpts,a0,p0,reproducible);
        shared_ptr<VoronoiQuadraticEnergy> spv = make_shared<VoronoiQuadraticEnergy>(numpts_normalized,boxScale,a0,p0,reproducible);

        //Changing box size to compress or strech
        //Dscalar stretch = 1.2;
        //BoxPtr newbox = make_shared<gpubox>(stretch*sqrt(numpts),stretch*sqrt(numpts));
        //spv->setBox(newbox);

        //Changing aspect ratio of the box
        //Dscalar aspectRatio = 1.5;
        //BoxPtr newbox = make_shared<gpubox>(aspectRatio*sqrt(numpts),sqrt(numpts)/aspectRatio);
        //spv->setBox(newbox);
        
        //for variety, we'll have cell division between two types of cells, with some applied surface tension between the types
        //...this section does the usual business of setting up the simulation
        vector<int> types(numpts_normalized);
        /*
        // Set initially all cell from different lineage
        vector<int> lineages(numpts);
        for (int ii = 0; ii < numpts; ++ii)
        {
            lineages[ii] = ii;
        }
        */
        //
        
        // Vector that stores parameters for each cell, that
        // is determined in the loop that set cell types;
        vector<Dscalar2> AP(numpts_normalized); // A0 and P0
        vector<Dscalar> NewKA(numpts_normalized); // Ks
        vector<Dscalar> NewKP(numpts_normalized); //Gammas
        vector<Dscalar> V0s(numpts_normalized); //Velocities
        vector<Dscalar> Drs(numpts_normalized); // Rotational diffusion
        Dscalar boxL = (Dscalar) sqrt(numpts_normalized)*boxScale;

        //Choose at random the mechanical parameters according to
        //random_param value
        //Putting the name after having the variables
        char datatextname[256];
        sprintf(datatextname,"./%dCells_%dSteps_initial-%d_div-%d_pressdep-%d_%d.txt",numpts_normalized,tSteps,init_type,division,pressure_dep,id);
        //netCDF databases require the same number of cells in every frame... the text databases lift that limitation, so are useful here
        //DatabaseTextVoronoi db1("../results/test_div_mixture.txt",0);
        DatabaseTextVoronoi db1(datatextname,0);
        // Random mixture of types 1 and 2
        if (init_type == 0)
        {
            for (int ii = 0; ii < numpts_normalized; ++ii)
            {
                if(ii < int(numpts_normalized*init_fraction))
                {
                    types[ii]=1;
                }
                else
                {
                    types[ii]=0;
                }
            }
            random_shuffle(types.begin(),types.end());
            for (int ii = 0; ii < numpts_normalized; ++ii)
            {
                //types[ii]=ii; // Label each initial cell to check lineage
                //types[ii]=noise.getInt(0,1); //Label randomly initially with 0 or 1 for a mixture 
                //if(types[ii]==0)
                if(types[ii]==0)
                // Solid-like cells, value taken from Sussman paper on
                // demixing
                {
                    AP[ii].x=a0;
                    AP[ii].y=p0_0;
                    NewKA[ii]=KAA_0;
                    NewKP[ii]=KPP_0;
                    V0s[ii] = v0;
                    Drs[ii] = Dr;
                }
                //else if(types[ii]==1)
                else 
                // Liquid-like cells, value taken from Sussman paper on
                // demixing
                {
                    AP[ii].x=alpha;
                    AP[ii].y=p0_1;
                    NewKA[ii]=KAA_1;
                    NewKP[ii]=KPP_1;
                    V0s[ii] = v0 + Dv0;
                    Drs[ii] = DrM;
                }
            }
        }
        // Left/Right split of type 1 and 2, resp.
        else if (init_type ==1)
        {
            ArrayHandle<Dscalar2> h_cp(spv->cellPositions);
            for (int ii = 0; ii < numpts_normalized; ++ii)
            {
                if (h_cp.data[ii].x < init_fraction*boxL)
                {
                    types[ii] = 1;
                    AP[ii].x=alpha;
                    AP[ii].y=p0_1;
                    NewKA[ii]=KAA_1;
                    NewKP[ii]=KPP_1;
                    V0s[ii] = v0;
                    Drs[ii] = DrM;
                }
                else 
                {
                    types[ii] = 0;
                    AP[ii].x=a0;
                    AP[ii].y=p0_0;
                    NewKA[ii]=KAA_0;
                    NewKP[ii]=KPP_0;
                    V0s[ii] = v0 + Dv0;
                    Drs[ii] = Dr;
                }
            }
        }
        // Central circular drop of type 2
        else if (init_type ==2)
        {
            ArrayHandle<Dscalar2> h_cp(spv->cellPositions);
            Dscalar radius = sqrt(init_fraction/PI)*boxL; //Choose size of circle as fraction of box size
            vector<Dscalar> center(2);
            center[0] = boxL/2; center[1] = boxL/2; // Center of the box
            for (int ii = 0; ii < numpts_normalized; ++ii)
            {
                if (pow(h_cp.data[ii].x-center[0],2) + pow(h_cp.data[ii].y-center[1],2) < pow(radius,2))
                {
                    types[ii] = 1;
                    AP[ii].x=alpha;
                    AP[ii].y=p0_1;
                    NewKA[ii]=KAA_1;
                    NewKP[ii]=KPP_1;
                    V0s[ii] = v0;
                    Drs[ii] = DrM;
                }
                else
                {
                    types[ii] = 0;
                    AP[ii].x=a0;
                    AP[ii].y=p0_0;
                    NewKA[ii]=KAA_0;
                    NewKP[ii]=KPP_0;
                    V0s[ii] = v0+ Dv0;
                    Drs[ii] = Dr;
                }
            }
        }
        spv->setCellType(types);
        //spv->setCellLineage(lineages);
        // No surface tension for now
        //spv->setSurfaceTension(gamma);
        //spv->setUseSurfaceTension(true);

        //Set cell preferences depending on the type
        spv->setCellPreferences(AP);
        //Set mechanical parameters depending on type
        spv->setModuli(NewKA,NewKP);
        //spv->setCellPreferencesUniform(1.0,p0);

        //Change velocity preferences to not be uniform
        //spv->setv0Dr(v0,1.0);
        spv->setCellMotility(V0s, Drs);

        //Change mechanical moduli
        //spv->setModuliUniform(KAA,KPP);

        SimulationPtr sim = make_shared<Simulation>();
        sim->setConfiguration(spv);
        sim->addUpdater(spp,spv);
        sim->setIntegrationTimestep(dt);
        //sim->setSortPeriod(initSteps/10);
        sim->setCPUOperation(!initializeGPU);
        sim->setReproducible(reproducible);

        Dscalar RunningAverage = 0; //Variable to have a running average at each time step
        vector<Dscalar> cellPressureInit;
        //perform some initialization timesteps
        //if program_switch = 2, save output file every so often
        for (int timestep = 0; timestep < initSteps+1; ++timestep)
            {
            sim->performTimestep();
            cellPressureInit = spv->computePressure();
            if(program_switch == 2 && timestep%((int)(1/dt))==0)
                {
                cout << timestep << endl;
                db1.WriteState(spv,p0_0,p0_1,KPP_0,KPP_1,alpha,KAA_0,KAA_1,Ph_1,lambda,divrate,Ph_adv,v0,DrM,init_fraction,RunningAverage, cellPressureInit);
                };
            };

        //to have a Voronoi model division, one needs to pass in various parameters.
        //the integer vector (of length 1) indicates the cell to divide
        //the Dscalar vector (of length 2) parameterizes the geometry of the cell division... see voronoiModelBase for details
        vector<int> cdtest(1); cdtest[0]=10;
        vector<Dscalar> dParams(2); dParams[0] = 3.0*PI/4.0-.1; dParams[1] = 0.5;
        int Ncells = spv->getNumberOfDegreesOfFreedom();

        //in this example, divide a cell every divisionTime tau
        int divisionTime = 20;
        Dscalar saveTime = 1000*dt; // Time for saving files in real time units, here almost every generation 
        Dscalar DivCut = 0.8; //On average division every 10/0.8 steps;
        Dscalar DeathCut = 0.8; //On average death every 10/0.8 steps;
        Dscalar divtime;
        Dscalar deathtime;
        //Dscalar divrate = 0.01/(2);
        Dscalar scale = 1; //Parameter for the sigmoid dependent rates.
        int relaxtime = 2; // Value for relaxation timesteps between events.
        Dscalar tau; // Next event time in Gillespie
        int timer = 0; //timer for Gillespie simulation
        // reaction is reaction identifier. Death of cell j if = 2j, and
        // division of cell j if = 2j+1
        vector<int> reaction(1); reaction[0]=0;
        //r1 and r2 random numbers for Gillespie.
        Dscalar r1;
        Dscalar r2;
        Dscalar Sa;
        //aa0 rate of having any event.
        Dscalar aa0b=0; // Define the sum for birth event
        Dscalar aa0d=0; // Define the sum for death event to normalize
        Dscalar aa0=0; // Define the sum for all events
        // Place holders for some debugging
        int output1;
        Dscalar output2;
        Dscalar cellarea; //Just a variable to avoid to many calls
        Dscalar cellperi; //Just a variable to avoid to many calls
        Dscalar areapref; //Just a variable to avoid to many calls
        Dscalar peripref; //Just a variable to avoid to many calls
        Dscalar rate; // Just a variable to avoid to many calls
        Dscalar pressure; //Pressure variable to compute the rate
        Dscalar K_cell; // Just a variable to avoid to many calls
        Dscalar G_cell; // Just a variable to avoid to many calls
        int type_temp; // A variable to avoid many calls

        //Displacements to perform a shear
        GPUArray<Dscalar2> displacements;
        Dscalar m = 0.1; //scale for the shear
        t1=clock();
        Dscalar TempPos = 0; // temp variable for the shear
        Dscalar AverageNc = 0; //Variable to compute the average cell number over the time of simulation
        //ArrayHandle<Dscalar2> h_area(spv->returnAreaPeri(),access_location::host,access_mode::read);
        for (int timestep = 0; timestep < tSteps; ++timestep)
        {
            sim->performTimestep();
            //Call to compute the pressures in each cell at
            //every time to allow saving
            vector<Dscalar> cellPressure = spv->computePressure();
            // Compute pressure using the global value. DO NOT USE
            //vector<Dscalar> cellPressureBis(Ncells);
            //ArrayHandle<Dscalar2> h_area(spv->returnAreaPeri(),access_location::host,access_mode::read);
            //ArrayHandle<Dscalar2> h_pref(spv->returnAreaPeriPreferences(),access_location::host,access_mode::read);
            //ArrayHandle<Dscalar2> h_moduli(spv->returnModuli(),access_location::host,access_mode::read);
            //Dscalar PressureBisSum=0;
            //Dscalar PressureSum=0;
            //for (int jj = 0; jj<Ncells; ++jj)
            //{
            //    cellarea = h_area.data[jj].x;
            //    cellperi = h_area.data[jj].y;
            //    areapref = h_pref.data[jj].x;
            //    peripref = h_pref.data[jj].y;
            //    K_cell = h_moduli.data[jj].x;
            //    G_cell = h_moduli.data[jj].y;
            //    cellPressureBis[jj] = -(K_cell*(cellarea-areapref) + G_cell*(cellperi-peripref)*cellperi/2/cellarea);
            //    PressureBisSum+=cellarea*cellPressureBis[jj];
            //    PressureSum+=cellarea*cellPressure[jj];
            //    cout << cellPressureBis[jj] << " " <<cellPressure[jj] << endl;
            //}
            //cout << "Force Sum: " << PressureBisSum << " " << PressureSum << endl; 
            //
            //
            //vector<Dscalar> cellPeri = spv->computePeri();
            if(shear)
            {
                displacements.resize(Ncells);
                //Performing shear
                ArrayHandle<Dscalar2> h_cellpos(spv->cellPositions,access_location::host,access_mode::read);
                ArrayHandle<Dscalar2> h_disp(displacements,access_location::host,access_mode::readwrite);
                for (int jj = 0; jj<Ncells; ++jj)
                {
                   TempPos = h_cellpos.data[jj].x;
                   //Shear the left half of the box
                   h_disp.data[jj].x = 0; 
                   h_disp.data[jj].y = m*(1-2/boxL*abs(TempPos-boxL/2)); 
                   //Left shear
                   //
                   //if(TempPos < boxL/2)
                   //{
                   //    h_disp.data[jj].x = 0;   
                   //    h_disp.data[jj].y = m;
                   //}
                   //else
                   //{
                   //    h_disp.data[jj].x = 0;
                   //    h_disp.data[jj].y = 0;
                   //}
                   //
                   //Simple shear, does not work
                   //
                   //h_disp.data[jj].x = 0;   
                   //h_disp.data[jj].y = m*(TempPos-boxL/2);   
                }
                spv->moveDegreesOfFreedom(displacements, 1.0);
            }
            if(division!=0)
            {
                if(program_switch >0 && gillespie!=0)
                {
                    if(timer==0)
                    {
                    //When timer is 0, initiate when and what is the next
                    //reaction
                        ArrayHandle<Dscalar2> h_area(spv->returnAreaPeri(),access_location::host,access_mode::read);
                        ArrayHandle<int> h_type(spv->cellType,access_location::host,access_mode::read);
                        ArrayHandle<Dscalar2> h_pref(spv->returnAreaPeriPreferences(),access_location::host,access_mode::read);
                        ArrayHandle<Dscalar2> h_moduli(spv->returnModuli(),access_location::host,access_mode::read);
                        r1 = noise.getRealUniform(0,1);
                        r2 = noise.getRealUniform(0,1);
                        vector<Dscalar> aa(2*Ncells);
                        //aa0 rate of having any event.
                        aa0b=0; // Define the sum for birth event
                        aa0d=0; // Define the sum for death event to normalize
                        aa0=0; // Define the sum for all events
                        for (int jj = 0; jj<Ncells; ++jj)
                        {
                            cellarea = h_area.data[jj].x;
                            cellperi = h_area.data[jj].y;
                            type_temp = h_type.data[jj];
                            areapref = h_pref.data[jj].x;
                            peripref = h_pref.data[jj].y;
                            K_cell = h_moduli.data[jj].x;
                            G_cell = h_moduli.data[jj].y;
                            pressure = -(K_cell*(cellarea-areapref) + G_cell*(cellperi-peripref)*cellperi/2/cellarea);
                            if(pressure_dep==0)
                            {
                                //Even number events are death
                                aa[2*jj] = divrate;
                                aa0d = aa0d + aa[2*jj];
                                //Odd number events are births
                                aa[2*jj+1] = divrate;
                                aa0b = aa0b + aa[2*jj+1];
                            }
                            else
                            {
                                //Even number events are death
                                //aa[2*jj] = exp(-h_area.data[jj].x/(lambda*a0));
                                //aa[2*jj] = divrate*exp(-(cellarea-a0)/(lambda*a0)); // Exponential rate
                                //aa[2*jj] = scale*divrate*(1-0.5*(1+lambda*(cellarea-a0)/sqrt(1+pow(lambda*(cellarea-a0),2)))); // Logistic rate
                                //aa[2*jj] = divrate*(1-0.5*(1+erf(lambda*(cellarea-0.5*areapref)))+(adv0+adv1)/2); // Error function rate 
                                switch(indexrate)
                                {
                                    case 0:
                                        rate = 1-erf(lambda*(cellarea-areapref))+(adv0+adv1)/2;
                                        break;
                                    case 1:
                                        rate = exp(-cellarea+areapref)+(adv0+adv1)/2;
                                        break;
                                    case 2:
                                        if(type_temp==0)
                                        {
                                            if(cellarea-areapref>0){rate = r0 + (adv0+adv1)/2;}
                                            else {rate = r0 +lambda*(-cellarea+areapref)+(adv0+adv1)/2;}
                                        }
                                        else
                                        {
                                            if(cellarea-areapref>0){rate = r0M + (adv0+adv1)/2;}
                                            else {rate = r0M +lambdaM*(-cellarea+areapref)+(adv0+adv1)/2;}
                                        }
                                        break;
                                    case 3:
                                        rate = r0;
                                        break;
                                    case 4:
                                        rate = r0;
                                        break;
                                    case 5:
                                        rate = r0;
                                        break;
                                    case 6:
                                        if(type_temp==0)
                                        {
                                            if(Ph_1-pressure<0)
                                            {
                                                rate = r0+lambda*(pressure-Ph_1);
                                            }
                                            else{rate = r0;}
                                        }
                                        else
                                        {
                                            if(Ph_2-pressure<0)
                                            {
                                                rate = r0M+lambdaM*(pressure-Ph_2);
                                            }
                                            else{rate = r0M;}
                                        }
                                        break;
                                    case 7:
                                        if(type_temp==0)
                                        {
                                            rate=0;
                                        }
                                        else
                                        {
                                            if(cellarea-areapref>0){rate = r0M + (adv0+adv1)/2;}
                                            else {rate = r0M+lambdaM*(-cellarea+areapref)+(adv0+adv1)/2;}
                                        }
                                        break;
                                    case 8:
                                        if(timestep < int(tSteps/4))
                                        {
                                            rate = r0; //This number to be changed by some input
                                        }
                                        else
                                        {
                                            rate = 0;
                                        }
                                        break;
                                    case 9:
                                        rate = r0;
                                        break;
                                    case 10:
                                        if(type_temp==0)
                                        {
                                            if(Ph_1-cellPressure[jj]<0)
                                            {
                                                rate = r0+lambda*(cellPressure[jj]-Ph_1);
                                            }
                                            else{rate = r0;}
                                        }
                                        else
                                        {
                                            if(Ph_2-cellPressure[jj]<0)
                                            {
                                                rate = r0M+lambdaM*(cellPressure[jj]-Ph_2);
                                            }
                                            else{rate = r0M;}
                                        }
                                        break;
                                }
                                aa[2*jj] = divrate*rate;
                                aa0d = aa0d + aa[2*jj];

                                //Odd number events are births
                                if(type_temp==0)
                                {
                                    //aa[2*jj] = exp(-h_area.data[jj].x/(lambda*a0));
                                    //aa[2*jj+1] = divrate*exp((cellarea-a0)/(lambda*a0));
                                    //aa[2*jj+1] = scale*divrate*0.5*(1+lambda*(cellarea-a0)/sqrt(1+pow(lambda*(cellarea-a0),2))); // Logistic rate
                                    //aa[2*jj+1] = divrate*(0.5*(1+erf(lambda*(cellarea-1.5*areapref)))+adv0); // Error function rate 
                                    switch(indexrate)
                                    {
                                        case 0:
                                            rate = 1+erf(lambda*(cellarea-areapref))+adv0;
                                            break;
                                        case 1:
                                            rate = exp(cellarea-areapref)+adv0;
                                            break;
                                        case 2:
                                            if(cellarea-areapref<0){rate = r0+adv0;}
                                            else {rate = r0+lambda*(cellarea-areapref)+adv0;}
                                            break;
                                        case 3:
                                            rate = exp(lambda*(K_cell*(cellarea-areapref)+Ph_1));
                                            break;
                                        case 4:
                                            rate = exp(lambda*(Ph_1-pressure));
                                            break;
                                        case 5:
                                            rate = exp(lambda*(G_cell*(cellperi-peripref) + Ph_1));
                                            break;
                                        case 6:
                                            if(Ph_1-pressure>0)
                                            {
                                                rate = r0-lambda*(pressure-Ph_1);
                                            }
                                            else{rate = r0;}
                                            break;
                                        case 7:
                                            rate = 0;
                                            break;
                                        case 8:
                                            if(timestep < int(tSteps/4))
                                            {
                                                rate = Ntarget/Ncells; //This number to be changed by some input
                                            }
                                            else
                                            {
                                                rate = 0;
                                            }
                                            break;
                                        case 9:
                                            rate = r0 + lambda_1*(cellarea-areapref)+lambda_2*(cellperi-peripref);
                                            break;
                                        case 10:
                                            if(Ph_1-cellPressure[jj]>0)
                                            {
                                                rate = r0-lambda*(cellPressure[jj]-Ph_1);
                                            }
                                            else{rate = r0;}
                                            break;
                                    }
                                    aa[2*jj+1] = divrate*rate;
                                    //aa[2*jj+1] = 0.5*divrate; //Constant division rate for birth.
                                }
                                else
                                {
                                    //aa[2*jj] = exp(-h_area.data[jj].x/(lambda*a0));
                                    //aa[2*jj+1] = divrate*exp((cellarea-a0)/(lambda*a0));
                                    //aa[2*jj+1] = scale*divrate*0.5*(1+lambda*(cellarea-a0)/sqrt(1+pow(lambda*(cellarea-a0),2))); // Logistic rate
                                    //aa[2*jj+1] = divrate*(0.5*(1+erf(lambda*(cellarea-1.5*areapref)))+adv1); // Error function rate 
                                    switch(indexrate)
                                    {
                                        case 0:
                                            rate = 1+erf(lambdaM*(cellarea-areapref))+adv1;
                                            break;
                                        case 1:
                                            rate = exp(cellarea-areapref)+adv1;
                                            break;
                                        case 2:
                                            if(cellarea-areapref<0){rate =r0M+ adv1;}
                                            else {rate = r0M+lambdaM*(cellarea-areapref)+adv1;}
                                            break;
                                        case 3:
                                            rate = exp(lambdaM*(K_cell*(cellarea-areapref)+Ph_2));
                                            break;
                                        case 4:
                                            rate = exp(lambdaM*(Ph_2-pressure));
                                            break;
                                        case 5:
                                            rate = exp(lambdaM*(G_cell*(cellperi-peripref) + Ph_2));
                                            break;
                                        case 6:
                                            if(Ph_2-pressure>0)
                                            {
                                                rate = r0M-lambdaM*(pressure-Ph_2);
                                            }
                                            else{rate = r0M;}
                                            break;
                                        case 7:
                                            if(cellarea-areapref<0){rate =r0M+ adv1;}
                                            else {rate = r0M+lambdaM*(cellarea-areapref)+adv1;}
                                            break;
                                        case 8:
                                            if(timestep < int(tSteps/4))
                                            {
                                                rate = Ntarget/Ncells; //This number to be changed by some input
                                            }
                                            else
                                            {
                                                rate = 0;
                                            }
                                            break;
                                        case 9:
                                            rate = r0M + lambda_1*(cellarea-areapref)+lambda_2*(cellperi-peripref);
                                            break;
                                        case 10:
                                            if(Ph_2-cellPressure[jj]>0)
                                            {
                                                rate = r0M-lambdaM*(cellPressure[jj]-Ph_2);
                                            }
                                            else{rate = r0M;}
                                            break;
                                    }
                                    aa[2*jj+1] = divrate*rate;
                                    //aa[2*jj+1] = 0.5*divrate; //Constant division rate for birth.
                                }
                                aa0b = aa0b + aa[2*jj+1];
                            }
                            //aa[2*jj+1] = divrate;

                            //if(jj%2 ==0)
                            //{
                            //    //Without pressure dependence,
                            //    //deathrate = divrate
                            //    if(pressure_dep==0)
                            //    {
                            //        aa[jj] = divrate;
                            //        //aa0d = aa0d + aa[jj];
                            //    }
                            //    //Exponential distribution of area for
                            //    //pressure dependence
                            //    else
                            //    {
                            //        aa[jj] = exp(-h_area.data[int(jj/2)].x/(lambda*a0));
                            //        aa0d = aa0d + aa[jj];
                            //    }
                            //}
                            ////Odd number events are births
                            //else
                            //{
                            //    aa[jj] = divrate;
                            //    //aa0b = aa0b + aa[jj];
                            //}
                        }
                        //Normalize death rate to have homeostasis on
                        //average. 
                        //for (int jj = 0; jj<Ncells; ++jj)
                        //{
                        //    if(pressure_dep=!0)
                        //    {
                        //        aa[2*jj]=Ncells*divrate*aa[2*jj]/aa0d;
                        //    }
                        //}
                        //Get the time to next event
                        //aa0 = aa0d+Ncells*divrate; 
                        aa0 = aa0d+aa0b; 
                        //aa0 = 2*Ncells*divrate; //The way it is set up now it is predetermined
                        tau = 1/aa0*log(1/r1);
                        //Get the identity of the event
                        Sa = aa[0];
                        while(Sa <= r2*aa0)
                        {
                            reaction[0]++;
                            Sa = Sa + aa[reaction[0]];
                        }
                        timer++;
                        //cout << "Time to next event : " << tau << endl;
                        //cout << "Id of next event : " << reaction[0] << " with rate : " << aa[reaction[0]] << endl;
                    }
                    else if(timer != 0 && timer < tau){timer++;}
                    else if (timer >= tau)
                    {
                        if(reaction[0]%2 !=0)
                        {
                            ArrayHandle<Dscalar2> h_area(spv->returnAreaPeri(),access_location::host,access_mode::read);
                            dParams[0] = noise.getRealUniform(0,PI);
                            reaction[0] = int(reaction[0]/2);
                            output1 = reaction[0];
                            output2 = h_area.data[reaction[0]].x;
                            //cout << "Birth of cell : " << output1 << " at timestep : " << timestep+initSteps << endl;
                            spv->cellDivision(reaction, dParams);
                            Ncells = spv->getNumberOfDegreesOfFreedom();
                            divCount++;
                        }
                        else
                        {
                            ArrayHandle<Dscalar2> h_area(spv->returnAreaPeri(),access_location::host,access_mode::read);
                            output1 = int(reaction[0]/2);
                            output2 = h_area.data[output1].x;
                            //cout << "Death of cell : " << output1 <<  " at timestep : " << timestep+initSteps << endl;
                            spv->cellDeath(int(reaction[0]/2));
                            Ncells = spv->getNumberOfDegreesOfFreedom();
                        }
                        // Reset timer and reaction id to zero
                        timer =-relaxtime;//Allow some steps of relaxation before choosing the events
                        reaction[0] =0;
                        // Reset the scaling based on the average area
                        // now in the system. This is to keep the p0 as
                        // perimeter/sqrt(area), a geometrical factor ->
                        // We don't want that in this simulation
                        /*
                        vector<Dscalar2> APdyn(Ncells); 
                        Dscalar meanA = numpts / (Dscalar) Ncells;
                        ArrayHandle<int> h_ct(spv->cellType,access_location::host,access_mode::read);
                        for (int ii = 0; ii < Ncells; ++ii)
                        {
                            if(h_ct.data[ii]==0)
                            {
                                APdyn[ii].x=1;
                                APdyn[ii].y=p0_0 * sqrt(meanA);
                            }
                            else if(h_ct.data[ii]==1)
                            {
                                APdyn[ii].x=alpha;
                                APdyn[ii].y=p0_1 * sqrt(meanA);
                            }
                        }
                        //Setting cell preferences depending on the type
                        spv->setCellPreferences(APdyn);
                        */
                    }
                }
                else
                {
                    if(program_switch >0 && timestep%((int)(divisionTime/dt))==0) 
                    {
                        ArrayHandle<Dscalar2> h_area(spv->returnAreaPeri(),access_location::host,access_mode::read);
                        ArrayHandle<int> h_type(spv->cellType,access_location::host,access_mode::read);
                        // Divisions and death based on random process
                        divtime = noise.getRealUniform(0,1);
                        deathtime = noise.getRealUniform(0,1);
                        //if(program_switch >0 && timestep%((int)(divisionTime/dt))==0)
                        if(program_switch >0 && divtime < DivCut)
                        {
                            cdtest[0] = noise.getInt(0,Ncells-1);
                            dParams[0] = noise.getRealUniform(0,PI);
                            if(h_type.data[cdtest[0]] == 1)
                            {
                                // Lower division rate of liquid-like cells
                                if(divtime < DivCut-adv)
                                {
                                    spv->cellDivision(cdtest,dParams);
                                    Ncells = spv->getNumberOfDegreesOfFreedom();
                                }
                            } else {
                                spv->cellDivision(cdtest,dParams);
                                Ncells = spv->getNumberOfDegreesOfFreedom();
                            }
                        }
                        //Adjust death rate to the average of the two rates,
                        //so that there is no decrease of population.
                        if(program_switch >0 && deathtime < DeathCut)
                        {
                            int deadIdx;
                            if (pressure_dep == 0)
                            {
                                //ADDING DEATH TO THE MIX. Possibility to kill new cell.
                                deadIdx = noise.getInt(0,Ncells-1);
                                //printf("killing cell %i\n", deadIdx);
                            }
                            else
                            {
                                //Choosing death cell based on area of cells.
                                //First create probability distrib based on area of
                                //cell, here exponential.
                                vector<Dscalar> Proba(Ncells);
                                for (int jj = 0; jj<Ncells; ++jj)
                                {
                                    Proba[jj] = exp(-h_area.data[jj].x/(0.8*a0));
                                }
                                // Now draw the index from this distribution, TO BE
                                // FINISHED
                                discrete_distribution<int> distribution (Proba.begin(),Proba.end());
                                if (reproducible) {deadIdx = distribution(noise.gen);}
                                else {deadIdx = distribution(noise.genrd);}
                                //cout << "death of cell " << deadIdx << " that has area " << h_area.data[deadIdx].x <<endl;
                            }
                            spv->cellDeath(deadIdx);
                            Ncells = spv->getNumberOfDegreesOfFreedom();
                            //suppose, for instance, you want to keep p_0/sqrt(<A>) constant...
                            //Dscalar meanA = numpts / (Dscalar) Ncells;
                            //ArrayHandle<int> h_ct(spv->cellType,access_location::host,access_mode::read);
                            //for (int ii = 0; ii < numpts; ++ii)
                            //{
                            //    if(h_ct.data[ii]==0)
                            //    {
                            //        AP[ii].y=p0_0 * sqrt(meanA);
                            //    }
                            //    else if(h_ct.data[ii]==1)
                            //    {
                            //        AP[ii].y=p0_1 * sqrt(meanA);
                            //    }
                            //}
                            ////Setting cell preferences depending on the type
                            //spv->setCellPreferences(AP);
                            //Dscalar scaledP0 = p0 * sqrt(meanA);
                            //spv->setCellPreferencesUniform(1.0,scaledP0);
                            //printf("Ncells = %i\t <A> = %f \t p0 = %f\n",Ncells,meanA,scaledP0);
                            /*
                            //An alternative would be to use something like the following to rescale the box size to keep <A> = 1,
                            //and not rescale the preferred perimeter
                            BoxPtr newbox = make_shared<gpubox>(sqrt(Ncells),sqrt(Ncells));
                            sim->setBox(newbox);
                            */
                        };
                    }
                }
            }
            AverageNc+=Ncells;
            if(program_switch == 2 && timestep%((int)(saveTime/dt))==0)
            {
                //cout << "Saving at : " << timestep << endl;
                RunningAverage = AverageNc/(timestep+1);
                db1.WriteState(spv,p0_0,p0_1,KPP_0,KPP_1,alpha,KAA_0,KAA_1,Ph_1,lambda, divrate,Ph_adv,v0,DrM,init_fraction,RunningAverage, cellPressure);//Need to put back cellPressure
            };
        };
        AverageNc = AverageNc/tSteps;

        t2=clock();
        cout << "final number of cells= " <<spv->getNumberOfDegreesOfFreedom() << endl;
        cout << "With division = " << division << endl;
        cout << "timestep time per iteration currently at " <<  (t2-t1)/(Dscalar)CLOCKS_PER_SEC/tSteps << endl << endl;
        cout << "Number of generations " << divCount/numpts << endl;
        cout << "Ph =  " << Ph_1 << endl;
        cout << "AverageNc =  " << AverageNc << endl;

        };

    //program_switch < 0 --> self-propelled vertex model
    /*
    if(program_switch <0)
        {
        //...this section does the usual business of setting up the simulation
        int Nvert = 2*numpts;
        //here we'll demonstrate the more awkward task of using a sequence of netCDF databases to record what's going on
        AVMDatabaseNetCDF ncdat(Nvert,dataname,NcFile::Replace);
        bool runSPV = false;

        EOMPtr spp = make_shared<selfPropelledCellVertexDynamics>(numpts,Nvert);
        shared_ptr<VertexQuadraticEnergy> avm = make_shared<VertexQuadraticEnergy>(numpts,1.0,4.0,reproducible,runSPV);
        avm->setCellPreferencesUniform(1.0,p0);
        avm->setv0Dr(v0,1.0);
        avm->setT1Threshold(0.04);

        //Labeling each cell with a type for lineages.
        vector<int> types(numpts);
        for (int ii = 0; ii < numpts; ++ii)
        { 
            types[ii]=ii;
        }
        avm->setCellType(types);

        SimulationPtr sim = make_shared<Simulation>();
        sim->setConfiguration(avm);
        sim->addUpdater(spp,avm);
        sim->setIntegrationTimestep(dt);
        sim->setSortPeriod(initSteps/10);
        sim->setCPUOperation(!initializeGPU);
        sim->setReproducible(reproducible);
        //initial time steps
        for (int timestep = 0; timestep < initSteps+1; ++timestep)
            {
            sim->performTimestep();
            if(program_switch < -1 && timestep%((int)(1/dt))==0)
                {
                cout << timestep << endl;
                ncdat.WriteState(avm);
                };
            };

        //in the vertex model cell division, an integer vector (of length 3) is used.
        //the first indicates the cell to divide, and the next two indicate the vertices at the CCW
        //start of the edge to get new vertices. See vertexModelBase for details
        vector<int> cdtest(3); cdtest[0]=10; cdtest[1] = 0; cdtest[2] = 2;
        avm->cellDivision(cdtest);

        t1=clock();
        int Nvertices = avm->getNumberOfDegreesOfFreedom();
        int Ncells = Nvertices/2;
        int fileidx=2;
        int divisionTime = 10; //In physical units
        int relaxTime = divisionTime/(10*dt); //Relaxation time for division in physical units
        int saveTime = tSteps/10; // Time for saving files in computing units, here a generation
        //int saveTime = 10;
        for (int timestep = 0; timestep < tSteps; ++timestep)
            {
            sim->performTimestep();
            if(program_switch <=-1 && timestep%((int)(divisionTime/dt))==0)
                {
                cdtest[0] = noise.getInt(0,Ncells);
                avm->cellDivision(cdtest);

                Nvertices = avm->getNumberOfDegreesOfFreedom();
                Ncells = Nvertices/2;
                int deadCell = noise.getInt(0,Ncells);
                //cout << "targeting cell " << deadCell << endl;
                Dscalar2 oldAP; oldAP.x=a0; oldAP.y = p0;
                vector<Dscalar2> newPrefs(Ncells,oldAP);
                newPrefs[deadCell].x = 0.0;
                newPrefs[deadCell].y = p0*0.1;
                avm->setCellPreferences(newPrefs);
                int CellVertices = 0;
                //run for up to one tau to see if the cell shrinks to a triangle...if not, restore
                //the area and perimeter preference to stop the simulation from becoming unstable
                // INTRODUCED RELAX SET TO 10 TAU TO SEE IF THIS MAKES A DIFFERENCE
                // IN NUMBER OF CELLS KILLED
                for (int tt =0; tt < relaxTime; ++tt)
                    {
                        {
                        ArrayHandle<int> cn(avm->cellVertexNum,access_location::host,access_mode::read);
                        CellVertices = cn.data[deadCell];
                        };
                    if(CellVertices==3) break;
                    sim->performTimestep();
                    };
                if(CellVertices ==3)
                    {
                    //if(program_switch ==-2)
                    //    {
                    //    char dataname2[256];
                    //    sprintf(dataname2,"../results/test%i.nc",fileidx);
                    //    fileidx +=1;
                    //    AVMDatabaseNetCDF ncdat2(avm->getNumberOfDegreesOfFreedom(),dataname2,NcFile::Replace);
                    //    ncdat2.WriteState(avm);
                    //    };
                    //cout << "killing cell " << deadCell << endl;
                    avm->cellDeath(deadCell);
                    }
                else
                    {//if the cell doesn't die, restore it's area preferences
                    Dscalar2 oldAP; oldAP.x=a0; oldAP.y = p0;
                    vector<Dscalar2> newPrefs(Ncells,oldAP);
                    avm->setCellPreferences(newPrefs);
                    };
                };
            if(program_switch == -2 && timestep%((int)(saveTime))==0)
                {
                //cout << timestep << endl;
                char dataname2[256];
                sprintf(dataname2,"../results/test%i.nc",fileidx);
                fileidx +=1;
                AVMDatabaseNetCDF ncdat2(avm->getNumberOfDegreesOfFreedom(),dataname2,NcFile::Replace);
                ncdat2.WriteState(avm);
                };
            };

        t2=clock();
        cout << "final number of vertices = " <<avm->getNumberOfDegreesOfFreedom() << endl;
        cout << "timestep time per iteration currently at " <<  (t2-t1)/(Dscalar)CLOCKS_PER_SEC/tSteps << endl;
        avm->reportMeanVertexForce();
        cout << "Mean q = " << avm->reportq() << endl;
        };
    */

    if(initializeGPU)
        cudaDeviceReset();
    return 0;
    };
