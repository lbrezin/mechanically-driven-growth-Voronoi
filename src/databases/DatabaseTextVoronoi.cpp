#define ENABLE_CUDA
#define _GLIBCXX_USE_CXX11_ABI 0
#include "DatabaseTextVoronoi.h"
/*! \file DatabaseTextVoronoi.cpp */

DatabaseTextVoronoi::DatabaseTextVoronoi(string fn, int mode)
    : BaseDatabase(fn,mode)
    {
    switch(Mode)
        {
        case -1:
            inputFile.open(filename);
            break;
        case 0:
            outputFile.open(filename);
            break;
        case 1:
            outputFile.open(filename, ios::app);
            break;
        default:
            ;
        };
    };

void DatabaseTextVoronoi::WriteState(STATE s, Dscalar p0_0, Dscalar p0_1, Dscalar gamma_1, Dscalar gamma_2, Dscalar alpha, Dscalar K1, Dscalar K2, Dscalar homeoP, Dscalar lambda, Dscalar divrate, Dscalar Ph_adv, Dscalar v0, Dscalar DrM, Dscalar init_fraction, Dscalar RunningAverage, vector<Dscalar> cellPressure, Dscalar time, int rec)
    {
    if (rec != -1)
        {
        printf("writing to the middle of text files not supported\n");
        throw std::exception();
        };
        
    int N = s->getNumberOfDegreesOfFreedom();
//    printf("saving %i cells\n",N);
    if (time < 0) time = s->currentTime;
    Dscalar x11,x12,x21,x22;
    s->returnBox().getBoxDims(x11,x12,x21,x22);

    //Add first two lines with mechanical parameters at the first time
    //it is called, just after performing a time step
    if (time==s->deltaT){
        outputFile << "p0_0 " << "p0_1 " << "gamma_1 " << "gamma_2 "<< "alpha " << "K1 " << "K2 " << "Ph  "<< "lambda  " << "divrate "<< "Ph_adv "<< "v0 "<< "DrM "<<"f0" << endl;
        outputFile << p0_0 << " " << p0_1 << " " <<gamma_1 << " " <<gamma_2 << " " <<alpha << " " <<K1 << " " <<K2 << " " << homeoP << " " << 
            lambda << " " << divrate << " " << Ph_adv << " " << v0 << " " << DrM << " " << init_fraction << endl;
    }
    //Changed tab to comma to make a csv file
    outputFile << "N, time, box, Average\n";
    outputFile << N << " " << time <<" " <<x11 <<" " <<x12<<" " <<x21<<" " <<x22 <<" " <<RunningAverage <<"\n";

    ArrayHandle<Dscalar2> h_pos(s->cellPositions,access_location::host,access_mode::read);
    ArrayHandle<int> h_ct(s->cellType,access_location::host,access_mode::read);
    ArrayHandle<int> h_nei(s->cellNeighborNum,access_location::host,access_mode::read);
    //ArrayHandle<int> h_cl(s->cellLineage,access_location::host,access_mode::read);
    //Adding area and perimeter to the output file.
    ArrayHandle<Dscalar2> h_area(s->returnAreaPeri(),access_location::host,access_mode::read);
    //ArrayHandle<int> h_itt(s->returnItt(),access_location::host,access_mode::read);
    for (int ii = 0; ii < N; ++ii)
        {
        int pidx = s->tagToIdx[ii];
        //outputFile <<h_pos.data[pidx].x <<" "<<h_pos.data[pidx].y <<" " <<h_ct.data[pidx] << " " << h_area.data[pidx].x << " " <<h_area.data[pidx].y << " " <<h_nei.data[pidx] <<" " << pidx << "\n";
        //outputFile <<h_pos.data[ii].x <<" "<<h_pos.data[ii].y <<" " <<h_ct.data[ii] << " " << h_area.data[ii].x << " " <<h_area.data[ii].y << " " <<h_nei.data[ii] <<" " << s->returnItt()[ii] << "\n";
        outputFile <<h_pos.data[pidx].x <<" "<<h_pos.data[pidx].y <<" " <<h_ct.data[pidx] << " " << h_area.data[pidx].x << " " <<h_area.data[pidx].y << " " <<h_nei.data[pidx] <<" " << s->returnItt()[ii] <<  " " << cellPressure[pidx] << "\n";
        };
    };

void DatabaseTextVoronoi::ReadState(STATE s, int rec, bool geometry)
    {
    printf("Reading from a text database currently not supported\n");
    throw std::exception();
    };
