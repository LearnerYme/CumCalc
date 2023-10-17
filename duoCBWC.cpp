/*
    Updated: Now read file lists (of input raw cumulants and centrality edges)
    Support up to 10 files.

    CBWC for U/L runs.
    Will do weighted average including both U/L runs and merge them together with specified bin edge.
    Yige Huang on Aug. 18, 2023

    A new version of cumulant calculating.
    Yige Huang on Dec. 22, 2022
    Based on Yu Zhang's code.
*/

#include <iostream>
#include <fstream>
#include <string>
#include <vector>

#include "TFile.h"
#include "TH1D.h"
#include "TTree.h"
#include "TGraphErrors.h"

#include "NpartLoader.h"
#include "CentDefinition.h"

using std::vector;
using std::string;

int main(int argc, char** argv){
    if (argc != 4) {
        std::cout << "[ERROR] Should have 3 arguments!\n";
        std::cout << " - Usage: ./duoCBWC FILE_LIST CENT_LIST OUTNAME\n";
    } else {
        std::cout << "[LOG] Input inforamtion: \n";
        std::cout << " - Terms file list: " << argv[1] << std::endl;
        std::cout << " - Cent. edge file list: " << argv[2] << std::endl;
        std::cout << "[LOG] Output file name: " << argv[3] << ".root" << std::endl;
    }

    vector<string> filelist;
    vector<string> centlist;
    TFile* tfs[10];
    CentDefinition* cds[10];
    string tmpStr;

    int cnt = 0;
    std::ifstream rfilelist;
    rfilelist.open(argv[1]);
    while(std::getline(rfilelist, tmpStr)) {
        filelist.push_back(tmpStr);
        tfs[cnt] = new TFile(tmpStr.c_str());
        cnt ++;
    }

    cnt = 0;
    std::ifstream rcentlist;
    rcentlist.open(argv[2]);
    while(std::getline(rcentlist, tmpStr)) {
        centlist.push_back(tmpStr);
        cds[cnt] = new CentDefinition();
        cds[cnt]->read_edge(tmpStr.c_str());
        cnt ++;
    }


    // CBWC
    const int MaxMult = 2000;
    const int LowEventCut = 10; // to avoid error caused by low event number
    std::cout << "[LOG] Now applying CBWC.\n";
    const int nCent = 9;
    NpartLoader* nDef = new NpartLoader();
    int* nPart = nDef->GetArray();

    std::cout << "[LOG] Initializing graphes.\n";
    const int nCums = 22; // cums for cumulant and correlation function (orders)
    const int nType = 3; // type for proton type
    const char* typeNames[nType] = {"Pro", "Pbar", "Netp"};
    TGraphErrors* tgs[3][nCums];
    TH1D* sCums[3][nCums][10];
    const char* cumNames[nCums] = {
        "C1", "C2", "C3", "C4", "C5", "C6",
        "R21", "R32", "R42", "R51", "R62",
        "k1", "k2", "k3", "k4", "k5", "k6",
        "k21", "k31", "k41", "k51", "k61"
    };
    const char* cumTitles[nCums] = {
        "C_{1}", "C_{2}", "C_{3}", "C_{4}", "C_{5}", "C_{6}",
        "C_{2}/C_{1}", "C_{3}/C_{2}", "C_{4}/C_{2}", "C_{5}/C_{1}", "C_{6}/C_{2}",
        "#kappa_{1}", "#kappa_{2}", "#kappa_{3}", "#kappa_{4}", "#kappa_{5}", "#kappa_{6}",
        "#kappa_{2}/#kappa_{1}", "#kappa_{3}/#kappa_{1}", "#kappa_{4}/#kappa_{1}", "#kappa_{5}/#kappa_{1}", "#kappa_{6}/#kappa_{1}"
    };

    for (int i=0; i<nType; i++){
        for (int j=0; j<nCums; j++){
            if (i == 2 && j >= 11){
                continue; // skip kappa for net proton
            }
            for (int k=0; k<cnt; k++) {
                tfs[k]->GetObject(
                    Form("%s%s", typeNames[i], cumNames[j]), 
                    sCums[i][j][k]
                );
            }
            tgs[i][j] = new TGraphErrors(nCent);
            tgs[i][j]->SetName(
                Form("%s_%s", typeNames[i], cumNames[j])
            );
            tgs[i][j]->SetTitle(cumTitles[j]);
            tgs[i][j]->SetMarkerStyle(20);
        }
    }

    const int LowMultCut = 1;
    double CentEvent[nType][nCent] = {0};
    double vSum[nType][nCent]; // value 
    double eSum[nType][nCent]; // error 

    // get number of events
    TH1D* hEntries[nType][10];
    for (int i=0; i<nType; i++){
        for (int j=0; j<cnt; j++) {
            tfs[j]->GetObject(Form("%shEntries", typeNames[i]), hEntries[i][j]);
        }
    }

    for (int i=0; i<nType; i++){ // i for proton types
        for (int j=0; j<nCums; j++){ // j for cumulant orders
            if (i == 2 && j >= 11){ // as well, skip kappa for net proton
                continue;
            }

            // loop 1, initialize the arrays to carry values and errors
            int nEvents = 0;
            for (int k=0; k<nCent; k++){ // k for centralities
                vSum[i][k] = 0;
                eSum[i][k] = 0;
                CentEvent[i][k] = 0;
            }
            
            // loop 2, sum up in different centrality bins
            for (int k=LowMultCut; k<=MaxMult; k++){ // here k for refmult3
                // inner loop: files
                for (int p=0; p<cnt; p++) {
                    nEvents = hEntries[i][p]->GetBinContent(k+1); // need plus 1 here, 0 is the first bin
                    int curCent = cds[p]->get_cent(k); // get current centrality
                    if (curCent >= 0 && nEvents >= LowEventCut){
                        if (
                            (j < 6) || // C1 ~ C6
                            (j >= 11 && j <= 16) // k1 ~ k6
                        ){ // for ratios, the value is just a dividing results, no need to do cbwc
                            vSum[i][curCent] += (sCums[i][j][p]->GetBinContent(k+1) * nEvents);
                        }
                        // but the errors should be calculated
                        eSum[i][curCent] += pow((sCums[i][j][p]->GetBinError(k+1) * nEvents), 2);
                        CentEvent[i][curCent] += nEvents;
                    }
                }
            }

            // loop 3, get the averaged results (by dividing number of events)
            for (int k=0; k<nCent; k++){
                // calculate values
                if (
                    (j < 6) || 
                    (j >= 11 && j <= 16)
                ){
                    vSum[i][k] /= CentEvent[i][k];
                } else if (j == 6) {
                    // vSum[i][k] = tgs[i][1]->GetPointY(k) * 1.0 / tgs[i][0]->GetPointY(k); // C2 / C1
                    vSum[i][k] = *(tgs[i][1]->GetY()+k) / *(tgs[i][0]->GetY()+k); // C2 / C1
                } else if (j == 7) {
                    vSum[i][k] = *(tgs[i][2]->GetY()+k) / *(tgs[i][1]->GetY()+k); // C3 / C2
                } else if (j == 8) {
                    vSum[i][k] = *(tgs[i][3]->GetY()+k) / *(tgs[i][1]->GetY()+k); // C4 / C2
                } else if (j == 9) {
                    vSum[i][k] = *(tgs[i][4]->GetY()+k) / *(tgs[i][0]->GetY()+k); // C5 / C1
                } else if (j == 10) {
                    vSum[i][k] = *(tgs[i][5]->GetY()+k) / *(tgs[i][1]->GetY()+k); // C6 / C2
                } else if (j >= 17) {
                    vSum[i][k] = *(tgs[i][j-5]->GetY()+k) / *(tgs[i][11]->GetY()+k); // k(2~6) / k1
                }
                // calculate errors
                eSum[i][k] = sqrt(eSum[i][k]) / CentEvent[i][k];
                // set points to graphes
                tgs[i][j]->SetPoint(k, nPart[k], vSum[i][k]);
                tgs[i][j]->SetPointError(k, 0.0, eSum[i][k]);
            }
        }
    }

    std::cout << "[LOG] Calculating finished, now saving.\n";
    TFile* tfout = new TFile(Form("%s.root", argv[5]), "recreate");
    tfout->cd();
    for (int i=0; i<nType; i++){
        for (int j=0; j<nCums; j++){
            if (i == 2 && j >= 11){
                continue;
            }
            tgs[i][j]->Write();
        }
    }
    for (int i=0; i<cnt; i++) {
        tfs[i]->Close();
    }
    tfout->Close();

    std::cout << "[LOG] All done!.\n";

    return 0;
}

