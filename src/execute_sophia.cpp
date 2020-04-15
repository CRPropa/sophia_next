#include <iostream>
#include <fstream>

#include "sophia_interface.h"

int main() {
    std::cout.precision(10);

    bool onProton = true;
    // bool onProton = false;
    double Ein = 1e9;  // GeV
    double eps = 1e-9;  // GeV
    bool declareChargedPionsStable = true;
    // bool declareChargedPionsStable = false;

    std::ofstream outfile;
    outfile.open("outData.csv");

    outfile << "eventID\t" << "partID\t" << "EGeV\n";

    int nEvent = 10000;
    for (int k = 0; k < nEvent; ++k) {

        sophia_interface SI;

        // std::cout << "Event # " << k + 1 << std::endl;

        sophiaevent_output seo = SI.sophiaevent(onProton, Ein, eps, declareChargedPionsStable);


        // std::cout << "\nevent result # " << (k + 1) << std::endl;
        int Nout = seo.Nout;
        // std::cout << "N = " <<Nout << std::endl;
        for (int i = 0; i < Nout; ++i) {
        //     std::cout << "ID = " << ID_sophia_to_PDG(seo.outPartID[i]) << std::endl;
            for (int j = 0; j < 5; ++j) {
                // if (j == 3) std::cout << seo.outPartP[j][i] << std::endl;
                if (j == 3) {
                    int eventID = k + 1;
                    outfile << eventID << "\t"
                            << ID_sophia_to_PDG(seo.outPartID[i]) << "\t"
                            << seo.outPartP[j][i] << "\n";
                }
            }
        }
        // std::cout << "\n---------------------------------" << std::endl;
    }
    outfile.close();
}
