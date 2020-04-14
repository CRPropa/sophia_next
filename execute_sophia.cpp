#include "sophia_interface.h"

int main() {

    std::cout.precision(10);
    int nEvent = 10;

    for (int k = 0; k < nEvent; ++k) {

        bool onProton = true;
        // bool onProton = false;
        double Ein = 1e9;  // GeV
        double eps = 1e-9;  // GeV
        bool declareChargedPionsStable = true;
        // bool declareChargedPionsStable = false;

        sophia_interface SI;

        std::cout << "Event # " << k + 1 << std::endl;

        sophiaevent_output seo = SI.sophiaevent(onProton, Ein, eps, declareChargedPionsStable);

        std::cout << "\nevent result # " << (k + 1) << std::endl;
        int Nout = seo.Nout;
        std::cout << "N = " <<Nout << std::endl;
        for (int i = 0; i < Nout; ++i) {
            std::cout << "ID = " << ID_sophia_to_PDG(seo.outPartID[i]) << std::endl;
            for (int j = 0; j < 6; ++j) {
                if (j == 3) std::cout << seo.outPartP[j][i] << std::endl;
            }
        }
        std::cout << "\n---------------------------------" << std::endl;
    }
}
