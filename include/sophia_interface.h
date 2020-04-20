#ifndef SOPHIA_INTERFACE_H
#define SOPHIA_INTERFACE_H

#include <iostream>
#include <cmath>
#include <string>
#include <algorithm>
#include <initializer_list>
#include <vector>
#include <functional>

#include "sophia_data.h"

struct sophiaevent_output {
    double outPartP[5][2000];
    int outPartID[2000];
    int Nout;
    double getPartP(int i, int j) const {
	return outPartP[i][j];
    }
    int getPartID(int i) const {
	return outPartID[i];
    }
};

struct DECPAR_zero_output {
    double P_out[5][10];
};

struct DECPAR_nonZero_output {
    int ND;
    int LL[10];
    double P_out[5][10];
};

struct PO_MSHELL_output {
    std::vector<double> P1;
    std::vector<double> P2;
};

struct dec_proc2_output {
    int IPROC;
    int IRANGE;
};

struct PO_ALTRA_output {
    double PX;
    double PY;
    double PZ;
    double P;
    double E;
};

struct PO_SELSX2_output {
    double XS1[2];
    double XS2[2];
    bool isRejected;
};

struct lund_get_output {
        double PX;
        double PY;
        double PZ;
        double EE;
        double XM;
        int IFL;
};

// - input:  function over which to integrate, integration limits A and B
// - output: 8-points Gauß-Legendre integral
static const double X[8] = {.0950125098, .2816035507, .4580167776, .6178762444, .7554044083, .8656312023, .9445750230, .9894009349};
static const double W[8] = {.1894506104, .1826034150, .1691565193, .1495959888, .1246289712, .0951585116, .0622535239, .0271524594};
template<typename Integrand>
double gaussInt(Integrand&& integrand, double A, double B) {
    const double XM = 0.5 * (B + A);
    const double XR = 0.5 * (B - A);
    double SS = 0.;
    for (int i = 0; i < 8; ++i) {
        double DX = XR * X[i];
        SS += W[i] * (integrand(XM + DX) + integrand(XM - DX));
    }
    return XR * SS;
}

int ID_sophia_to_PDG(int sophiaID);

class sophia_interface {
public:

    void debug(std::string, bool stopProgram = false);
    void debugNonLUND(std::string, bool stopProgram = false);

    // JETSET
    int N;
    int K[5][4000];
    double P[5][4000];
    double V[5][4000];

    void LU2ENT(int KF1, int KF2, double PECM);  // prototype qq interaction method;
    void LUEXEC();
    void LUDECY(int IP);
    void LUDECY_setupPartonShowerEvolution(int IP, int NSAV, int MMAT, int ND, int MMIX);
    void LUSTRF(int IP);
    void LUINDF(int IP);
    void LUPREP(int IP);
    void LUPREP_checkFlavour(int IP);
    double LUZDIS(int KFL1, int KFL2, double PR);
    std::vector<double> LUPTDI(int KFL);
    std::vector<int> LUKFDI(int KFL1, int KFL2);
    void LUEDIT(int MEDIT);
    void LUSHOW(int IP1, int IP2, double QMAX);
    void LUBOEI(int NSAV);
    void LUJOIN(int NJOIN, int IJOIN[]);
    void LUERRM(int MERR, std::string CHMESS);
    void LUROBO(double THE, double PHI, double BEX, double BEY, double BEZ);
    void LUDBRB(int IMIN, int IMAX, double THE, double PHI, double DBX, double DBY, double DBZ, bool skip=false);
    int KLU(int I, int J);
    double RLU(bool isCalledByRNDM = false);  // defalut, internal RNG. To be removed later.
    double ULMASS(int KF);
    double ULANGL(double X, double Y);
    double PLU(int I, int J);
    int LUCHGE(int KF);
    int LUCOMP(int KF);

    // SOPHIA
    int np;
    double p[5][2000];
    int LLIST[2000];

    sophiaevent_output sophiaevent(bool onProton, double Ein, double eps, bool declareChargedPionsStable=false);
    void eventgen(int L0, double E0, double eps, double theta);
    void gamma_h(double Ecm, int ip1, int Imode);
    void DECSIB();
    DECPAR_zero_output DECPAR_zero(double P0[5], int ND, int LL[10]);
    DECPAR_nonZero_output DECPAR_nonZero(int LA, double P0[5]);
    PO_MSHELL_output PO_MSHELL(double PA1[4], double PA2[4], double XM1, double XM2);
    std::vector<int> valences(int ip);
    void check_event(int Ic, double Esum, double PXsum, double PYsum, double PZsum, int IQchr, int IQbar);
    int dec_res2(double eps_prime, int IRESMAX, int L0);
    void RES_DECAY3(int IRES, int IPROC, int IRANGE, double s, int L0);
    void PROC_TWOPART(int LA, int LB, double AMD, double cosAngle);
    int dec_inter3(double eps_prime, int L0);
    double sample_s(double eps, int L0, double Ein);
    double functs(double s, int L0);
    double crossection(double x, int NDIR, int NL0);
    dec_proc2_output dec_proc2(double x, int IRES, int L0);
    double singleback(double x);
    double twoback(double x);
    double scatterangle(int IRES, int L0);
    double probangle(int IRES, int L0, double z);
    PO_ALTRA_output PO_ALTRA(double GA, double BGX, double BGY, double BGZ, double PCX, double PCY, double PCZ, double EC);
    std::vector<double> PO_TRANS(double XO, double YO, double ZO, double CDE, double SDE, double CFE, double SFE);
    PO_SELSX2_output PO_SELSX2(double XMIN[2], double XMAX[2], double AS1, double AS2);
    double PO_RNDBET(double GAM, double ETA);
    double PO_RNDGAM(double ETA);
    void lund_frag(double SQS);
    void lund_put(int I, int IFL, double PX, double PY, double PZ, double EE);
    lund_get_output lund_get(int I);
    int ICON_PDG_SIB(int ID);
    double PO_XLAM(double X, double Y, double Z);
    double RNDM();
    double Pl(double x, double xth, double xmax, double alpha);
    double Ef(double x, double th, double w);
    double breitwigner(double sigma_0, double Gamma, double DMM, double eps_prime);
};

#endif
