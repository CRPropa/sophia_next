#ifndef SOPHIA_INTERFACE_H
#define SOPHIA_INTERFACE_H

/**
    This is a C++ version of SOPHIA.
    Translated by Mario Hoerbe (mario.hoerbe@ruhr-uni-bochum.de)
    2020
*/

/**
c*****************************************************************************
c**!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!***
c**!!              IF YOU USE THIS PROGRAM, PLEASE CITE:                 !!***
c**!! A.M"ucke, Ralph Engel, J.P.Rachen, R.J.Protheroe and Todor Stanev, !!***
c**!!  1999, astro-ph/9903478, to appear in Comp.Phys.Commun.            !!***
c**!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!***
c*****************************************************************************
c** Further SOPHIA related papers:                                         ***
c** (1) M"ucke A., et al 1999, astro-ph/9808279, to appear in PASA.        ***
c** (2) M"ucke A., et al 1999, to appear in: Proc. of the                  ***
c**      19th Texas Symposium on Relativistic Astrophysics, Paris, France, ***
c**      Dec. 1998. Eds.: J.~Paul, T.~Montmerle \& E.~Aubourg (CEA Saclay) ***
c** (3) M"ucke A., et al 1999, astro-ph/9905153, to appear in: Proc. of    ***
c**      19th Texas Symposium on Relativistic Astrophysics, Paris, France, ***
c**      Dec. 1998. Eds.: J.~Paul, T.~Montmerle \& E.~Aubourg (CEA Saclay) ***
c** (4) M"ucke A., et al 1999, to appear in: Proc. of 26th Int.Cosmic Ray  ***
c**      Conf. (Salt Lake City, Utah)                                      ***
c*****************************************************************************
*/


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

    // debug functions outputting the central variables (i.e. particle status).
    // can be opted to stop the program at where they were called for debugging
    // these functions must not be used if no debugging is going on
    void debug(std::string message, bool stopProgram = false);
    void debugNonLUND(std::string message, bool stopProgram = false);

    /*
        The seubsequent functions belong to JETSET v7.4 which is imbedded into what later was called SOPHIA.
        The central variables, which for the time being are implemented as globals,
        are listed directly below - none of these are used in the "SOPHIA" function cluster
    */
    int N;              // number of current particles in the event
    int K[5][4000];     // their current color for QCD purposes
    double P[5][4000];  // their current 5-momenta (px, py, pz, E, m0) in GeV
    double V[5][4000];  // their vertices for QCD purposes

    // prototype qq interaction method (not used in p/gamma event generation);
    // This is a separate entry to JETSET
    // Purpose: to store two partons/particles in their CM frame, with the first along the +z axis. 
    void LU2ENT(int KF1, int KF2, double PECM);

    // Purpose: to administrate the fragmentation and decay chain.
    void LUEXEC();

    // Purpose: to handle the decay of unstable particles.
    void LUDECY(int IP);

    // Set up for parton shower evolution from jets.
    // This function is only called by LUDECY. After this function was called, LUDECY returns.
    void LUDECY_setupPartonShowerEvolution(int IP, int NSAV, int MMAT, int ND, int MMIX);
    
    // Purpose: to handle the fragmentation of an arbitrary colour singlet 
    // jet system according to the Lund string fragmentation model. 
    // NOTE: this is kond-of the main method of JETSET. It is the only
    // function in which there are still goto pointers left
    void LUSTRF(int IP);

    // Purpose: to handle the fragmentation of a jet system (or a single jet)
    // according to independent fragmentation models.
    void LUINDF(int IP);

    // Purpose: to rearrange partons along strings, to allow small systems 
    // to collapse into one or two particles and to check flavours. 
    void LUPREP(int IP);

    // Serves as a local function of LUPREP. Solves goto->320 statement.
    // When this function is called, LUPREP returns.
    void LUPREP_checkFlavour(int IP);

    // Purpose: to generate the longitudinal splitting variable Z.
    double LUZDIS(int KFL1, int KFL2, double PR);

    // Purpose: to generate transverse momentum according to a Gaussian.
    // output: vector = (PX, PY)
    std::vector<double> LUPTDI(int KFL);


    // Purpose: to generate a new flavour pair and combine off a hadron.
    // Default flavour values. Input consistency checks.
    // Returns vector (KFL3, KF).
    std::vector<int> LUKFDI(int KFL1, int KFL2);

    // Purpose: to perform global manipulations on the event record,
    // in particular to exclude unstable or undetectable partons/particles.
    void LUEDIT(int MEDIT);

    // Purpose: to generate timelike parton showers from given partons. 
    void LUSHOW(int IP1, int IP2, double QMAX);

    // Purpose: to modify event so as to approximately take into account
    // Bose-Einstein effects according to a simple phenomenological
    // parametrization.
    void LUBOEI(int NSAV);

    // Purpose: to connect a sequence of partons with colour flow indices,
    // as required for subsequent shower evolution (or other operations).
    void LUJOIN(int NJOIN, int IJOIN[]);

    // Purpose: to inform user of errors in program execution.
    // VERY IMPORTANT NOTE: in the original version, all error messages have been set
    // to silently fail!!! This option was set with MSTU[20] = 0. While operating SOPHIA,
    // silent fails actually happend now and then!    
    void LUERRM(int MERR, std::string CHMESS);

    // Purpose: to perform rotations and boosts.
    void LUROBO(double THE, double PHI, double BEX, double BEY, double BEZ);

    // Entry of LUROBO for specific range and double precision boost. 
    void LUDBRB(int IMIN, int IMAX, double THE, double PHI, double DBX, double DBY, double DBZ, bool skip=false);

    // Purpose: to provide various integer-valued event related data.
    int KLU(int I, int J);

    // Purpose: to generate random numbers uniformly distributed between
    // 0 and 1, excluding the endpoints.
    // defalut, internal RNG. Useful for debugging and to be replaced later.
    double RLU(bool isCalledByRNDM = false);

    // Purpose: to give the mass of a particle/parton.
    double ULMASS(int KF);

    // Purpose: to reconstruct an angle from given x and y coordinates.
    double ULANGL(double X, double Y);

    // Purpose: to provide various real-valued event related data.
    double PLU(int I, int J);
    
    // Purpose: to give three times the charge for a particle/parton.
    int LUCHGE(int KF);

    // Purpose: to compress the standard KF codes for use in mass and decay 
    // arrays; also to check whether a given code actually is defined.   
    int LUCOMP(int KF);

    /*
        These are functions which are responsible for doing the "low-energy part" of SOPHIA or interface to LUND JETSET.
        Some methods are written by the authors of the original SOPHIA publications
        and others are directly taken from SIBYLL or PHOJET
    */
    int np;             // number of current particles
    double p[5][2000];  // 5-momenta of current particles (px, py, pz, E, m0) in GeV
    int LLIST[2000];    // particle IDs of current particles (SIBYLL convention)

    // ****************************************************************************
    //    SOPHIAEVENT
    // 
    //    interface between Sophia and CRPropa
    //    simulate an interaction between p/n of given energy and a photon
    // 
    //    Eric Armengaud, 2005
    //    modified & translated from FORTRAN to C++: Mario Hoerbe, 2020
    // *******************************
    //  onProton = primary particle is proton or neutron
    //  Ein = input energy of primary nucleon in GeV (SOPHIA standard energy unit)
    //  eps = input energy of target photon in GeV (SOPHIA standard energy unit)
    //  declareChargedPionsStable = pi+-0 are set to be stable particles. See array IDB for details
    //  OutPartP = list of 4-momenta + rest masses of output particles (neutrinos approx. 0)
    //  OutPartID = ID list of output particles (PDG IDs)
    //  Nout = number of output particles
    // ****************************************************************************
    // this is the only function a user should run
    sophiaevent_output sophiaevent(bool onProton, double Ein, double eps, bool declareChargedPionsStable=false);

    // *******************************************************
    // ** subroutine for photopion production of            **
    // ** relativistic nucleons in a soft photon field      **
    // ** subroutine for SOPHIA inVersion 1.2               **
    // ****** INPUT ******************************************
    //  E0 = energy of incident proton (in lab frame) [in GeV]
    //  eps = energy of incident photon [in GeV] (in lab frame)
    //  theta = angle between incident proton and photon [in degrees]
    //  L0 = code number of the incident nucleon
    // ****** OUTPUT *************************************************
    //  P(2000,5) = 5-momentum of produced particles 
    //  LLIST(2000) = code numbers of produced particles
    //  NP = number of produced particles
    // ***************************************************************
    // ** Date: 20/01/98       **
    // ** correct.:19/02/98    **
    // ** change:  23/05/98    **
    // ** last change:06/09/98 **
    // ** authors: A.Muecke    **
    // **          R.Engel     **
    // **************************
    void eventgen(int L0, double E0, double eps, double theta);
    
    // **********************************************************************
    // 
    //      simple simulation of low-energy interactions (R.E. 03/98)
    // 
    //      changed to simulate superposition of reggeon and pomeron exchange 
    //      interface to Lund / JETSET 7.4 fragmentation
    //                                                   (R.E. 08/98)
    // 
    //      input: ip1    incoming particle
    //                    13 - p
    //                    14 - n
    //             Ecm    CM energy in GeV
    //             Imode  interaction mode
    //                    0 - multi-pion fragmentation
    //                    5 - fragmentation in resonance region
    //                    1 - quasi-elastic / diffractive interaction 
    //                        (p/n-gamma  --> n/p rho)
    //                    4 - quasi-elastic / diffractive interaction 
    //                        (p/n-gamma  --> n/p omega)
    //                    2 - direct interaction (p/n-gamma  --> n/p pi)
    //                    3 - direct interaction (p/n-gamma  --> delta pi)
    //             IFBAD control flag
    //                   (0  all OK,
    //                    1  generation of interaction not possible)
    // 
    // **********************************************************************
    void gamma_h(double Ecm, int ip1, int Imode);
    
    // ***********************************************************************
    // Decay all unstable particle in SIBYLL
    // decayed particle have the code increased by 10000
    // (taken from SIBYLL 1.7, R.E. 04/98)
    // ***********************************************************************
    void DECSIB();

    // the routine generates a phase space decay of a particle
    // with 5-momentum P0(1:5) into ND particles of codes LL(1:nd)
    // (taken from SIBYLL 1.7, muon decay corrected, R.E. 04/98)    
    DECPAR_zero_output DECPAR_zero(double P0[5], int ND, int LL[10]);

    // ***********************************************************************
    // This subroutine generates the decay of a particle
    // with ID = LA, and 5-momentum P0(1:5)
    // into ND particles of 5-momenta P(j,1:5) (j=1:ND)
    // (taken from SIBYLL 1.7, muon decay corrected, R.E. 04/98)
    // ***********************************************************************
    DECPAR_nonZero_output DECPAR_nonZero(int LA, double P0[5]);

    // ********************************************************************
    // rescaling of momenta of two partons to put both on mass shell
    // input:       PA1,PA2   input momentum vectors
    //              XM1,2     desired masses of particles afterwards
    //              P1,P2     changed momentum vectors
    // (taken from PHOJET 1.12, R.E. 08/98)
    // ********************************************************************
    PO_MSHELL_output PO_MSHELL(double PA1[4], double PA2[4], double XM1, double XM2);
    
    // valence quark composition of various particles  (R.E. 03/98)
    // (with special treatment of photons)
    std::vector<int> valences(int ip);

    // check energy-momentum and quantum number conservation (R.E. 08/98) 
    void check_event(int Ic, double Esum, double PXsum, double PYsum, double PZsum, int IQchr, int IQbar);

    // *****************************************************************************
    // decides which resonance with ID=IRES in list takes place at eps_prime
    // ** Date: 20/01/98   **
    // ** author: A.Muecke **
    // **********************
    int dec_res2(double eps_prime, int IRESMAX, int L0);

    // ********************************************************
    // RESONANCE AMD with code number IRES  INTO  M1 + M2
    // PROTON ENERGY E0 [in GeV] IN DMM [in GeV]
    // E1,E2 [in GeV] are energies of decay products
    // LA,LB are code numbers of decay products
    // P(1,1:5),P(2,1:5) are 5-momenta of particles LA,LB;
    // resulting momenta are calculated in CM frame;
    // cosAngle is cos of scattering angle in CM frame
    // ********************************************************
    // ** Date: 20/01/98    **
    // ** correct.:28/04/98 **
    // ** author: A.Muecke  **
    // **********************    
    void RES_DECAY3(int IRES, int IPROC, int IRANGE, double s, int L0);

    // ***********************************************************
    // 2-particle decay of CMF mass AMD INTO  M1 + M2
    // nucleon energy E0 [in GeV];
    // E1,E2 [in GeV] are energies of decay products
    // LA,LB are code numbers of decay products
    // P1(1:5),P2(1:5) are 5-momenta of particles LA,LB;
    // resulting momenta are calculated in CM frame;
    // costheta is cos of scattering angle in CM frame
    // this program also checks if the resulting particles are
    // resonances; if yes, it is also allowed to decay a
    // mass AMD < M1 + M2 by using the width of the resonance(s)
    // ***********************************************************
    // ** Date: 20/01/98    **
    // ** correct.:19/02/98 **
    // ** author: A.Muecke  **
    // **********************
    void PROC_TWOPART(int LA, int LB, double AMD, double cosAngle);
    
    // *** decides which process takes place at eps_prime *********
    // (0) multipion production (fragmentation)                 ***
    // (1) diffractive scattering: N\gamma --> N \rho           ***
    // (2) direct pion production: N\gamma --> N \pi            ***
    // (3) direct pion production: N\gamma --> \Delta \pi       ***
    // (4) diffractive scattering: N\gamma --> N \omega         ***
    // (5) fragmentation in resonance region                    ***
    // (6) excitation/decay of resonance                        ***
    // ************************************************************
    // ** Date: 15/04/98   **
    // ** author: A.Muecke **
    // **********************
    int dec_inter3(double eps_prime, int L0);

    // ***********************************************************************
    //  samples distribution of s: p(s) = (s-mp^2)sigma_Ngamma
    //  rejection for s=[sth,s0], analyt.inversion for s=[s0,smax]
    // ***********************************************************************
    // ** Date: 20/01/98   **
    // ** author: A.Muecke **
    // **********************
    double sample_s(double eps, int L0, double Ein);

    // a function returning a scalaer related to the CM energy s.
    double functs(double s, int L0);

    // calculates crossection of Nucleon-gamma-interaction
    // (see thesis of J.Rachen, p.45ff and corrections
    // report from 27/04/98, 5/05/98, 22/05/98 of J.Rachen)
    // ** Date: 20/01/98   **
    // ** correct.:27/04/98**
    // ** update: 23/05/98 **
    // ** author: A.Muecke **    
    double crossection(double x, int NDIR, int NL0);

    // decide which decay with ID=IPROC of resonance IRES takes place ***
    // Date: 20/01/98   **
    // correct.: 27/04/98*
    // author: A.Muecke **
    dec_proc2_output dec_proc2(double x, int IRES, int L0);

    // single pion production
    double singleback(double x);

    // two pion production
    double twoback(double x);

    // *******************************************************************
    // This routine samples the cos of the scattering angle for a given **
    // resonance IRES and incident nucleon L0; it is exact for          **
    // one-pion decay channel and if there is no                        **
    // other contribution to the cross section from another resonance   **
    // and an approximation for an overlay of resonances;               **
    // for decay channels other than the one-pion decay a isotropic     **
    // distribution is used                                             **
    // *******************************************************************
    // ** Date: 16/02/98   **
    // ** author: A.Muecke **
    // **********************    
    double scatterangle(int IRES, int L0);

    // *******************************************************************
    // probability distribution for scattering angle of given resonance **
    // IRES and incident nucleon L0 ;                                   **
    // z is cosine of scattering angle in CMF frame                     **
    // *******************************************************************    
    double probangle(int IRES, int L0, double z);

    // *********************************************************************
    // arbitrary Lorentz transformation
    // (taken from PHOJET v1.12, R.E. 08/98)
    // *********************************************************************
    PO_ALTRA_output PO_ALTRA(double GA, double BGX, double BGY, double BGZ, double PCX, double PCY, double PCZ, double EC);

    // *********************************************************************
    // rotation of coordinate frame (1) de rotation around y axis
    //                              (2) fe rotation around z axis
    // (taken from PHOJET v1.12, R.E. 08/98)
    // *********************************************************************
    std::vector<double> PO_TRANS(double XO, double YO, double ZO, double CDE, double SDE, double CFE, double SFE);

    // *********************************************************************
    // select x values of soft string ends using PO_RNDBET
    // (taken from PHOJET v1.12, R.E. 08/98)
    // *********************************************************************
    PO_SELSX2_output PO_SELSX2(double XMIN[2], double XMAX[2], double AS1, double AS2);
    
    // *********************************************************************
    // random number generation from beta
    // distribution in region  0 < X < 1.
    // F(X) = X**(GAM-1.)*(1.-X)**(ETA-1)*GAMM(ETA+GAM) / (GAMM(GAM*GAMM(ETA))
    // (taken from PHOJET v1.12, R.E. 08/98)
    // *********************************************************************
    double PO_RNDBET(double GAM, double ETA);
    
    // *********************************************************************
    // random number selection from gamma distribution
    // F(X) = ALAM**ETA*X**(ETA-1)*EXP(-ALAM*X) / GAM(ETA)
    // (taken from PHOJET v1.12, R.E. 08/98)
    // Note1: if you compare this function with its FORTRAN original version,
    // have a look at how GOTOs work there and how doubles are being cast into
    // integers and vice versa. The redundant code solves the GOTO pointers.
    // this was definetely one of the most exciting functions to translate!
    // Note2 ALAM would have been an argument to this function but is always called with ALAM=1.
    // *********************************************************************
    double PO_RNDGAM(double ETA);

    // interface to Lund/Jetset fragmentation (R.E. 08/98)
    void lund_frag(double SQS);

    // store initial configuration into Lund common block (R.E. 08/98)
    void lund_put(int I, int IFL, double PX, double PY, double PZ, double EE);
    
    // read final states from Lund common block (R.E. 08/98)
    lund_get_output lund_get(int I);
    
    // convert PDG particle codes to SIBYLL particle codes (R.E. 09/97)
    int ICON_PDG_SIB(int ID);
    
    //    auxiliary function for two/three particle decay mode
    //    (standard LAMBDA**(1/2) function)
    //    (taken from PHOJET 1.12, R.E. 08/98)
    double PO_XLAM(double X, double Y, double Z);

    // This is the RNG called by all other non-JETSET routines.
    // The same as RLU from JETSET except for another seed.
    double RNDM();

    double Pl(double x, double xth, double xmax, double alpha);
    double Ef(double x, double th, double w);
    double breitwigner(double sigma_0, double Gamma, double DMM, double eps_prime);
};

#endif