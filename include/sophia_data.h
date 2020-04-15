#ifndef SOPHIA_DATA_H
#define SOPHIA_DATA_H

/* SOPHIA data block */

// block data DATDEC
extern int FRES[49];
extern double XLIMRES[49];
extern double AMRESp[9];
extern double AMRESn[9];
extern int IDBRES1p[9];
extern int IDBRES2p[9];
extern int IDBRES3p[9];
extern int IDBRES1n[9];
extern int IDBRES2n[9];
extern int IDBRES3n[9];
extern double CBRRES1p[18];
extern double CBRRES2p[36];
extern double CBRRES3p[26];
extern double CBRRES1n[18];
extern double CBRRES2n[36];
extern double CBRRES3n[22];
extern int KDECRES1p[90];
extern int KDECRES2p[180];
extern int KDECRES3p[130];
extern int KDECRES1n[90];
extern int KDECRES2n[180];
extern int KDECRES3n[110];
extern double RESLIMp[36];
extern double RESLIMn[36];
extern int ELIMITSp[9];
extern int ELIMITSn[9];
extern double BGAMMAp[9];
extern double RATIOJp[9];
extern double WIDTHp[9];
extern double BGAMMAn[9];
extern double RATIOJn[9];
extern double WIDTHn[9];
extern double CBR[102];
extern double AM[49];
extern double AM2[49];

/* 
    The subsequent array IDB is used by SOPHIA and not JETSET.
    It declares particles stable (entry >= 1) or not (entry = 0), if applicable.
    Remember, that in FORTRAN, array numbers are being counted from 1 onwards
    (as opposed to from 0 onwards in modern languages).
    Hence, by setting: IDB[5] = 0, IDB[6] = 0, IDB[7] = 0 declares all pions stable.

    scheme: *** 'particle name': (SOPHIA/SIBYLL id , PDG id) ***
    
        'gamma': (1,22),
        'e+': (2,-11),
        'e-': (3,11),
        'mu+': (4,-13),
        'mu-': (5,13),
        'pi0': (6,111),
        'pi+': (7,211),
        'pi-': (8,-211),
        'K+': (9,321),
        'K-': (10,-321),
        'K_L0': (11,130),
        'K_S0': (12,310),
        'p+': (13,2212),
        'n0': (14,2112),
        'nu_e': (15,12),
        'nu_ebar': (16,-12),
        'nu_mu': (17,14),
        'nu_mubar': (18,-14),
        'pbar-': (19,-2212),
        'nbar0': (20,-2112),
        'K0': (21,311),
        'Kbar0': (22,-311),
        'eta': (23,221),
        "eta'": (24,331),
        'rho+': (25,213),
        'rho-': (26,-213),
        'rho0': (27,113),
        'K*+': (28,323),
        'K*-': (29,-323),
        'K*0': (30,313),
        'K*bar0': (31,-313),
        'omega': (32,223),
        'phi': (33,333),
        'Sigma+': (34,3222),
        'Sigma0': (35,3212),
        'Sigma-': (36,3112),
        'Xi0': (37,3322),
        'Xi-': (38,3312),
        'Lambda0': (39,3122),
        'Delta++': (40,2224),
        'Delta+': (41,2214),
        'Delta0': (42,2114),
        'Delta-': (43,1114),
        'Sigma*+': (44,3224),
        'Sigma*0': (45,3214),
        'Sigma*-': (46,3114),
        'Xi*0': (47,3324),
        'Xi*-': (48,3314),
        'Omega-': (49,3334),
*/

extern int IDB[49];
extern int KDEC[612];
extern int LBARP[49];
extern int ICHP[49];
extern int ISTR[49];
extern int IBAR[49];

// block data PARAM_INI: This block data contains default values of the parameters used in fragmentation
// parameters of flavor formation
extern double PAR[8];



/* JETSET data block */



// block data LUDATA. Purpose: to give default values to parameters and particle and decay data. 

// LUDAT1, containing status codes and most parameters. 
extern int MSTU[200];
extern double PARU[200];
extern int MSTJ[200];
extern double PARJ[200];

// LUDAT2, with particle data and flavour treatment parameters.
extern int KCHG[3][500];
extern double PMAS[4][500];
extern double PARF[2000];

// LUDAT3, with particle decay parameters and data.
extern int MDCY[3][500];
extern int MDME[2][2000];
extern double BRAT[2000];
extern int KFDP[5][2000];

#endif