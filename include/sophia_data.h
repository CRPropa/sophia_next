#ifndef SOPHIA_DATA_H
#define SOPHIA_DATA_H

/* 
    SOPHIA supports 49 particles, listed in the subsequently.

    scheme: 'particle name' : (SOPHIA/SIBYLL id , PDG id)
    
        'gamma'   : (1, 22),
        'e+'      : (2, -11),
        'e-'      : (3, 11),
        'mu+'     : (4, -13),
        'mu-'     : (5, 13),
        'pi0'     : (6, 111),
        'pi+'     : (7, 211),
        'pi-'     : (8, -211),
        'K+'      : (9, 321),
        'K-'      : (10, -321),
        'K_L0'    : (11, 130),
        'K_S0'    : (12, 310),
        'p+'      : (13, 2212),
        'n0'      : (14, 2112),
        'nu_e'    : (15, 12),
        'nu_ebar' : (16, -12),
        'nu_mu'   : (17, 14),
        'nu_mubar': (18, -14),
        'pbar-'   : (19, -2212),
        'nbar0'   : (20, -2112),
        'K0'      : (21, 311),
        'Kbar0'   : (22, -311),
        'eta'     : (23, 221),
        "eta'"    : (24, 331),
        'rho+'    : (25, 213),
        'rho-'    : (26, -213),
        'rho0'    : (27, 113),
        'K*+'     : (28, 323),
        'K*-'     : (29, -323),
        'K*0'     : (30, 313),
        'K*bar0'  : (31, -313),
        'omega'   : (32, 223),
        'phi'     : (33, 333),
        'Sigma+'  : (34, 3222),
        'Sigma0'  : (35, 3212),
        'Sigma-'  : (36, 3112),
        'Xi0'     : (37, 3322),
        'Xi-'     : (38, 3312),
        'Lambda0' : (39, 3122),
        'Delta++' : (40, 2224),
        'Delta+'  : (41, 2214),
        'Delta0'  : (42, 2114),
        'Delta-'  : (43, 1114),
        'Sigma*+' : (44, 3224),
        'Sigma*0' : (45, 3214),
        'Sigma*-' : (46, 3114),
        'Xi*0'    : (47, 3324),
        'Xi*-'    : (48, 3314),
        'Omega-'  : (49, 3334),
*/



/* SOPHIA data block */



//--------------------------------------------------------------------------------------
// these arrays contain particle data for the 49 particles known to SOPHIA.
// they serve as an early dictionary: if inputting the SIBYLL particle ID,
// you get the particle data associated with it
// remember, that in FORTRAN, array numbers are being counted from 1 onwards and not 0
extern double AM[49];  // particle masses in GeV
extern int IBAR[49];   // baryon numbers
extern int ICHP[49];   // charges
extern int IDB[49];    // if entry > 0 -> particle is declared stable. all default numbers > 0 in this array are the particle IDs to remind the reader of the function of that array entry

//--------------------------------------------------------------------------------------
// used by DECPAR_nonZero only which is a function that executes particle decays (of the 49 particles avilable)
extern double CBR[102];  // containing numbers in [0,1]. I seems like these are probabilities for decay channels
extern int KDEC[612];    // length is 6x larger than length of CBR. Entries correspond to particle IDs of decay products (in units of 6 possible particles per decay). These decays are likely to occur with probability CBR
extern int LBARP[49];    // contains SIBYLL anti-particle IDs. Due to the nature of how particles are listed, some anti-particles are contained already. Also, some listed particles might not have an anti-particle

//--------------------------------------------------------------------------------------
// used RES_DECAY3 only
extern int KDECRES1p[90];
extern int KDECRES2p[180];
extern int KDECRES3p[130];
extern int KDECRES1n[90];
extern int KDECRES2n[180];
extern int KDECRES3n[110];

//--------------------------------------------------------------------------------------
// used by PROC_TWOPART only
extern int FRES[49];
extern double XLIMRES[49];

//--------------------------------------------------------------------------------------
// used by crossection only
extern double AMRESp[9];
extern double AMRESn[9];
extern double BGAMMAp[9];
extern double BGAMMAn[9];
extern double RATIOJp[9];
extern double RATIOJn[9];
extern double WIDTHp[9];
extern double WIDTHn[9];

//--------------------------------------------------------------------------------------
// used by DEC_PROC2 only
extern double CBRRES1p[18];
extern double CBRRES2p[36];
extern double CBRRES3p[26];
extern double CBRRES1n[18];
extern double CBRRES2n[36];
extern double CBRRES3n[22];
extern int IDBRES1p[9];
extern int IDBRES2p[9];
extern int IDBRES3p[9];
extern int IDBRES1n[9];
extern int IDBRES2n[9];
extern int IDBRES3n[9];
extern int ELIMITSp[9];
extern int ELIMITSn[9];
extern double RESLIMp[36];
extern double RESLIMn[36];



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