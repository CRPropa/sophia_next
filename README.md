# SOPHIA++

This is a C++ translation of the original SOPHIA code. SOPHIA stands for Simulations of Photo-Hadronic Interactions in Astrophysics and was published by A. Mücke et al in 1999: https://arxiv.org/pdf/astro-ph/9903478. SOPHIA in return makes use of JETSET v7.4 written by Torbjorn Sjostrand https://arxiv.org/abs/hep-ph/9508391 which thus with this distribution is, for its most part, translated to C++ as well.

## Why translate such an old code?
There are a variety of reasons why a C++ version of SOPHIA might be of interest:
* SOPHIA still is well-used, well-"proven" and run on many systems. It is known to reliably simulate nucleon/gamma-interactions.
* Everyone interested in the topic should have easy, open, and central access to SOPHIA.
* You may not want to have to make a call to the original FORTRAN77 code each time you generate an event.
* You may wish to run SOPHIA (or your project embedding it) on multiple CPUs. At this point your FORTRAN version of SOPHIA would be a 1-CPU bottleneck. This version aims to become suitable for multi-threading and will do so quite soon.
* The old, FORTRAN version of SOPHIA features a stand-alone random number generator. However, it does not actually generate random numbers but rather sample deterministic values from an internal number generation procedure thus making the code seemingly random but globally deterministic. When out of the prototype phase, this version will invoke actually random, SOPHIA-unrelated techniques for random number generation.

## What is the status of this project?
Actually, it's ready to use! However, a lot of the internal code needs more modernisation in every regard and can not be multi-threaded for the time being. This is work in progress, thus the code shall become more expressive and maybe be improved by further details some day.

## How easy is it to get running?
Super straight-forward on your Linux system and with no special dependencies. Just clone this repository to ```yourDirectory``` and in it run in your terminal
```
mkdir build
cd build
cmake ..
make
```
you can now modify the execute_sophia.cpp file to your needs. Make sure you run ```make``` again before running the executable via ```./run_sophia```.

## Does this code produce the "right" results?
This is a univeristy-world product and, to be fair, we can not exclude that bugs might occur and try to offer support on any issue that may arise. However, there are two kinds of approaches by which we aim to prove the correctness of this code by:
1) comparisons to the original SOPHIA code
2) tests against modern literature

The former is achieved by using an exact C++ copy of the "deterministic", internal, original SOPHIA random number generator and have it applied to this code. If the codes are supposed to be identical, then they should in each calculation step produce and pass the exact same numbers. This is the case for at least each 100.000 nucleon-gamma events along a large variety of relevant input parameter constellations. While running these parameter constellations, the runtime was measured on an average laptop system:
![SOPHIA runtime performance](https://github.com/CRPropa/sophia_next/blob/master/test/SOPHIA_performance.png)

The latter is achieved by comparing the energy spectra of produced secondary particles to peer-reviewed literature of people who have achieved their energy spectra from the FORTRAN version of SOPHIA. A very useful paper was published by Kelner et al. in 2008 https://arxiv.org/pdf/0803.0688.pdf. The results for the most meaningful secondary particles are:
![SOPHIA_secondarySpectra](https://github.com/CRPropa/sophia_next/blob/master/test/SOPHIA%2B%2B_spectra.png)

## Acknowledgement
Thanks a lot to those who contribute to this repository, thus helping to make science more open.

### About the author
I am a PhD student conducting research on high-energy astrophysics at Bochum University, Germany and Oxford University, UK. I also develop software for the CRPropa3 cosmic ray propagation framework https://github.com/CRPropa/CRPropa3. I value open research over local knowledge and free access over paywalls because I see the most potential for progress in the empowerment of the individual. In this sence, if at least one person finds this repository useful, I consider this project a success.
