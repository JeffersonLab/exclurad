#Exclurad
Fortran code for calculation of radiative corrections to exclusive electroproduction of pions on a nucleon. Current version computes corrections to the unpolarized coincidence cross section and beam asymmetry (aka the fifth structure function). Distinctive features are 
a) Covariant technique of cancellation of the infrared divergence leads to independence of the parameter that splits soft and hard regions of brem photons and 
b) Integration over the brem photon phase space is exact, not relying on the peaking approximation. 
MAID and AO are used to model the reaction mechanism. The code is extendable to any exclusive electron scattering with two-body breakup, such as p(e,e'K)Lambda, d(e,e'p)n, 3He(e,e'p)d, etc.

##References
```
QED RADIATIVE CORRECTIONS IN PROCESSES OF EXCLUSIVE PION ELECTROPRODUCTION.
By A. Afanasev (Jefferson Lab), I. Akushevich (Duke U.), V. Burkert, K. Joo (Jefferson Lab) 
Published in Phys.Rev.D66:074004,2002
e-Print Archive: hep-ph/0208183 
```
## Links
https://arxiv.org/pdf/hep-ph/0208183.pdf

# Getting Started

## Clone the repository
```
git clone https://github.com/JeffersonLab/exclurad.git
```

## Build the executable
```
cd exclurad/exclurad
scons-2.7

```
## Example input file (input.dat):

```
3       !  1: AO 2: maid98  3: maid2007
0       !  0: Full, 1: Factorizable and Leading log  
1.645   !  bmom - lepton momentum
0.0     !  tmom - momentum per nucleon
1       !  lepton - 1 electron, 2 muon
2       !  ivec - detected hadron (1) p, (2) pi+
0.05      !  vcut - cut on inelasticity (0.) if no cut, negative -- v


10 ! no. of points
1.232 1.232 1.232 1.232 1.232 1.232 1.232 1.232 1.232 1.232 ! W values 
0.4   0.4   0.4   0.4   0.4   0.4   0.4   0.4   0.4   0.4 ! Q^2 values
0.    0.    0.    0.    0.    0.    0.    0.    0.    0. ! Cos(Theta)
10.   50.   90.   130.  170.  190.  230.  270.  310.  350. ! phi values

```

## Example usage:
```
./build/exclurad.exe < input.dat

```

