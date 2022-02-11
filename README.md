# Getting Started

## clone the repository
git clone https://github.com/JeffersonLab/exclurad.git

## build the executable
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

