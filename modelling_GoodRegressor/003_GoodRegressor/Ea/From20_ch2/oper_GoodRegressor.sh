#!/bin/sh

#SBATCH -p F144cpu
#SBATCH -N 144
#SBATCH -n 18432
#SBATCH -c 1
#SBATCH -t 24:00:00

#mpicxx GoodRegressor.cpp -o GoodRegressor.x -std=c++11

rm -f *.pr.txt
cat > basic_properties_GoodRegressor.txt << EOF
DataFilename(TSV;ForDataNameEachLine,PrepareColumn"Name";Otherwise,Ignored)	input_oxygenionconductor_Ea.txt
OutputFilename(*_[setY].txt)	output_input_oxygenionconductor_Ea
OutputPrecision	20
tFTestMesh	1000
NVariableLimit(IfNegative,Limit)	-20 -19I -18I -17I -16I -15I -14I -13I -12I -11I -10I -9I -8I -7I -6I -5I -4I -3I -2I -1I
PickUpPLevel	0.05
Verbose(No0Yes1)	0
TargetCPUTimeBeforeSwap(sec) 1000
TrainRatio	0.8
R2Threshold 0.05
MaxBFGSIterations 10
DuplicateX 1
NameDuplicateHandling(Ignore-1;TakeMin0;TakeMax1;TakeMostHit2,WorkOnlyWhenColumn"Name"IsPrepared;If2,ProvideMinMaxIncrementForHistogram)	-1
NTestRobustness 100
setY(Name/Min/Max:IfNoLimit,Min/Max=*)	Ea	*	*
setX	OM	AveElectronegativityFromO	StdDvElectronegativityFromO	SkewElectronegativityFromO	AveMass	StdDvMass	SkewMass	AveValence	StdDvValence	SkewValence	AveRadius06	StdDvRadius06	SkewRadius06	AveXensity	StdDvXensity	SkewXensity	AveMolDensity	StdDvMolDensity	SkewMolDensity	AveFilling	StdDvFilling	SkewFilling	AveBulkModulus	StdDvBulkModulus	SkewBulkModulus	AveShearModulus	StdDvShearModulus	SkewShearModulus	AvePoissonRatio	StdDvPoissonRatio	SkewPoissonRatio	AveThermalConductivity	StdDvThermalConductivity	SkewThermalConductivity	AveThermalExpansion	StdDvThermalExpansion	SkewThermalExpansion	AveOuterPrincipal	StdDvOuterPrincipal	SkewOuterPrincipal	AveOuterAzimuthal	StdDvOuterAzimuthal	SkewOuterAzimuthal	AveDebyeTemp	StdDvDebyeTemp	SkewDebyeTemp	AveCoor	StdDvCoor	SkewCoor	AveRadius	StdDvRadius	SkewRadius	(OM*AveElectronegativityFromO)	(OM/AveElectronegativityFromO)	(OM*AveMass)	(OM/AveMass)	(OM*AveValence)	(OM/AveValence)	(OM*AveRadius06)	(OM/AveRadius06)	(OM*AveXensity)	(OM/AveXensity)	(OM*AveMolDensity)	(OM/AveMolDensity)	(OM*AveFilling)	(OM/AveFilling)	(OM*AveBulkModulus)	(OM/AveBulkModulus)	(OM*AveShearModulus)	(OM/AveShearModulus)	(OM*AvePoissonRatio)	(OM/AvePoissonRatio)	(OM*AveThermalConductivity)	(OM/AveThermalConductivity)	(OM*AveThermalExpansion)	(OM/AveThermalExpansion)	(OM*AveOuterPrincipal)	(OM/AveOuterPrincipal)	(OM*AveOuterAzimuthal)	(OM/AveOuterAzimuthal)	(OM*AveDebyeTemp)	(OM/AveDebyeTemp)	(OM*AveCoor)	(OM/AveCoor)	(OM*AveRadius)	(OM/AveRadius)	(AveElectronegativityFromO*AveMass)	(AveElectronegativityFromO/AveMass)	(AveElectronegativityFromO*AveValence)	(AveElectronegativityFromO/AveValence)	(AveElectronegativityFromO*AveRadius06)	(AveElectronegativityFromO/AveRadius06)	(AveElectronegativityFromO*AveXensity)	(AveElectronegativityFromO/AveXensity)	(AveElectronegativityFromO*AveMolDensity)	(AveElectronegativityFromO/AveMolDensity)	(AveElectronegativityFromO*AveFilling)	(AveElectronegativityFromO/AveFilling)	(AveElectronegativityFromO*AveBulkModulus)	(AveElectronegativityFromO/AveBulkModulus)	(AveElectronegativityFromO*AveShearModulus)	(AveElectronegativityFromO/AveShearModulus)	(AveElectronegativityFromO*AvePoissonRatio)	(AveElectronegativityFromO/AvePoissonRatio)	(AveElectronegativityFromO*AveThermalConductivity)	(AveElectronegativityFromO/AveThermalConductivity)	(AveElectronegativityFromO*AveThermalExpansion)	(AveElectronegativityFromO/AveThermalExpansion)	(AveElectronegativityFromO*AveOuterPrincipal)	(AveElectronegativityFromO/AveOuterPrincipal)	(AveElectronegativityFromO*AveOuterAzimuthal)	(AveElectronegativityFromO/AveOuterAzimuthal)	(AveElectronegativityFromO*AveDebyeTemp)	(AveElectronegativityFromO/AveDebyeTemp)	(AveElectronegativityFromO*AveCoor)	(AveElectronegativityFromO/AveCoor)	(AveElectronegativityFromO*AveRadius)	(AveElectronegativityFromO/AveRadius)	(AveMass*AveValence)	(AveMass/AveValence)	(AveMass*AveRadius06)	(AveMass/AveRadius06)	(AveMass*AveXensity)	(AveMass/AveXensity)	(AveMass*AveMolDensity)	(AveMass/AveMolDensity)	(AveMass*AveFilling)	(AveMass/AveFilling)	(AveMass*AveBulkModulus)	(AveMass/AveBulkModulus)	(AveMass*AveShearModulus)	(AveMass/AveShearModulus)	(AveMass*AvePoissonRatio)	(AveMass/AvePoissonRatio)	(AveMass*AveThermalConductivity)	(AveMass/AveThermalConductivity)	(AveMass*AveThermalExpansion)	(AveMass/AveThermalExpansion)	(AveMass*AveOuterPrincipal)	(AveMass/AveOuterPrincipal)	(AveMass*AveOuterAzimuthal)	(AveMass/AveOuterAzimuthal)	(AveMass*AveDebyeTemp)	(AveMass/AveDebyeTemp)	(AveMass*AveCoor)	(AveMass/AveCoor)	(AveMass*AveRadius)	(AveMass/AveRadius)	(AveValence*AveRadius06)	(AveValence/AveRadius06)	(AveValence*AveXensity)	(AveValence/AveXensity)	(AveValence*AveMolDensity)	(AveValence/AveMolDensity)	(AveValence*AveFilling)	(AveValence/AveFilling)	(AveValence*AveBulkModulus)	(AveValence/AveBulkModulus)	(AveValence*AveShearModulus)	(AveValence/AveShearModulus)	(AveValence*AvePoissonRatio)	(AveValence/AvePoissonRatio)	(AveValence*AveThermalConductivity)	(AveValence/AveThermalConductivity)	(AveValence*AveThermalExpansion)	(AveValence/AveThermalExpansion)	(AveValence*AveOuterPrincipal)	(AveValence/AveOuterPrincipal)	(AveValence*AveOuterAzimuthal)	(AveValence/AveOuterAzimuthal)	(AveValence*AveDebyeTemp)	(AveValence/AveDebyeTemp)	(AveValence*AveCoor)	(AveValence/AveCoor)	(AveValence*AveRadius)	(AveValence/AveRadius)	(AveRadius06*AveXensity)	(AveRadius06/AveXensity)	(AveRadius06*AveMolDensity)	(AveRadius06/AveMolDensity)	(AveRadius06*AveFilling)	(AveRadius06/AveFilling)	(AveRadius06*AveBulkModulus)	(AveRadius06/AveBulkModulus)	(AveRadius06*AveShearModulus)	(AveRadius06/AveShearModulus)	(AveRadius06*AvePoissonRatio)	(AveRadius06/AvePoissonRatio)	(AveRadius06*AveThermalConductivity)	(AveRadius06/AveThermalConductivity)	(AveRadius06*AveThermalExpansion)	(AveRadius06/AveThermalExpansion)	(AveRadius06*AveOuterPrincipal)	(AveRadius06/AveOuterPrincipal)	(AveRadius06*AveOuterAzimuthal)	(AveRadius06/AveOuterAzimuthal)	(AveRadius06*AveDebyeTemp)	(AveRadius06/AveDebyeTemp)	(AveRadius06*AveCoor)	(AveRadius06/AveCoor)	(AveRadius06*AveRadius)	(AveRadius06/AveRadius)	(AveXensity*AveMolDensity)	(AveXensity/AveMolDensity)	(AveXensity*AveFilling)	(AveXensity/AveFilling)	(AveXensity*AveBulkModulus)	(AveXensity/AveBulkModulus)	(AveXensity*AveShearModulus)	(AveXensity/AveShearModulus)	(AveXensity*AvePoissonRatio)	(AveXensity/AvePoissonRatio)	(AveXensity*AveThermalConductivity)	(AveXensity/AveThermalConductivity)	(AveXensity*AveThermalExpansion)	(AveXensity/AveThermalExpansion)	(AveXensity*AveOuterPrincipal)	(AveXensity/AveOuterPrincipal)	(AveXensity*AveOuterAzimuthal)	(AveXensity/AveOuterAzimuthal)	(AveXensity*AveDebyeTemp)	(AveXensity/AveDebyeTemp)	(AveXensity*AveCoor)	(AveXensity/AveCoor)	(AveXensity*AveRadius)	(AveXensity/AveRadius)	(AveMolDensity*AveFilling)	(AveMolDensity/AveFilling)	(AveMolDensity*AveBulkModulus)	(AveMolDensity/AveBulkModulus)	(AveMolDensity*AveShearModulus)	(AveMolDensity/AveShearModulus)	(AveMolDensity*AvePoissonRatio)	(AveMolDensity/AvePoissonRatio)	(AveMolDensity*AveThermalConductivity)	(AveMolDensity/AveThermalConductivity)	(AveMolDensity*AveThermalExpansion)	(AveMolDensity/AveThermalExpansion)	(AveMolDensity*AveOuterPrincipal)	(AveMolDensity/AveOuterPrincipal)	(AveMolDensity*AveOuterAzimuthal)	(AveMolDensity/AveOuterAzimuthal)	(AveMolDensity*AveDebyeTemp)	(AveMolDensity/AveDebyeTemp)	(AveMolDensity*AveCoor)	(AveMolDensity/AveCoor)	(AveMolDensity*AveRadius)	(AveMolDensity/AveRadius)	(AveFilling*AveBulkModulus)	(AveFilling/AveBulkModulus)	(AveFilling*AveShearModulus)	(AveFilling/AveShearModulus)	(AveFilling*AvePoissonRatio)	(AveFilling/AvePoissonRatio)	(AveFilling*AveThermalConductivity)	(AveFilling/AveThermalConductivity)	(AveFilling*AveThermalExpansion)	(AveFilling/AveThermalExpansion)	(AveFilling*AveOuterPrincipal)	(AveFilling/AveOuterPrincipal)	(AveFilling*AveOuterAzimuthal)	(AveFilling/AveOuterAzimuthal)	(AveFilling*AveDebyeTemp)	(AveFilling/AveDebyeTemp)	(AveFilling*AveCoor)	(AveFilling/AveCoor)	(AveFilling*AveRadius)	(AveFilling/AveRadius)	(AveBulkModulus*AveShearModulus)	(AveBulkModulus/AveShearModulus)	(AveBulkModulus*AvePoissonRatio)	(AveBulkModulus/AvePoissonRatio)	(AveBulkModulus*AveThermalConductivity)	(AveBulkModulus/AveThermalConductivity)	(AveBulkModulus*AveThermalExpansion)	(AveBulkModulus/AveThermalExpansion)	(AveBulkModulus*AveOuterPrincipal)	(AveBulkModulus/AveOuterPrincipal)	(AveBulkModulus*AveOuterAzimuthal)	(AveBulkModulus/AveOuterAzimuthal)	(AveBulkModulus*AveDebyeTemp)	(AveBulkModulus/AveDebyeTemp)	(AveBulkModulus*AveCoor)	(AveBulkModulus/AveCoor)	(AveBulkModulus*AveRadius)	(AveBulkModulus/AveRadius)	(AveShearModulus*AvePoissonRatio)	(AveShearModulus/AvePoissonRatio)	(AveShearModulus*AveThermalConductivity)	(AveShearModulus/AveThermalConductivity)	(AveShearModulus*AveThermalExpansion)	(AveShearModulus/AveThermalExpansion)	(AveShearModulus*AveOuterPrincipal)	(AveShearModulus/AveOuterPrincipal)	(AveShearModulus*AveOuterAzimuthal)	(AveShearModulus/AveOuterAzimuthal)	(AveShearModulus*AveDebyeTemp)	(AveShearModulus/AveDebyeTemp)	(AveShearModulus*AveCoor)	(AveShearModulus/AveCoor)	(AveShearModulus*AveRadius)	(AveShearModulus/AveRadius)	(AvePoissonRatio*AveThermalConductivity)	(AvePoissonRatio/AveThermalConductivity)	(AvePoissonRatio*AveThermalExpansion)	(AvePoissonRatio/AveThermalExpansion)	(AvePoissonRatio*AveOuterPrincipal)	(AvePoissonRatio/AveOuterPrincipal)	(AvePoissonRatio*AveOuterAzimuthal)	(AvePoissonRatio/AveOuterAzimuthal)	(AvePoissonRatio*AveDebyeTemp)	(AvePoissonRatio/AveDebyeTemp)	(AvePoissonRatio*AveCoor)	(AvePoissonRatio/AveCoor)	(AvePoissonRatio*AveRadius)	(AvePoissonRatio/AveRadius)	(AveThermalConductivity*AveThermalExpansion)	(AveThermalConductivity/AveThermalExpansion)	(AveThermalConductivity*AveOuterPrincipal)	(AveThermalConductivity/AveOuterPrincipal)	(AveThermalConductivity*AveOuterAzimuthal)	(AveThermalConductivity/AveOuterAzimuthal)	(AveThermalConductivity*AveDebyeTemp)	(AveThermalConductivity/AveDebyeTemp)	(AveThermalConductivity*AveCoor)	(AveThermalConductivity/AveCoor)	(AveThermalConductivity*AveRadius)	(AveThermalConductivity/AveRadius)	(AveThermalExpansion*AveOuterPrincipal)	(AveThermalExpansion/AveOuterPrincipal)	(AveThermalExpansion*AveOuterAzimuthal)	(AveThermalExpansion/AveOuterAzimuthal)	(AveThermalExpansion*AveDebyeTemp)	(AveThermalExpansion/AveDebyeTemp)	(AveThermalExpansion*AveCoor)	(AveThermalExpansion/AveCoor)	(AveThermalExpansion*AveRadius)	(AveThermalExpansion/AveRadius)	(AveOuterPrincipal*AveOuterAzimuthal)	(AveOuterPrincipal/AveOuterAzimuthal)	(AveOuterPrincipal*AveDebyeTemp)	(AveOuterPrincipal/AveDebyeTemp)	(AveOuterPrincipal*AveCoor)	(AveOuterPrincipal/AveCoor)	(AveOuterPrincipal*AveRadius)	(AveOuterPrincipal/AveRadius)	(AveOuterAzimuthal*AveDebyeTemp)	(AveOuterAzimuthal/AveDebyeTemp)	(AveOuterAzimuthal*AveCoor)	(AveOuterAzimuthal/AveCoor)	(AveOuterAzimuthal*AveRadius)	(AveOuterAzimuthal/AveRadius)	(AveDebyeTemp*AveCoor)	(AveDebyeTemp/AveCoor)	(AveDebyeTemp*AveRadius)	(AveDebyeTemp/AveRadius)	(AveCoor*AveRadius)	(AveCoor/AveRadius)
SubFilterX(IfNone,"-";{Name/Min/Max:IfNoLimit,Min/Max=*}*NSubFilterX)	-
Functypes
1 ^1 pow a= 0 b= 1
2 ^-3 pow a= 0 b= -3
3 ^-2 pow a= 0 b= -2
4 ^-1 pow a= 0 b= -1
5 ^-0.5 pow a= 0 b= -0.5
6 ^-0.333333 pow a= 0 b= -0.333333
7 ^0.333333 pow a= 0 b= 0.333333
8 ^0.5 pow a= 0 b= 0.5
9 ^2 pow a= 0 b= 2
10 ^3 pow a= 0 b= 3
11 log10 log a= 0 b= 10 c= 1
12 log10^-3 log a= 0 b= 10 c= -3
13 log10^-2 log a= 0 b= 10 c= -2
14 log10^-1 log a= 0 b= 10 c= -1
15 log10^2 log a= 0 b= 10 c= 2
16 log10^3 log a= 0 b= 10 c= 3
17 10^1 exp a= 0 b= 1 c= 1 d= 10 e= 1
18 10^-3 exp a= 0 b= 1 c= 1 d= 10 e= -3
19 10^-2 exp a= 0 b= 1 c= 1 d= 10 e= -2
20 10^-1 exp a= 0 b= 1 c= 1 d= 10 e= -1
21 10^2 exp a= 0 b= 1 c= 1 d= 10 e= 2
22 10^3 exp a= 0 b= 1 c= 1 d= 10 e= 3
23 exp^1 exp a= 0 b= 1 c= 1 d= 2.71828 e= 1
24 exp^-3 exp a= 0 b= 1 c= 1 d= 2.71828 e= -3
25 exp^-2 exp a= 0 b= 1 c= 1 d= 2.71828 e= -2
26 exp^-1 exp a= 0 b= 1 c= 1 d= 2.71828 e= -1
27 exp^2 exp a= 0 b= 1 c= 1 d= 2.71828 e= 2
28 exp^3 exp a= 0 b= 1 c= 1 d= 2.71828 e= 3
29 erf erf a= 0 b= 1 c= 1
30 erf-3 erf a= 0 b= 0.001 c= 1
31 erf-2 erf a= 0 b= 0.01 c= 1
32 erf-1 erf a= 0 b= 0.1 c= 1
33 erf1 erf a= 0 b= 10 c= 1
34 erf2 erf a= 0 b= 100 c= 1
35 erf3 erf a= 0 b= 1000 c= 1
36 erf^2 erf a= 0 b= 1 c= 2
37 erf-3^2 erf a= 0 b= 0.001 c= 2
38 erf-2^2 erf a= 0 b= 0.01 c= 2
39 erf-1^2 erf a= 0 b= 0.1 c= 2
40 erf1^2 erf a= 0 b= 10 c= 2
41 erf2^2 erf a= 0 b= 100 c= 2
42 erf3^2 erf a= 0 b= 1000 c= 2
43 erf^3 erf a= 0 b= 1 c= 3
44 erfoffset-0.05 erf a= -0.05 b= 1 c= 1
45 erfoffset-0.1 erf a= -0.1 b= 1 c= 1
46 erfoffset-0.5 erf a= -0.5 b= 1 c= 1
47 erfoffset-0.875 erf a= -0.875 b= 1 c= 1
48 erfoffset-1 erf a= -1 b= 1 c= 1
49 erfoffset-5 erf a= -5 b= 1 c= 1
50 erfoffset-10 erf a= -10 b= 1 c= 1
51 erfoffset-50 erf a= -50 b= 1 c= 1
52 erfoffset-100 erf a= -100 b= 1 c= 1
53 erfoffset-1000 erf a= -1000 b= 1 c= 1
54 erfoffset-10000 erf a= -10000 b= 1 c= 1
55 erfoffset-0.05^2 erf a= -0.05 b= 1 c= 2
56 erfoffset-0.1^2 erf a= -0.1 b= 1 c= 2
57 erfoffset-0.5^2 erf a= -0.5 b= 1 c= 2
58 erfoffset-0.875^2 erf a= -0.875 b= 1 c= 2
59 erfoffset-1^2 erf a= -1 b= 1 c= 2
60 erfoffset-5^2 erf a= -5 b= 1 c= 2
61 erfoffset-10^2 erf a= -10 b= 1 c= 2
62 erfoffset-50^2 erf a= -50 b= 1 c= 2
63 erfoffset-100^2 erf a= -100 b= 1 c= 2
64 erfoffset-1000^2 erf a= -1000 b= 1 c= 2
65 erfoffset-10000^2 erf a= -10000 b= 1 c= 2
66 erfoffset-0.875^3 erf a= -0.875 b= 1 c= 3
67 erfoffset-0.875^4 erf a= -0.875 b= 1 c= 4
68 sin sin a= 0 b= 1 c= 1 d= 1
69 sinhpi sin a= 0 b= 1.57079632679 c= 1 d= 1
70 sinpi sin a= 0 b= 3.14159265359 c= 1 d= 1
71 sin-3 sin a= 0 b= 0.001 c= 1 d= 1
72 sinhpi-3 sin a= 0 b= 0.00157079632679 c= 1 d= 1
73 sinpi-3 sin a= 0 b= 0.00314159265359 c= 1 d= 1
74 sin-2 sin a= 0 b= 0.01 c= 1 d= 1
75 sinhpi-2 sin a= 0 b= 0.0157079632679 c= 1 d= 1
76 sinpi-2 sin a= 0 b= 0.0314159265359 c= 1 d= 1
77 sin-1 sin a= 0 b= 0.1 c= 1 d= 1
78 sinhpi-1 sin a= 0 b= 0.157079632679 c= 1 d= 1
79 sinpi-1 sin a= 0 b= 0.314159265359 c= 1 d= 1
80 sin1 sin a= 0 b= 10 c= 1 d= 1
81 sinhpi1 sin a= 0 b= 15.7079632679 c= 1 d= 1
82 sinpi1 sin a= 0 b= 31.4159265359 c= 1 d= 1
83 sin2 sin a= 0 b= 100 c= 1 d= 1
84 sinhpi2 sin a= 0 b= 157.079632679 c= 1 d= 1
85 sinpi2 sin a= 0 b= 314.159265359 c= 1 d= 1
86 sin3 sin a= 0 b= 1000 c= 1 d= 1
87 sinhpi3 sin a= 0 b= 1570.79632679 c= 1 d= 1
88 sinpi3 sin a= 0 b= 3141.59265359 c= 1 d= 1
89 cos cos a= 0 b= 1 c= 1 d= 1
90 coshpi cos a= 0 b= 1.57079632679 c= 1 d= 1
91 cospi cos a= 0 b= 3.14159265359 c= 1 d= 1
92 cos-3 cos a= 0 b= 0.001 c= 1 d= 1
93 coshpi-3 cos a= 0 b= 0.00157079632679 c= 1 d= 1
94 cospi-3 cos a= 0 b= 0.00314159265359 c= 1 d= 1
95 cos-2 cos a= 0 b= 0.01 c= 1 d= 1
96 coshpi-2 cos a= 0 b= 0.0157079632679 c= 1 d= 1
97 cospi-2 cos a= 0 b= 0.0314159265359 c= 1 d= 1
98 cos-1 cos a= 0 b= 0.1 c= 1 d= 1
99 coshpi-1 cos a= 0 b= 0.157079632679 c= 1 d= 1
100 cospi-1 cos a= 0 b= 0.314159265359 c= 1 d= 1
101 cos1 cos a= 0 b= 10 c= 1 d= 1
102 coshpi1 cos a= 0 b= 15.7079632679 c= 1 d= 1
103 cospi1 cos a= 0 b= 31.4159265359 c= 1 d= 1
104 cos2 cos a= 0 b= 100 c= 1 d= 1
105 coshpi2 cos a= 0 b= 157.079632679 c= 1 d= 1
106 cospi2 cos a= 0 b= 314.159265359 c= 1 d= 1
107 cos3 cos a= 0 b= 1000 c= 1 d= 1
108 coshpi3 cos a= 0 b= 1570.79632679 c= 1 d= 1
109 cospi3 cos a= 0 b= 3141.59265359 c= 1 d= 1
BetaRegression(IfNo,"none")	none
logit probit cloglog cauchit nloglog

####functype####
pow	(x+a)^b
log	(log_b (x+a))^c
sin	(sin(((b(x+a))^c))^d
cos	(cos(((b(x+a))^c))^d
tan	(tan(((b(x+a))^c))^d
sinh	(sinh(((b(x+a))^c))^d
cosh	(cosh(((b(x+a))^c))^d
tanh	(tanh(((b(x+a))^c))^d
exp	(d^(c*((x+a)^b)))^e
erf	erf(b(x+a))^c
erfc	erf(b(x+a))^c
abs	|x+a|^b

####betatype####
logit probit cloglog cauchit nloglog

EOF

srun ./GoodRegressor.x 
