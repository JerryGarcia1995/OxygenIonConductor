**Main database: see database/oxygen_ion_conductor_dataset.csv**

Field	Description:

References --	Family name of the first author, journal title, volume, initial page number, and published year.

DOI	-- Digital Object Identifier of the source publication, enabling traceability.

Year --	Publication year of the experimental report.

Class -- Structural class of the material (e.g., perovskite, scheelite, apatite).

Formula -- Reported chemical composition.

Parsed formula -- Standardized composition represented by A_s–B_s(–C_s) site classification. For example, in the perovskite oxide {\rm La}_{0.95}{\rm Sr}_{0.05}{\rm Ga}_{0.95}{\rm Mg}_{0.05}O_{3-\delta}, La and Sr occupy the A_s-site, while Ga and Mg occupy the B_s-site of the perovskite lattice.10
E_a / A -- Activation energy in meV and prefactor in K\bulletS\bullet{\rm cm}^{-1}. When there are two distinct linear regions, values for the low-temperature regime are shown. 

E_{a,HT} / A_{HT} -- Activation energy and prefactor values for the high-temperature regime (if applicable).

T^\ast -- Transition temperature separating low- and high-temperature regimes (if applicable).

Type of \sigma_T -- Conductivity type, distinguishing bulk conductivity from total conductivity (bulk and grain boundary contributions).

Measurement -- Experimental technique employed (e.g., two-probe AC, impedance spectroscopy).

Measurement temperature range -- The range of temperatures over which experimental measurements were performed.

Source -- Origin of the data within the publication (figure or table reference).

Plot type -- Axis configuration used in the original Arrhenius representation (e.g., \log_e{\sigma_T}-\sfrac{1000}{T})

**[Reference]**

Note. Will be uploaded soon in arXiv.

**Symbolic regression modelling package: see modelling_GoodRegressor**

How to compile:

For example,

g++ GoodDesigner.cpp -o GoodDesigner.x -std=c++11

g++ GoodCurator.cpp -o GoodCurator.x -std=c++11

mpicxx GoodRegressor_slow.cpp -o GoodRegressor.x -std=c++11

Note. GoodDesigner and GoodRegressor need Eigen library (see https://libeigen.gitlab.io/): add -I/path/to/Eigen if you have your own Eigen directory.

Note. As the filename suggests, there is a faster version (~x100 fast in node-hour scale) of "GoodRegressor", which is not open yet.

Each directory contains input files and a source file (same in "CodeOnly")

001_GoodDesigner: From compositions and structures, the module draws out features (descriptor candidates).

002_GoodCurator: The module adds the first-order (simple) interactions between terms of features. 

003_GoodRegressor: The module generates symbolic regeession models. In this example, you can see 20 symbolic regression models for Ea, and 10 for Log10A.

004_GoodDesignerPost: The module generates the final stackin-ensembled models and finds the important features/interactions. On top of this, it generates LaTeX/Mathematica codes of the symbolic regression formulae for users' further post-process.

**[Reference]**

S.-H. Jang, GoodRegressor: A General-Purpose Symbolic Regression Framework for Physically Interpretable Materials Modeling, arXiv:2510.18325
