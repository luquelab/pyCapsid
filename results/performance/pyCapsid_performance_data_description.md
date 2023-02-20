# File descriptors
The list below describes the columns in `results_table_toni_wip.csv`:
+ `PDB`: Protein Data Bank (PDB) identification, for example, 7kq5, unless indicated otherwise (for example, an entry from the electron microscopy database, like EMD-21123).
+ `Capsid`: Name associated with the capsid or protein shell.
+ `Abbreviation`: Abbreviated version of the capsid name.
+ `Residues`: Number of amino acid residues in the capsid.
+ `MCPs`: Number of major capsid proteins (MCPs) in the capsid.
+ `CPs`: Number of capsid proteins (CPs) including MCPs and minor and auxiliary proteins (or reinforcement proteins).
+ `Reinforcement`: Axes of symmetry were the reinforcement proteins are placed.
+ `Resolution`: Resolution of empirical capsid reconstruction. The units are in angstroms (Å).
+ `T-number`: T-number associated at the distribution of hexamers and pentamers (not including the variations associated with minor proteins).
+ `Spring_constant`: Value of the spring constant obtained by fitting the normal mode analysis of the capsid elastic model to the thermal fluctuations (B-factors) of the empirical capsid reconstruction. The units are in kT per Å^2.
+ `Spring_stderr`: Standard error associated with the fitted value of the spring constant.
+ `Spring_R^2`: Correlation coefficient associated with the fitting leading to the spring constant values.
+ `CC_uENM`: Corrlelation coefficient between the empirical and the simulated B-factors using the unified elasticn networ model (uENM).
+ `CC_bGNM`: Corrlelation coefficient between the empirical and the simulated B-factors using the beta Gaussian network model (bGNM).
+ `CC_ANM`: Corrlelation coefficient between the empirical and the simulated B-factors using the anysotropic network model (ANM).
+ `Clusters`: Number of quasi-rigid clusters obtained for the uENM model.
+ `Outliers`: Number of amino acid residues not assigned to clusters.
+ `CQS`: Clusteer quality score obtained for the best selected capsid clustering using the uENM model.
+ `Icos_class`: Icosahedral classification of the distribution of quasi-rigid clusters obtained in the capsid.
+ `Memory`: Memory usage obtained in the computational performance benchmarking. Units: Mebibytes (MiB), that is, 2^20 bytes. It is equal to 2^20/10^6 MegaBytes (MB) (Conversion factor: 1.048576).
+ `Runtime`: Runtime obtained in the computational performance benchmarking. Units: seconds (sec).
+ `pyCap_bench`: Runtime?? obtained for pyCapsid in the benchmarking with ProDy. Units: seconds (sec)?? <font color = 'red'> Confirm! </font>
+ `ProDy_bench`: Runtime?? obtained for ProDy in for benchmarking. Units: seconds (sec)?? <font color = 'red'> Confirm! </font>
