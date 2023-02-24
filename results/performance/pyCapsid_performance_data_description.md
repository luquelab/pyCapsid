# General information
This file provides complementary information for the data file `pyCapsid_performance_data.csv`:

## Description of columns
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
+ `pyCap_runtime`: Runtime obtained for pyCapsid in the benchmarking with ProDy. Units: seconds (sec).
+ `ProDy_runtime`: Runtime obtained for ProDy for benchmarking. Units: seconds (sec).
+ `pyCap_memory`: Memory peak obtained for pyCapsid in the benchmarking with ProDy. Units: Mebibytes (MiB).
+ `ProDy_memory`: Runtime obtained for ProDy for benchmarking. Units: seconds (sec). Units: Mebibytes (MiB).

## Notes
+ For the general performanc of pyCapsid, the unified elastic network model (UENM) was used. This is the default model in pyCapsid due to its otpimal performance.
+ For benchmarking with ProDy the model used was the anysotropic network model (ANM) because ProDy does not offer UENM as option. Both pyCapsid and ProDy offer ANM as on option.
+ The differences between the runtime and memory of pyCapsid in the general performance versus the benchmarking performance in the comparison with ProDy are associated with using different elastic models and focusing on the NMA results (no quasi-rigid domain decomposition, which is not available in ProDy).
+ The computational performance was estimated in a HPC core (CSRC cinci cluster) with Intel Xeon CPU E5-2650 v4 at 2.20GHz CPUs and 128GB RAM.



