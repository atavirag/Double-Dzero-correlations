# Double-Dzero-correlations
Macros to run the double D⁰ correlations analysis. 

Ingredients
- 
- **AnalysisResults.root** files containing the results of task-d0 and correlator-d-meson-pairs for data and MC.
- **AO2D.root** file containing o2d0pair tree (data) and also O2d0pairmcinfo (MC).

Analysis steps
- 
1. *Load tree using **loadInvMass.cxx**. Ensure that fname contains the data's AO2D (line 49), fname_MC contains the MC's AO2D (line 52) and that the directory's name is correct in line 301. It will give a file called **fitFromTree.root**as an output.* **NOT NEEDED IF WE USE THE FILES FROM STEP 2**
2. We will need to perform a fit using all candidates that pass in our correlator to get the values of the parameters we will fix in the 2D fit. We will start by using **SeparateMassInPtRanges.cxx** for the data and **SeparateMassInPtRangesMc.cxx** for the data and MC, respectively. Both files can be modified using **config-SeparateMassInPtRanges.json**. As the name suggests, this will separate the invariant-mass plots in the correlator in our desired pT ranges.
3. Once our invariant masses are separated in pT, we can use the **HFInvMassFitter class**, run through **runMassFitter.C**, to get the 1D result. The fitter version uploaded here was slightly modified to produce a workspace containing values such as the mean, the sigma, the tau, etc., all of which will be used later.
4. If the MC is new, it is important to create the efficiency map using **efficiencyMap.cxx**. It is run using **config-efficiencyMap.json** and only needs the AnalysisResults file of the MC as input. It produces ***Eff_times_Acc_Map.root**.
5. Now, we have all the necessary files to run the 2D fitter. It contains a class, **InvMassFitter2D**, and a running file, **runInvMassFitter2D.cpp**. Again, make sure that the data's AO2D and the tree directory are correct. Give it all the previous files, and run!

And that's it :)
