# Double-Dzero-correlations
Macros and instructions to run the double D⁰ correlations analysis. More info and first results can be found in this presentation: https://indico.cern.ch/event/1474633/contributions/6209450/attachments/2964901/5215937/DDMeasurement_update_11-11-24.pdf.

## Ingredients
- **AnalysisResults.root** files containing the results of task-d0 and correlator-d-meson-pairs for data and MC.
- **AO2D.root** file containing o2d0pair tree (data) and also O2d0pairmcinfo (MC).

## BDT trainings
TODO: add info about trainings here

## Analysis steps

### 1. Get the correct pT ranges
We will start by using ``SeparateMassInPtRanges.cxx`` and ``SeparateMassInPtRangesMc.cxx`` for the data and MC, respectively. Both files can be modified using the same config file, called ``config-SeparateMassInPtRanges.json``. As the name suggests, this will separate the invariant-mass plots in the correlator in our desired pT ranges. For now, this analysis is meant to be performed in one single pT bin, ranging from 1 to 24 GeV/c.

### 2. Perform the 1D fit
Once our invariant masses are separated in pT, we can use the ``HFInvMassFitter`` class, run through ``runMassFitter.C`` using ``config_massfitter.json``, to get the 1D result. The fitter version uploaded here was slightly modified to produce a workspace containing values such as the mean, the sigma, the tau, etc., all of which will be used later. This workspace is stored in a file called ``workspace_massfitter.root``.

**NOTE:** the workspace file will be created in the folder where the code is stored. Remember to move it/change the name manually if needed!!

### 3. (Only if needed) Get the efficiency map
 If the MC is new, it is important to create an efficiency map. Here, there are two ways to proceed, depending on how old the AnalysisResults file is:
 - For older files, we will use ``efficiencyMap.cxx``. It is run using ``config-efficiencyMap.json`` and only needs the AnalysisResults file of the MC as input. It produces ``Eff_times_Acc_Map.root``. It divides hPtVsYReco and hPtVsYGen
 - For newer files, created after the merging of the PR #8452, we will need to use ``efficiencyMap_PVContrib.cxx`` (config file: ``config-efficiencyMap_PVContrib.json``). In this case, we will use a TH3F as input, hPtVsYVsNContrib. In this case, the info from the nContrib axis (the number of primary vertex contributors) is used as a weight before computing the division of the hPtVsYReco and hPtVsYGen projections.

### 4. (Only if needed) Calculate the integrated efficiency
4. We will also need to run ``efficiency.py config_efficiency.yml`` to get the integrated efficiency. This uses the code available in O2Physics.

### 5. Run the 2D fitter
6. Now, we have all the necessary files to run the 2D fitter. It contains a class, ``InvMassFitter2D``, and a running file, ``runInvMassFitter2D.cpp``. Again, ensure that the data's AO2D and the tree directory are correct. Give it all the previous files, and run!

And that's it :)
