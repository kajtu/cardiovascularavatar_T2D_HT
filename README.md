# Cardiovascular Avatar T2D HT
These are the scripts to reproduce the results in the publication *"Hemodynamic effects of hypertension and type 2 diabetes - insights through a 4D flow MRI-based personalized cardiovascular model"* (DOI: [10.1113/JP284652](https://doi.org/10.1113/JP284652), Tunedal et al 2023).

If you use this implementation in your academic projects, please cite this paper.

The cardiovascular model file with all model equations is located in Modelfiles/avatar_HEALTH_syms_fast.m.

# Requirements
The model is created in the [AMICI toolbox](https://doi.org/10.1093/bioinformatics/btab227) in MATLAB. To compile the model, MATLAB 2017b or earlier is needed, but to run the already compiled model any later matlab verison works. The provided compiled model is compiled on Windows, but will not work on Linux or macOS. To re-compile: run GenerateModels in the folder Modelfiles. To compile the model, you need a valid C-compiler (such as xcode on Mac or MinGW on Windows. Run mex -setup to check if you have an installed compiler in matlab) and the MATLAB Symbolic Math Toolbox.

To run the other scripts, the following MATLAB toolboxes are needed: Statistics and Machine Learning Toolbox, Optimization Toolbox, Parallel Computing Toolbox, Symbolic Math Toolbox, Signal Processing Toolbox. The code was created with Matlab R2021a. Earlier matlab versions might not be compatible with many of the scripts.

Additionally, the following toolboxes are included in the folder Requirements and are needed to reproduce the results:
* For simulation, the [AMICI toolbox](https://doi.org/10.1093/bioinformatics/btab227) is needed. 
* For optimization with eSS, the [MEIGO toolbox](https://doi.org/10.1186/1471-2105-15-136) is needed.
* For optimization with MCMC, the [PESTO toolbox](https://doi.org/10.1093/bioinformatics/btx676) is needed.
* For spider plots, the [spider_plot](https://github.com/NewGuy012/spider_plot/releases/tag/19.4) function is needed.
* For the color palette magma, [MatPlotLib Perceptually Uniform Colormaps](https://www.mathworks.com/matlabcentral/fileexchange/62729-matplotlib-perceptually-uniform-colormaps) is needed.
* For Shapiro-Wilk test, the function [swtest](https://www.mathworks.com/matlabcentral/fileexchange/13964-shapiro-wilk-and-shapiro-francia-normality-tests) is needed (in the folder Statistical_tests).

# Re-create result figures
Re-create the main results: Run the script *createFigures_hypertensionT2D*. This will take a while, since it takes time to save the figures in high resolution. 

Re-create the subject-specific predictions: run *createFigures_hypertensionT2D_predictions*. This will take several minutes, since all parameter sets will be simulated for each subject. It also takes time to save the figures in high resolution.

Figure 6 - correlations - was created with the script plot_parameterCorrelations in the folder Plotting. However, the data needed to re-create the plot is not publically available due to ethical restrictions.
Similarly, Table 1 - cohort characteristics - was created with the script cohortTable in the folder Plotting, but the data is not publically available due to ethical restrictions.

# Optimization
All scripts for parameter estimation are found in the folder Optimzation.

To re-do the optimization to find the best fit to data, go to the Optimizaion folder and run runParamEstimationSeveral_ESS_asNSC. 
The optimization ran in 8 batches of 10 subjects at a time at the Swedish NSC. The script runParamEstimationSeveral_ESS_asNSC runs the optimization in the same way as it was done at the NSC, but you need to change the variable patrange manually. On a normal workstation, this could take days or weeks.

To re-run the MCMC sampling, run runParamEstimationSeveral_MCMC_asNSC. The MCMC sampling ran in 3 batches at the swedish NSC: subject 1-32, subject 33-64, and subject 65-80. The script runParamEstimationSeveral_ESS_asNSC runs the optimization in the same way as it was done at the NSC - but you need to manually change the variables patrange and numberOfSubjects. On a normal workstation, this could take days or weeks.

The parameter estimation for different training data for SBP and DBP to calculate of model sensitivity to that data was done using the scripts *runHEALTH_changeBP_SBPmax.sh* etc, which runs *EstimateParametersESS_HEALTH_changeBP.m* several times for different data input of SBP and DBP for subject 5 (the HT+T2D subject with best fit to data) and saves the results in the folder Parameters\BPsensitivity.

# Simulation
An example of using the model to do a simple simulation can be found in the script *simpleModelSimulation.m*.
The model file is found in the folder Modelfiles (*avatar_HEALTH_syms_fast.m*).

## Author
Kajsa Tunedal (kajsa.tunedal@liu.se)

## License
The MIT License (MIT)

Copyright (c) 2023 Kajsa Tunedal

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

