# Sensitivity
This is the repository for the paper: How to choose a model sensitivity analysis method as a function of computation costs and the number of parameters (Ye et al.)

Instructions:

To replicate the figures, download the 'data' folder into your working directory and run 6_Gridplot.R and 7_traceplot.R

To check or replicate the results in data, run the other scripts based on the order of the numbers (from 1 to 5), and you should get the same outputs as in folder 'Example_data'. For scripts '2_Kriging.R', '3_BASS.R', '4_AKMCS.R', be sure to set the parameter 'Check_switch' to 1, otherwise the running time will be extremely long!

#####
Detailed explanation of each script: 

1_Sobol.R uses the 'sensobol' package to analyze the standard Sobol' analysis.

2_Kriging.R uses the 'GPfit' package to perform the Sobol' analysis based on the Kriging emulation method.

3_BASS.R uses the 'BASS' package to perform the Sobol' analysis based on the Bayesian Adaptive Spline Surface (BASS) emulation method.

4_AKMCS.R uses an adaptive Kriging combined with Monte Carlo Sampling (AKMCS) method based on the Kriging emulation method.

5_traceplot_data.R prepares the data of different methods' results under different random seeds in one of the test model.

6_Gridplot.R makes the figure 3 and 4 in the paper.

7_traceplot.R makes the figure 2 in the paper (using the data generated by 5_traceplot_data.R.

[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/yhaochen/Sensitivity/HEAD?labpath=Gaussian_process.ipynb)
