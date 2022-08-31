# Ye-etal_2022

**How to choose a model sensitivity analysis method as a function of computation costs and the number of parameters**

Haochen Ye<sup>1\*</sup>, Robert E. Nicholas<sup>2,3</sup>, Vivek Srikrishnan<sup>4</sup> and Klaus Keller <sup>5</sup>

<sup>1 </sup> Department of Geosciences, the Pennsylvania State University, University Park, PA, USA

<sup>2 </sup> Earth and Environmental System Institute, the Pennsylvania State University, University Park, PA, USA

<sup>3 </sup> Department of Meteorology and Atmospheric Science, the Pennsylvania State University, University Park, PA, USA
  
<sup>4 </sup> Department of Biological and Environmental Engineering, Cornell University, Ithaca, NY, USA
   
<sup>5 </sup> Thayer School of Engineering, Dartmouth College, Hanover, NH, USA
  
\* corresponding author: Haochen Ye hxy46@psu.edu

## Abstract

Complex models with many parameters and large uncertainties are common in many sectors including environmental and earth system modeling. Sensitivity analysis reveals how the variation of input model parameters and their interactions influences the model outputs. Careful sensitivity analysis helps researchers understand the robustness of the model and where to improve the model. Choosing a sound sensitivity analysis method can pose challenges. It can become computationally infeasible to evaluate the large number of model runs required by standard sensitivity analysis methods in high-dimension parameter spaces. Common approaches to tackle these challenges are to (i) emulate the expensive model using a computationally faster approximation and (ii) to sample adaptively. These approaches are, however, not a panacea, as they also involve computational costs and can introduce approximation errors. Here we compare the performances of four methods through a perfect-model experiment by sampling a wide range of model computational costs and numbers of model parameters. We find that the standard Sobol’ method is the fastest for computationally cheap models; the Kriging emulation method is the fastest for low-dimensional models with medium running time; the Adaptive Kriging combined with Monte Carlo Sampling (AKMCS) is the fastest for high-dimensional computationally expensive models; the Bayesian Adaptive Spline Surface (BASS) method is the fastest for the remaining high-dimensional models. Our results provide guidance on how to choose a sensitivity analysis method for a range of model parameters and model run-times. This guidance can inform choices of the use of computational resources, especially when the model of interest is high-dimensional or computationally expensive. 

## Journal reference

To be updated

## Reproduce my experiment

Running all the scripts will take an extremely long time because the experiments for high-dimensional emulators are slow (this can take up to several weeks). In order to avoid this long waiting time (and the heavy computational burden), we provide an option to reproduce the experiments for low-dimensional models. The scripts in this repo activate this optional low-dimensional test experiment by default. Detailed explanation of how to replicate the experiment is shown below.

Before starting, make sure the required R packages are installed under R version 4.1.0 (the version we use in the paper). These libraries are: sensobol, lhs, GPfit, BASS, plot.matrix, RColorBrewer

Step 1: Run 1_Sobol.R

Step 2: Run 2_Kriging.R, 3_BASS.R, 4_AKMCS.R. These three scripts can be run separately (if computational condition allows) as they are independent from each other.

Step 3: Run 5_traceplot_data.R

Note 1: For scripts '2_Kriging.R', '3_BASS.R', '4_AKMCS.R', be sure to make sure the parameter 'Check_switch' is 1. This only tests the experiments in 2D and 5D scenarios for these emulation-based methods. The "1_Sobol.R" is fast enough to run through all the test dimensions.

Note 2: Users should expect exactly the same outputs as in "Example_data" folder **except for** ave_eval_time files and every files start with "T_" in each subfolder. This is because these files record the computational wall time, which is different based on different machines and random noises. However their relative magnitudes should be still similar.

## Reproduce my figures

To reproduce the figures, users should use the full outputs in the 'data' folder. 

Step 1: Run 6_traceplot.R, and users should see Figure_2 in the "Figures" folder.

Step 2: Run 7_Gridplot.R, and users should see Figure_3 and Figure_4 in the "Figures" folder.


## Detailed explanation of each script: 

1_Sobol.R uses the 'sensobol' package to analyze the standard Sobol' analysis.

2_Kriging.R uses the 'GPfit' package to perform the Sobol' analysis based on the Kriging emulation method.

3_BASS.R uses the 'BASS' package to perform the Sobol' analysis based on the Bayesian Adaptive Spline Surface (BASS) emulation method.

4_AKMCS.R uses an adaptive Kriging combined with Monte Carlo Sampling (AKMCS) method based on the Kriging emulation method.

5_traceplot_data.R prepares the data of different methods' results under different random seeds in one of the test model.

6_traceplot.R makes the figure 2 in the paper (using the data generated by 5_traceplot_data.R).

7_Gridplot.R makes the figure 3 and 4 in the paper.


