# Sensitivity

Sobol_Morris.R analyzes a simple two-variable reliability problem with Sobol' method and Morris method. For Sobol' method, we include random Monte Carlo (MC) sampling, Latin Hypercube Sampling (LHS) and quasi MC sampling (Sobol' sequence). The detailed explanations of plotting and sampling methods are described in the script comments.

Sobol_Morris.ipynb is a Jupyter notebook illustration of the code with outputs.

Sobol_sensobol.R analyzes the same problem with sensobol R package.

Sobol_SALib.py analyzes the same problem with SALib Python library. (this library can only use Saltelli sampling)

Gaussian_process.R is an illustration example of Gaussian process model fitting.

5D_Sobol.R, 5D_Kriging.R, 5D_AKMCS.R run a 5D test problem by these three methods. Then plot_5D organizes these results and make plots.

Similarly, 42D_Sobol.R, 42D_Kriging.R run a 42D test problem.

[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/yhaochen/Sensitivity/HEAD?labpath=Gaussian_process.ipynb)
