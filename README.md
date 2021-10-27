# Sensitivity

Sobol_Morris.R analyzes a simple two-variable reliability problem with Sobol' method and Morris method. For Sobol' method, we include random Monte Carlo (MC) sampling, Latin Hypercube Sampling (LHS) and quasi MC sampling (Sobol' sequence). The detailed explanations of plotting and sampling methods are described in the script comments.

Sobol_Morris.ipynb is a Jupyter notebook illustration of the code with outputs.

Sobol_sensobol.R analyzes the same problem with sensobol R package.

Sobol_SALib.py analyzes the same problem with SALib Python library. (this library can only use Saltelli sampling)

[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/yhaochen/Sensitivity/HEAD?labpath=Gaussian_process.ipynb)
