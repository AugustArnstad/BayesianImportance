# BayesianImportance
**Bayesian Variable Importance for GLMMs using INLA**

`BayesianImportance` is an R package designed to compute Bayesian variable importance metrics for Generalized Linear Mixed Models (GLMMs) utilizing the Integrated Nested Laplace Approximation (INLA) methodology.

## Features
- **Bayesian Variable Importance Computation**: Allows for the quantification of the importance of predictors in GLMMs in a Bayesian framework.
- **INLA Integration**: Leverages the computational advantages of INLA, a popular method for Bayesian inference for latent Gaussian models.
- **Support for Various GLMMs**: Compatible with a wide range of generalized linear mixed models.
- **Extensible**: Designed with the modern R user in mind, offering a range of utilities to further expand upon the base functionality.

## Installation
To install the latest version of `BayesianImportance` from GitHub, use the following command:
```R
# If not already installed, install the 'devtools' package
if(!require(devtools)) install.packages("devtools")

# Install BayesianImportance
devtools::install_github("AugustArnstad/BayesianImportance")
``` 

## Usage
To compute the Bayesian variable importance for your GLMMs, follow the basic structure:

```R
library(BayesianImportance)
# ... [Your GLMM setup and data]
# result <- bayesian_importance(your_model)
```

Detailed examples and tutorials will be made available in future releases.

## Documentation
Further documentation and function references can be found within the package. Use the standard R help and documentation commands to access detailed information about each function.

## Contributing
Contributions to `BayesianImportance` are welcome. Please ensure that you adhere to standard coding practices and include tests for any new features. Open a pull request with details about your changes.

## License
Just me

## Acknowledgements
INLA team for their outstanding work on the INLA methodology and R package.
My counsellor Stefanie Muff, Associate Professor, Department of Mathematical Sciences, NTNU Trondheim

