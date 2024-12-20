# KSD-on-Riemannian-Manifolds

The file _KSD-on-Riemannian-Manifold.R_ contains the code for computing the Kernel Stein Discrepancy (KSD) and Minimum Kernel Stein Discrepancy Estimator (MKSDE) of widely-used distribution families on commonly-encountered manifolds, and also contain the code for conducting the composite goodness-of-fit test. The corresponding algorithms are presented in the paper:

Theory and Applications of Kernel Stein Discrepancy on Riemannian Manifolds

## 1. Arguments

All functions in the repository requires two arguments

### _X_
a 3D array that contains all the samples, where each X[,,i] is a sample matrix; 

### _model_
a list that incorporates all other parameters, including
- _model$space_: character string, specifying which manifold is the sample space
> - "Stiefel": the Stiefel manifold $\mathcal{V}_ r(N)$, including the sphere $\mathcal{V}_ 1(N)=\mathbb{S}^{N-1}$ and the special rotation group $\mathcal{V}_ {N-1}(N)=\text{SO}(N)$
> - "Grassmann": the Grassmann manifold $\mathcal{G}_ r(N)$

- _model$family_: character string, specifying which distribution family is being used
> - "MF": the matrix Fisher family $p(X)\propto\exp(F^T X)$
> - "MB": the matrix Bingham family $p(X)\propto\exp(X^T A X)$
> - "MFB": the matrix Fisher-Bingham family $p(X)\propto \exp(X^T A X + F^T X) $
  
- _model$kernel_: character string, specifying which kernel is being used
> - "Gaussian": the Gaussian kernel $\kappa(x,y)=\exp(-\frac{\tau}{2}\Vert x-y\Vert^2)$
> - "InverseQ": the inverse quadratic kernel $\kappa(x,y)= (\beta+\Vert x-y\Vert^2)^{-\gamma}$

- _model$F_: matrix, the parameter $F$ of the MF (or MFB) family, with default zero matrix
- _model$A_: matrix, the parameter $A$ of the MB (or MFB) family, with default zero matrix
- _model$tau_: scalar, the parameter $\tau$ of the Gaussian kernel, with default $1$
- _model$beta_: scalar, the parameter $\beta$ of the inverse quadratic kernel, with default $1$
- _model$gamma_: scalar, the parameter $\gamma$ of then inverse quadratic kernel, with default $1$

- _model$nprime_: integrer, the number $n'$ of generations in the goodness-of-fit test, with default $10000$

### Dimensions of the data _X_
Once _X_ is inputed, its dimensions will be automatically recorded in _model_ as following variables:
> - _model$N1_: integer, the col number of each sample matrix
> - _model$N2_: integer, the row number of each sample matrix
> - _model$r_: integer, the rank of each sample matrix
> - _model$n_: integer, number of samples

## 2. Functions

### _KSD(X,model)_
Compute the values of empirical KSD. The output is a list, including the values of $U^w_ p$, $V^w_ p$ and the matrix $(\kappa^w_ p(x_ i,x_ j))_{i j}$.

### _MKSDE(X,model,GoF)_
compute the MKSDE of given samples, if _GoF=TRUE_, then also output the p values of the composite goodness-of-fit test

### Other functions
