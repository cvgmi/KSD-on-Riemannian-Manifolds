# KSD-on-Riemannian-Manifolds

This repository contains the code for computing the Kernel Stein Discrepancy (KSD) and Minimum Kernel Stein Discrepancy Estimator (MKSDE) of widely-used distribution families on commonly-encountered manifolds, and also contain the code for conducting the composite goodness-of-fit test. The corresponding algorithms are presented in the paper:

Theory and Applications of Kernel Stein Discrepancy on Riemannian Manifolds

## 1. Input data and Parameters

### Data
- _X_: a 3D array, where each X[,,i] represents a sample matrix

### Parameters
- _model_: a list that incorporates all other parameters




### Manifolds
- Stiefel manifold $\mathcal{V}_ r(N)$: including the sphere $\mathcal{V}_ 1(N)=\mathbb{S}^{N-1}$  and the rotational group $\mathcal{V}_ {{N-1}}(N)=\text{SO}(N)$.
- Grassmann manifold $\mathcal{G}_ r(N)$

### Distribution families
- Matrix Fisher family $\text{MF}(F)$
- Matrix Bingham family $\text{MB}(A)$
- Matrix Fisher-Bingham family $\text{MFB}(A,F)$

### Kernel
- Gaussian kernel: $\kappa(x,y)$


## 2. Functions

- _KSD(X,model)_: compute the values of empirical KSD
- _MKSDE(X,model,GoF)_: compute the MKSDE of given samples, if _GoF=TRUE_, then also output the p values of the composite goodness-of-fit test

