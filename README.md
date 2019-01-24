# Overcomplete Independent Component Analysis via SDP

This project contains implmentation of the novel overcomplete independent component analysis algorithm (OverICA) introduced in our recent paper:

Anastasia Podosinnikova, Amelia Perry, Alexander Wein, Francis Bach, Alexandre d'Aspremont, David Sontag. Overcomplete Independent Component Analysis via SDP. In Proceedings of the 22nd International Conference on Artificial Intelligence and Statistics (AISTATS), 2019.

Please cite this paper if you use this code for your research.

The algorithm is developed for the overcomplete setting where the latent dimension exceeds the observed dimention, although it also works in the complete and undercomplte settings. This algorithm estimates the overcomplete mixing matrix; standard techniques can be used for the estimation of the latent representations. In short, the algorith consists of two parts: (a) estimation of a basis of the subspace where the (modified) atoms live (obtained from generalized covariances or the fourth-order cumulant) and (b) estimation of the atoms via SDP and subsequent deflation. This project also contains the scripts necessary for reproduction of all the experiments from the paper and includes an implementation of the Fourier PCA algorithm.

## Dependencies:
* To reproduce the phase transition plot (Fig. 1), please install the [cvx optimization toolbox](http://cvxr.com/cvx/) first.
* To reproduce the synthetic experiments (Fig. 2, 3, 6, 7), please [contact](http://homes.esat.kuleuven.be/~delathau/) Lieven De Lathauwer to obtain the FOOBI code.
* To reproduce the real data experiment (Fig. 3), please download the [CIFAR-10](https://www.cs.toronto.edu/~kriz/cifar.html) dataset.

## Quick Start:
* Make sure your Matlab recognizes a C++ complier: ```mex -setup```.
* Add all dependent pahts: ```install.m```.
* If necessary to recomplie the .mex files: ```install(1)```.
* Reproduce experiments: ```reproduce_*```.
* The results of the experiments will be saved to the ```expres/``` directory.

To run the OverICA algorithm, execute ```ds = overica(X, k);``` where ```X``` is the ```p```-times-```n``` data matrix with obervations in columns and ```k``` is the desired latent dimension. The output is the ```p```-times-```k``` mixing matrix. Note that ```k``` should not exceed ```p^2/4``` for our guarantees to hold. Additionally, one can specify the details of implementation with ```ds = overica(X, k, opts);``` where ```opts``` is the structure with the following possible (optional) fields:
* ```opts.sub``` can take either ```'quad'``` or ```'gencov'``` (default) values that define whether the quadricovariance or generalized covariance implementations will be used for the first step of subspace estimation;
* ```opts.s``` defines the number ```(s*k)``` of generalized covariances to be constructed (default ```s=5```);
* ```opts.t``` sets the paramter ```t``` for generalized covariance (we found that ```t=0.1``` often works well in practice; the default is ```t = 0.05/sqrt( max( max( abs( cov(X') ) ) ) )```);
* ```opts.sdp``` can take either ```'clust'``` or ```'ada'``` or ```'semiada'``` (default) values that define whehter the clustering-based or adaptive or semi-adaptive deflation will be used.

## Questions and Feedback:
Please don't hesitate to contact me with questions regarding this code, usage of this algorithm, or bug reports.
