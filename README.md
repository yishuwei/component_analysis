# Component Analysis for Multi-Channel Data

Collection of a few simple, flexible, computationally efficient tools for deriving components from multi-channel data that meet user-specified temporal, spectral, spatial, etc. criteria.

Theoretical Backgrounds
======

Decomposition based on generalized eigenvectors
------

Generalized eigenvector decomposition underlies many well-known blind source separation algorithms, such as maximum/minimum noise fraction (MNF), temporal decorrelation separation (TDSEP), common spatial pattern (CSP), canonical correlation analysis (CCA), etc. References:

Parra et al. Recipes for the linear analysis of EEG. NeuroImage 28 (2005). DOI: 10.1016/j.neuroimage.2005.05.032

de Cheveigne & Parra. Joint decorrelation, a versatile tool for multichannel data analysis. NeuroImage 98 (2014). DOI: 10.1016/j.neuroimage.2014.05.068

Adaptive spatial filters based on LCMV
------

Linearly constrained minimum variance (LCMV) beamforming is a widely used method for source localization. There are many flavors of LCMV-based spatial filters.  Reference:

Sekihara & Nagarajan. Adaptive Spatial Filters for Electromagnetic Brain Imaging. Springer, 2008.

Remarks on Covariance Estimation
======

The decomposition/beamforming functions here take the signal variance-covariance matrix (or its inverse) as input. This allows flexibility, inasmuch as users may first apply advanced covariance estimation techniques (e.g. regularization or shrinkage) before invoking the decomposition/beamforming routines.

It is important to carefully think about the best way to estimate the variance-covariance matrix given one's situation, for it can crucially affect the quality of the derived components. Note that many packages implicitly apply some covariance estimation strategy without informing the users. Simply relying on the generic defaults of the packages is not recommended, as better results are often attained when specific signal and noise characteristics are taken into account. Reference:

Engemann & Gramfort. Automated model selection in covariance estimation and spatial whitening of MEG and EEG signals. NeuroImage 108 (2015). DOI: 10.1016/j.neuroimage.2014.12.040
