# LR_SEDfit
This repo contains a very project specific Spectral Energy Distribution (SED) fitting code optimised for refitting the SEDs of redshift ~3-4 galaxies
from the FOURstar Galaxy Evolution Survey (ZFOURGE). A key difference when compared with SED fits provided by the ZFOURGE team is that here, 
galaxy redshift is fixed rather than being included as a free parameter. Thus, LR_SEDfit is ideally used for galaxies that already have very well
constrained redshifts (i.e. from spectroscopic detection of emission lines). Indeed, this is the primary usage of this code by our team. 

Other key differences when compared to previous fits are as follows:
*SEDs are fit using Binary Population and Spectral Synthesis (BPASS) models which produce larger quantities of ionizing radiation
*No assumptions are made regarding the star-formation history of galaxies
*Default dust attenuation curve is that of Reddy et al. 2016, however other popular attenuation curves can be selected as well
*NOTE THAT EMISSION LINES ARE NOT INCLUDED IN FITS (however by default LRSED_fit ignores photometry possibly contaminated by Ly-alpha and OIII 5007)

Regarding the second point, star-formation histories provided by SED fits produced by LR_SEDfit can be considered non-parametric. The fitting method
of LR_SEDfit is __linear regression__. In a nutshell, LR_SEDfit loops through each of the possible BPASS models available (differing only in metallicity)
and for each model loops through a range of possible dust attenuation values (E(B-V)). For each model+E(B-V) combination, LR_SEDfit then determines
the linear combination of each of the 50 different aged SEDs included in the given BPASS metallicity model that minimized the cost function: 

<img src="/cost_func.png" width="700">

Where *m* is the number of photometric bands being considered,*X* is the grid of model photometric fluxes, *Y* are the observed fluxes for a
given object, and *$\sigma_Y$* are the associated measurement errors on the observed fluxes. Here the function *h($\Theta$, X)* describes a 1D vector
that represents the output SED where:

<img src="/htheta.png" width=500>

The gradients function is then defined as the partial derivatives of the cost function with respect to each element of *$\Theta$*. For each model and E(B-V) combination, the optimal $\Theta$ is determined using *scipy.optimize.minimize* with the truncated Newton algorithm. 

It should be noted that this implementation does not allow for mixing of metallicities between different aged subpopulations and also assumes that
the dust attenuation of all subpopulations is the same.

Usage information coming soon...
