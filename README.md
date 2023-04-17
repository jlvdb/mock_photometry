# mock_photometry

The package enables computing photometry realisations for simulated galaxy
catalogs (typically from n-body simulations). It implements an analytical
photometry model that factors in observation conditions, such as point-spread
function (PSF) and magnitude limits, but does not attempt to simulate the far
more complex image detection process:

1. Compute the effective galaxy size by integrating the galaxy light profiles.
2. Compute aperture sizes based on the intrinsic galaxy and PSF sizes.
3. Compute the photometric uncertainty from a magnitude limit.
4. Compute an additional uncertainty base on the aperture size by comparing to
   the PSF size.
5. Perturb the photometry (flux or magnitude) by adding a Gaussian random flux
   component based on the photometric uncertainty.

> When using this code in published work, please cite  
> van den Busch et al. (2020), A&A 642, A200 (arXiv:2007.01846)


## Installation

This package requires a python version >= 3.9. To install, run

    pip install 'mock_photometry @ git+https://github.com/jlvdb/mock_photometry.git'

or clone and switch into repository and run

    pip install .


## Usage


Import the package with

```python
import mock_photometry as mp
```

### Galaxy size

Start by computing the galaxy size, assuming there is a bulge and a disk
component, with sizes `r_disk` and `r_bulge` (typically given as half-light
radius). The relative importance of both components is goverend by the
bulge-to-total flux fraction `bulge_frac`. Either of the two components can be
disabled by choising `bulge_frac=0` or `bulge_frac=1` accordingly.

The main free parameter is `flux_frac`, which governs up to which radii the
Sersic profiles are integrated. The default is `0.5`, which corresponds to the
half-light radius.

```python
r_eff = mp.effective_radius_vectorized(
    r_disk, r_bulge, bulge_frac, flux_frac=0.5)
```

This operation is quite slow and can be parallelised easily, for example by
using the wrapper function `mp.effective_radius_parallel` instead.

### Aperture

Next compute an aperture, which governs the size contribution to the photometric
errors. The model assumes that galaxies are elliptical, with major axis $a$ and
minor axis $b$. The ellipticity must be provided as B/A, here `ba_ratio`, which
can also be set to `1.0` or some effective value. The function approximates the
convolution of the intrinsic galaxy size with the PSF as
$\sqrt{r^2 + r_\mathsf{PSF}^2}$.

The main free parameters are the PSF size and the `scale` parameter. The PSF
should be provided in terms of its Gaussian standard deviation, in the same
units as the galaxy size (see above). The scale parameter can be used to scale
the galaxy size to adjust the signal-to-noise curve.

```python
aper_a, aper_b = mp.aperture_size(
    psf_sigma, r_eff, ba_ratio, scale=1.0)
```

### Uncertainty from aperture size

If the galaxy is resolved, its aperture size should larger than that of a point
source, i.e. the PSF. In that case, the photometric error must be scaled by an
area correction $\sqrt{\frac{\pi a b}{\pi r_\mathsf{PSF}^2}}$.

```python
sn_weight = mp.aperture_sn_weight(psf_sigma, aper_a, aper_b)
```

### Generate photometry realisation

Finally compute a noise realisation of the perfect, simulated model fluxes based
on the observational magnitude limit and the area correction from above.

The magnitude limit is specified as magnitude value and conficence, i.e. a
2-sigma limit of 25 mag is `mag_limit=25` and `confidence=2.0`. The default
zeropoint is assumed to be `0.0`.


```python
flux, flux_error = mp.flux_realisation(
    model_flux, mag_limit=25, confidence=2.0,
    zeropoint=0.0, sn_weight=sn_weight, seed=None)
```

This random operation can be seeded for reproducible results. There is also a
corresponding function to work with magnitudes, if the catalogue does not
contain fluxes.
