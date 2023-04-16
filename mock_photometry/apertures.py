"""
This module provides functions to compute an analytic aperture from given
galaxy effective ratio, ellpiticities and observational PSF size. The required
convolution of the galaxy image with the PSF is approximated by assuming
both have Gaussian light profiles, which simplifies the convolution to a
addition in quadrature of galaxy and PSF "sigmas".
"""

from future import __annotations__

import numpy as np

from mock_photometry.effective_radius import _TFloatArray


def aperture_size(
    psf_sigma: float,
    r_eff: _TFloatArray,
    ba_ratio: _TFloatArray = 1.0,
    scale: float = 1.0
) -> tuple[_TFloatArray, _TFloatArray]:
    """
    Computes an analytic galaxy-image aperture for given observing conditions.

    The aperture is computed from the effective radii of the galaxies which is
    convolved with the PSF, assuming that both PSF and galaxy light-profiles are
    Gaussian.

    Parameters
    ----------
    psf_sigma : float
        PSF size give as sigma of a Gaussian.
    r_eff : float or array
        Effective size of the galaxies (typically half-light radius), given in
        same units as psf_sigma.
    ba_ratio : float or array (optional)
        Galaxy ellipticity given as minor-to-major axis ratio.
    scale : float (optional)
        Scaling factor applied to effective galaxy size before PSF convolution.

    Returns
    -------
    aper_a : float or array
        Major axis of aperture, given in same units as psf_sigma.
    aper_b : float or array
        Minor axis of aperture, given in same units as psf_sigma.
    """
    intr_a = r_eff * scale
    intr_b = intr_a * ba_ratio
    aper_a = np.sqrt(intr_a**2 + psf_sigma**2)
    aper_b = np.sqrt(intr_b**2 + psf_sigma**2)
    return aper_a, aper_b


def aperture_sn_weight(
    psf_sigma: float,
    aper_a: _TFloatArray,
    aper_b: _TFloatArray,
) -> _TFloatArray:
    """
    Computes an aperture-size correction factor for the photometric signal-to-
    noise by comparing to the PSF size.

    Parameters
    ----------
    psf_sigma : float
        PSF size give as sigma of a Gaussian.
    aper_a : float or array
        Major axis of aperture, given in same units as psf_sigma.
    aper_b : float or array
        Minor axis of aperture, given in same units as psf_sigma.

    Returns
    -------
    sn_weight
    """
    return np.sqrt(psf_sigma**2 / (aper_a * aper_b))
