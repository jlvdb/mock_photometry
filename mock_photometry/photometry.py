"""
This module provides functions to compute photometry realisations based on
limiting magnitudes and galaxy aperture sizes.
"""

from future import __annotations__

import warnings

import numpy as np

from mock_photometry.effective_radius import _TFloatArray


def compute_flux_error(
    mag_limit: _TFloatArray,
    confidence: float,
    zeropoint: _TFloatArray = 0.0,
    sn_weight: _TFloatArray = None
) -> _TFloatArray:
    """
    Computes the flux error for a given magnitude limit.

    The magnitude limit is specified as magnitude value and conficence, i.e. a
    2-sigma limit of 25 mag is `mag_limit=25` and `confidence=2.0`. Additionally
    an aperture size correction can be applied.

    Parameters
    ----------
    mag_limit : float or array
        Observational magnitude limit, can be a global value or per object.
    confidence: float
        Sigma confidence of the magnitude limit.
    zeropoint : float or array (optional)
        Magnitude zeropoint of the magnitude limit.
    sn_weight : float or array (optional)
        Aperture size correction weight (see `apertures.aperture_sn_weight`).

    Returns
    -------
    flux_error : float or array
        Flux error for the given magnitude limit and aperture size.
    """
    mag_lim_corr = mag_limit - zeropoint
    flux_error = 10 ** (-0.4 * mag_lim_corr) / confidence
    if sn_weight is not None:
        flux_error = flux_error / sn_weight
    return flux_error


def draw_flux_realisation(
    model_flux: _TFloatArray,
    flux_error: _TFloatArray,
    seed: int = None
) -> _TFloatArray:
    """
    Compute a Gaussian flux realisation from a given (perfect) model flux and
    flux error.

    Parameters
    ----------
    model_flux : float or array
        Model flux.
    flux_error: float
        Flux error.
    seed : int (optional)
        Random generator seed.

    Returns
    -------
    realisation : float or array
        Model flux perturbed with Gaussian noise based on flux error.
    """
    rng = np.random.default_rng(seed)
    return rng.normal(model_flux, flux_error)


def compute_flux_realisation(
    model_flux: _TFloatArray,
    mag_limit: _TFloatArray,
    confidence: float,
    zeropoint: _TFloatArray = 0.0,
    sn_weight: _TFloatArray = None,
    seed: int = None
) -> tuple[_TFloatArray, _TFloatArray]:
    """
    Shorthand for calling `compute_flux_error` and `draw_flux_realisation`.

    Parameters
    ----------
    model_flux : float or array
        Model flux.
    mag_limit : float or array
        Observational magnitude limit, can be a global value or per object.
    confidence: float
        Sigma confidence of the magnitude limit.
    zeropoint : float or array (optional)
        Magnitude zeropoint of the magnitude limit.
    sn_weight : float or array (optional)
        Aperture size correction weight (see `apertures.aperture_sn_weight`).
    seed : int (optional)
        Random generator seed.

    Returns
    -------
    flux : float or array
        Model flux perturbed with Gaussian noise based on flux error.
    flux_error : float or array
        Flux error for the given magnitude limit and aperture size.
    """
    flux_error = compute_flux_error(mag_limit, confidence, zeropoint, sn_weight)
    flux = draw_flux_realisation(model_flux, flux_error, seed)
    return flux, flux_error


def compute_magnitude_realisation(
    model_mag: _TFloatArray,
    mag_limit: _TFloatArray,
    confidence: float,
    zeropoint: _TFloatArray = 0.0,
    sn_weight: _TFloatArray = None,
    seed: int = None,
    sn_detect: float = 1.0,
    non_detect_mag: float = np.nan
) -> tuple[_TFloatArray, _TFloatArray]:
    """
    Compute a magnitude realision for given observing conditions, see
    `compute_flux_realisation`.
    
    Additionally applies a signal-to-noise cut and assigns a placeholder value
    for objects that have fluxes below this limit, the error will be set to the
    magnitude limit.

    Parameters
    ----------
    model_mag : float or array
        Model magnitude.
    mag_limit : float or array
        Observational magnitude limit, can be a global value or per object.
    confidence: float
        Sigma confidence of the magnitude limit.
    zeropoint : float or array (optional)
        Magnitude zeropoint of the magnitude limit.
    sn_weight : float or array (optional)
        Aperture size correction weight (see `apertures.aperture_sn_weight`).
    seed : int (optional)
        Random generator seed.
    sn_detect : float (optional)
        A signal-to-noise threshold below which galaxies are assigned a
        the `non_detect_mag` placeholder value (default is 1.0).
    non_detect_mag : float
        The placeholder value for objects failing the signal-to-noise cut
        (default is NaN).

    Returns
    -------
    magnitude : float or array
        Model magnitude perturbed with non-Gaussian noise based on flux error.
    mag_error : float or array
        Magnitude error for the given magnitude limit and aperture size.
    """
    # compute flux realisation and error from magnitudes
    model_flux = 10 ** (-0.4 * (model_mag - zeropoint))
    flux, flux_error = compute_flux_realisation(
        model_flux, mag_limit, confidence, zeropoint, sn_weight, seed)
    snr = np.maximum(flux / flux_error, 0.0)  # ignore negative signal-to-noise

    # convert back to magnitudes
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        magnitude = -2.5 * np.log10(flux) + zeropoint
        mag_error = 2.5 / np.log(10.0) / snr

    # apply signal-to-noise cut and handle objects with non-positive flux
    cut_failed = (snr < sn_detect) | np.isnan(magnitude)
    magnitude[cut_failed] = non_detect_mag
    mag_error[cut_failed] = mag_limit

    return magnitude, mag_error
