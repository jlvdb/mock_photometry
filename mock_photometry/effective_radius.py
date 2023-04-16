"""
This module provides functions to compute the effective radius of galaxies given
a set properties. These are the radius of the bulge and disk component
(interpreted as the effective, i.e. half-light radius in arbitrary units) and
the bulge fraction (the ratio of bulge-to-total flux/luminosity).

The bulge component is assumed to be described by a Sersic profile with index
n=4 and the disk component by a Sersic profile with index n=1.
"""

from __future__ import annotations

import multiprocessing
from typing import TypeAlias, TypeVar

import numpy as np
import pandas as pd
from numpy.typing import NDArray
from scipy.optimize import root_scalar
from scipy.special import gammainc

_TFloatArray: TypeAlias = float | NDArray[np.float_]
_TRadius = TypeVar("_TRadius", _TFloatArray)


def sersic_b_n(index: int) -> float:
    """
    Computes the constant b_n for a Sersic profile of given index.
    """
    return 1.9992 * index - 0.3271


def flux_fraction_sersic(
    radius: _TRadius,
    r_eff: _TRadius,
    index: int
) -> _TRadius:
    """
    Computes the fraction of the total flux/luminosity emitted from within a
    given radius for a Sersic profile.

    Parameters
    ----------
    radius : float or array
        Radius from within which the flux is integrated.
    r_eff : float or array
        Effective radius of the Sersic profile.
    index : int
        Index of the Sersic profile.

    Returns
    -------
    flux_frac : float or array
        Fraction of flux emitted from within radius in the same units as inputs.
    """
    r_scaled = radius / r_eff
    return gammainc(2 * index, sersic_b_n(index) * r_scaled**(1/index))


def flux_fraction_within_radius(
    radius: _TRadius,
    r_disk: _TRadius,
    r_bulge: _TRadius,
    bulge_frac: _TFloatArray
) -> _TRadius:
    """
    Computes the fraction of the total flux/luminosity emitted from within a
    given radius for a galaxy with bulge and disk component.

    Parameters
    ----------
    radius : float or array
        Radius from within which the flux is integrated.
    r_disk : float or array
        Effective radius of the n=1 Sersic disk component, zero to disable.
    r_bulge : float or array
        Effective radius of the n=4 Sersice bulge component, zero to disable.
    bulge_frac : float or array
        Ratio of bulge-to-total flux, must be within [0, 1].

    Returns
    -------
    flux_frac : float or array
        Fraction of flux emitted from within radius in the same units as inputs.
    """
    # compute fraction of total flux emmitted by disk from within radius
    if r_disk == 0.0 or bulge_frac == 1.0:
        disk_flux_frac = 0.0
    else:
        disk_frac = (1.0 - bulge_frac)
        disk_flux_frac = disk_frac * flux_fraction_sersic(radius, r_disk, 1)

    # compute fraction of total flux emmitted by bulge from within radius
    if r_bulge == 0.0 or bulge_frac == 0.0:
        bulge_flux_frac = 0.0
    else:  # evaluate the integrated Sersic n=4 profile
        bulge_flux_frac = bulge_frac * flux_fraction_sersic(radius, r_bulge, 4)

    return disk_flux_frac + bulge_flux_frac


def effective_radius(
    r_disk: float,
    r_bulge: float,
    bulge_frac: float,
    flux_frac: float = 0.5
) -> float:
    """
    Compute the radius within which a certain fraction of the total flux is
    emitted.
    
    With the default value, this corresponds to the effective radius, the radius
    from within which half of the flux is emitted. This function essentially
    computes the radius R as root of
        `flux_fraction_sersic(R, r_disk, r_bulge, bulge_frac) - flux_frac`
    using Brent's method.

    Parameters
    ----------
    r_disk : float or array
        Effective radius of the n=1 Sersic disk component, zero to disable.
    r_bulge : float or array
        Effective radius of the n=4 Sersice bulge component, zero to disable.
    bulge_frac : float or array
        Ratio of bulge-to-total flux, must be within [0, 1].
    flux_frac : float
        The fraction of the total flux for which the corresponding radius is
        computed, must be within (0, 1).

    Returns
    -------
    radius : float
        The radius within which the fraction of the total flux is emitted.
    """
    # guess a bracketing radii for the root finder
    test_radii = max(r_disk, r_bulge) * np.logspace(-5, 5, 11)
    test_frac = flux_fraction_sersic(test_radii, r_disk, r_bulge, bulge_frac)
    idx_min = np.argmin(np.abs(test_frac - flux_frac))
    r_bracket = (test_radii[idx_min-1], test_radii[idx_min+1])

    def root_function(radius: float):
        return flux_fraction_sersic(
            radius, r_disk, r_bulge, bulge_frac) - flux_frac

    solution = root_scalar(
        root_function, method="brentq", bracket=r_bracket, maxiter=100)
    return solution.root


effective_radius_vectorized = np.vectorize(effective_radius)


def effective_radius_parallel(
    r_disk: NDArray[np.float_],
    r_bulge: NDArray[np.float_],
    bulge_frac: NDArray[np.float_],
    flux_frac: float = 0.5,
    threads: int = None
) -> NDArray[np.float_]:
    """
    A wrapper for `effective_radius` that distributes the data from the input
    arrays on parallel threads using multiprocessing and checks if inputs are
    valid.
    """
    r_disk = np.asarray(r_disk)
    r_bulge = np.asarray(r_bulge)
    bulge_frac = np.asarray(bulge_frac)
    if threads is None:
        threads = multiprocessing.cpu_count()

    if np.any(r_disk < 0.0):
        raise ValueError("r_disk must not be negative")
    if np.any(r_bulge < 0.0):
        raise ValueError("r_bulge must not be negative")
    if flux_frac <= 0.0 or flux_frac >= 1.0:
        raise ValueError("flux_frac must be within (0, 1)")

    # split in to chunks and process in parallel
    with multiprocessing.Pool(threads) as pool:
        args = (
            np.array_split(r_disk, threads),
            np.array_split(r_bulge, threads),
            np.array_split(bulge_frac, threads),
            [flux_frac] * threads)
        r_eff = np.concatenate(pool.starmap(effective_radius_vectorized, args))
    return r_eff