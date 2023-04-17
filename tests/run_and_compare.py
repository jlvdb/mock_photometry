import astropandas as apd
import numpy as np

import sys, os
sys.path.insert(0, os.path.dirname(os.getcwd()))
import mock_photometry


if __name__ == "__main__":

    expect = apd.read_fits("result.fits")

    r_eff = mock_photometry.effective_radius_parallel(
        r_disk=expect["disk_length"], r_bulge=expect["bulge_length"],
        bulge_frac=expect["bulge_fraction"], flux_frac=0.6)
    assert np.allclose(r_eff, expect["R_E"])

    aper_a, aper_b = mock_photometry.aperture_size(
        1.0, r_eff, expect["disk_axis_ratio"], scale=2.0)
    assert np.allclose(aper_a, expect["aper_a_sdss_r"])
    assert np.allclose(aper_b, expect["aper_a_sdss_r"]*expect["aper_ba_ratio_sdss_r"])

    sn_weight = mock_photometry.aperture_sn_weight(1.0, aper_a, aper_b)
    assert np.allclose(sn_weight, expect["sn_factor_sdss_r"])

    magnitude, mag_error = mock_photometry.magnitude_realisation(
        expect["sdss_r_true"], 24.0, confidence=2.0, sn_weight=sn_weight,
        seed=12345, sn_detect=3.0, non_detect_mag=99)
    assert np.allclose(magnitude, expect["sdss_r_true_obs"])
    mask = mag_error != 24.0
    assert np.allclose(mag_error[mask], expect["sdss_r_true_obserr"][mask])

    print("OK")
