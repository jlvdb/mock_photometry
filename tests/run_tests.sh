#!/usr/bin/env bash

if [ ! -d table_tools ]
then
    git clone https://github.com/jlvdb/table_tools.git
else
    cd table_tools
    git pull
    cd ..
fi

python .mocks_extended_object_sn \
    -i mice.fits \
    -o temp.fits \
    --bulge-ratio bulge_fraction \
    --bulge-size bulge_length \
    --disk-size disk_length \
    --ba-ratio disk_axis_ratio \
    --flux-frac 0.6 \
    --psf 1.0 \
    --filters sdss_r \
    --scale 2.0 \
    --threads 6

python .merge_data.py mice.fits temp.fits result.fits

python .mocks_photometry_realisation \
    -i result.fits \
    -o temp.fits \
    --filters sdss_r_true \
    --limits 24 \
    --significance 2.0 \
    --sn-detect 3.0 \
    --sn-factors sn_factor_sdss_r \
    --seed 12345

python .merge_data.py result.fits temp.fits result.fits

rm temp.fits

python run_and_compare.py
