# meds_inject
Inject fake signals into MEDS files

examples
========

```bash
seed=8212
cp DES0003-3832_r_meds-y3v02.fits DES0003-3832_r_meds-y3v02-g1-0.02.fits
meds-inject --file DES0003-3832_r_meds-y3v02-g1-0.02.fits --g1 -0.02 --g2 0.0 --seed ${seed}

cp DES0003-3832_r_meds-y3v02.fits DES0003-3832_r_meds-y3v02-g1+0.02.fits
meds-inject --file DES0003-3832_r_meds-y3v02-g1-0.02.fits --g1 0.02 --g2 0.0 --seed ${seed}
```

Requirements
============
numpy
galsim
fitsio
tqdm
