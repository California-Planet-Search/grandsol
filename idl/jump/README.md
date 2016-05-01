This directory contains supplementary IDL routines that illustrate periodic features in flat field images obtained with the HIRES middle (green) chip. These periodic features in the flat are probably caused by non-uniform pixel sizes, which in turn are caused by limitations of the CCD fabrication process. Larger pixels appear up to 0.6% brighter, while smaller pixels appear up to 0.6% fainter.

Use the following IDL commands to generate figures. See Issue #23 for sample output.
```
fits_read, 'j550225.fits', im, hd, exten=2
imsum = total(float(im), 1, /double)
plot_jump_flat, imsum
plot_jump_phase, imsum
```

See Issue #23 for sample plots.
