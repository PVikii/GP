# Galaxy photometry v2  

## Description 

GPv2 is a wrapper built around the "photutils.isophote" package, which is the Python equivalent of the IRAF ellipse task. This IRAF task is used for computing SBPs of galaxies and was implemented in the method used by Vaduvescu et al. (2006) and Lian et al. (2015). These methods are based on the iterated execution of the IRAF ellipse task with adjusted input parameters.  The GPv2 provides an  automatic way to apply that widely used technique and to obtain the physical parameters of the galaxies.  

GPv2 attempts  fitting  the galaxy profiles in 4 step. In the first step, the script generates an estimation of the fitting ellipses allowing variable centers, ellipticity and position angle (PA). Based on these results, more accurate parameters can be defined. In the second step the ellipses are re-fitted with fixed central coordinates. In the last two steps the calculations for step two are repeated while attempting to add additional constraints, fixing the ellipticity and later the PA. The  program provides  verification plots  after each step and  returns the final isophotal profile plot in magnitude.  For more details see GP_v2.pdf.

## Getting started

The GPV2 was tested in  Python 3.10.9 environment with the pacakge versions:  matplotlib 3.7.1, astropy 5.2.2, scipy 1.8.1, photutils 1.6.0, numpy 1.22.3 and argparse 1.1 . 

For installing a specific verion of the package use ```sudo apt-get intall python3-astropy==5.0```

## Usage

The GPv2 has several input parameters options from which the  required imput parameters are 'fits_name', 'target_name', 'zp': 

e.g.: ```galaxy_photometry_v2.1.py -zp 26.97 -fits_name  nostars.fits -target_name "LEDA 2308331"  ```

But most of the cases the png name is different than the name of the fits and  some of the initial parameters are chaged  for a better approximation for the galaxy  shape or size. So one can see a more complex example :

```galaxy_photometry_v2.1.py -zp 26.97 -fits_name  nostars.fits  -png_name coadd_cropped_2MASXJ0856_2.0 -target_name "LEDA 2308331" -bkg_scale 4 -sma 5 -x0 238 -y0 237.8   -chk_infl  -no_fix_ellip  -ellip 0.44  -pa 2.5 ```

The full list of options can be seen from: ```galaxy_photometry_v2.py -h ```

## License

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the  GNU General Public License for more details.

You should have received a copy of the GNU General Public License  along with this program.  If not, see <https://www.gnu.org/licenses/>.

Scientific Use: If you use this program for scientific purposes or include results derived from it in a publication, you must cite the following paper:
Pinter V., Vaduvescu, O. Deep K-band surface brightness photometry of dE galaxies, Astronomy & Astrophysics, Volume 675, id. A90, 38 pp., 07/2023, doi: 10.1051/0004-6361/202244858
