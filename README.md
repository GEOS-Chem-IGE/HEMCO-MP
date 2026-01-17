HEMCO-MP
========

This repository (https://github.com/GEOS-Chem-IGE/HEMCO-MP) contains the source code for
HEMCO-MP, a microplastics-enabled version of the GEOS-Chem Harmonized Emissions Component
(HEMCO). It is forked from the [official HEMCO repository](https://github.com/geoschem/HEMCO).

HEMCO is a software component for computing (atmospheric) emissions from different
sources, regions, and species on a user-defined grid. It can combine, overlay, and update
a set of data inventories ('base emissions') and scale factors, as specified by the user
through the HEMCO configuration file. Emissions that depend on environmental variables and
non-linear parameterizations are calculated in separate HEMCO extensions. HEMCO can be run
in standalone mode or coupled to an atmospheric model. A more detailed description of
HEMCO is given in Keller et al. (2014) and Lin et al (2021).

HEMCO has been coupled to several atmospheric and Earth System Models, and can be coupled
with or without using the Earth System Modeling Framework (ESMF). A detailed description
of HEMCO coupled with other models is given in Lin et al. (2021).

The initial version of the microplastics module code was written by Yanxu Zhang and is
described in Fu et al. (2023).

> [!NOTE]
> This repository is based on [version 3.6.2 of the official HEMCO code](https://github.com/geoschem/HEMCO/tree/3.6.2).


Documentation
-------------

### Online user's manual

Installation and usage instructions for HEMCO are posted online at [hemco.readthedocs.io](https://hemco.readthedocs.io/en/3.6.2/).

### References

Fu, Y., Pang, Q., Ga, S. L. Z., Wu, P., Wang, Y., Mao, M., Yuan, Z., Xu, X. Liu, K., Wang, X., Li, D., & Zhang, Y. (2023). Modeling atmospheric microplastic cycle by GEOS-Chem: An optimized estimation by a global dataset suggests likely 50 times lower ocean emissions. *One Earth*, *6*(6), 705-714. https://doi.org/10.1016/j.oneear.2023.05.012

Keller, C. A., Long, M. S., Yantosca, R. M., Da Silva, A. M., Pawson, S., & Jacob, D. J. (2014). HEMCO v1. 0: A versatile, ESMF-compliant component for calculating emissions in atmospheric models. *Geoscientific Model Development*, *7*(4), 1409-1417. https://doi.org/10.5194/gmd-7-1409-2014

Lin, H., Jacob, D. J., Lundgren, E. W., Sulprizio, M. P., Keller, C. A., Fritz, T. M., Eastham, S. D., Emmons, L. K., Campbell, P. C., Baker, B., Saylor, R. D., and Montuoro, R. (2021). Harmonized Emissions Component (HEMCO) 3.0 as a versatile emissions component for atmospheric models: application in the GEOS-Chem, NASA GEOS, WRF-GC, CESM2, NOAA GEFS-Aerosol, and NOAA UFS models. *Geoscientific Model Development*, *14*(9), 5487-5506. https://doi.org/10.5194/gmd-14-5487-2021


Support
-------

Please use [the Github issue tracker attached to this repository](https://github.com/GEOS-Chem-IGE/HEMCO-MP/issues) to report bugs or technical issues related to the plastics module code in HEMCO-MP. Issues with other modules or the broader HEMCO code should be reported on the [official HEMCO issue tracker](https://github.com/geoschem/HEMCO/issues/new/choose).


License
-------

HEMCO is distributed under the MIT license. Please see the license documents LICENSE.txt and AUTHORS.txt in the root folder.
