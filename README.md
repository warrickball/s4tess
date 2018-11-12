# A synthetic sample of short-cadence solar-like oscillators for TESS (S4TESS)

This repository contains the scripts used to compute the mock
catalogue of TESS targets by Ball et al. (2018, submitted).

    @ARTICLE{s4tess,
       author = {{Ball}, W.~H. and {Chaplin}, W.~J. and {Schofield}, M. and
                   {Miglio}, A. and {Bossini}, D. and {Davies}, G.~R. and
                   {Girardi}, L.},
        title = "{A synthetic sample of short-cadence solar-like oscillators for TESS}",
      journal = {ApJS, accepted,},
         year = 2018
      }

The output
files from various states of the project are available
[here](https://zenodo.org/record/1470155) on
[Zenod](https://zenodo.org).

The complete pipeline makes use of

* [MESA](http://mesa.sourceforge.net) (r7385)
* [GYRE](https://bitbucket.org/rhdtownsend/gyre/wiki/Home)
(specifically, [this fork](https://bitbucket.org/warrickball/gyre/commits/b24514643680e7b5c16a564cbe4bb16178e0610f?at=Edbeta_dx)
of v5.2)
* [AADG3](https://warrickball.github.io/AADG3) (v3.0.0a)
* [TOMSO](https://tomso.readthedocs.io) (v0.0.9)

as well as standard Python packages [NumPy](http://www.numpy.org/),
[SciPy](https://www.scipy.org/scipylib/index.html),
[Astropy](http://www.astropy.org/) and
[Matplotlib](https://matplotlib.org/).
