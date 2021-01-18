# A synthetic sample of short-cadence solar-like oscillators for TESS (S4TESS)

This repository contains the scripts used to compute the mock
catalogue of TESS targets by
[Ball et al. (2018, ApJS, 239, 34)](https://ui.adsabs.harvard.edu/abs/2018ApJS..239...34B).

    @ARTICLE{s4tess,
           author = {{Ball}, Warrick H. and {Chaplin}, William J. and {Schofield}, Mathew
                     and {Miglio}, Andrea and {Bossini}, Diego and {Davies}, Guy R.
                     and {Girardi}, L{\'e}o},
            title = "{A Synthetic Sample of Short-cadence Solar-like Oscillators for TESS}",
          journal = {ApJS},
         keywords = {stars: oscillations: including pulsations,
                     Astrophysics - Solar and Stellar Astrophysics},
             year = 2018,
            month = dec,
           volume = {239},
           number = {2},
              eid = {34},
            pages = {34},
              doi = {10.3847/1538-4365/aaedbc},
    archivePrefix = {arXiv},
           eprint = {1809.09108},
     primaryClass = {astro-ph.SR},
           adsurl = {https://ui.adsabs.harvard.edu/abs/2018ApJS..239...34B},
          adsnote = {Provided by the SAO/NASA Astrophysics Data System}
    }

The output files from various states of the project are available
[on Zenodo](https://zenodo.org/record/1470155).

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
