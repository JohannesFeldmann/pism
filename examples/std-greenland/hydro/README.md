std-greenland/hydro/
===========

These subglacial hydrology runs are documented in the paper

*  E. Bueler, W. Van Pelt (2014) _Mass-conserving subglacial hydrology in the_
   _Parallel Ice Sheet Model_. Geoscientific Model Development Discussions
   7 (4) pp. 4705–4775, doi:10.5194/gmdd-7-4705-2014

These runs start with input data from files generated one level up (directory
`examples/std-greenland/`), especially from script `gridseq.sh`.  Thus these
runs are based on the well-documented runs in Chapter 1 of the PISM User's
Manual, specifically those in the "grid sequencing" section.

Do

    $ ln -s ../pism_Greenland_5km_v1.1.nc
    $ ln -s ../Greenland_5km_v1.1.nc

Run a grid sequencing like in Chapter 1 of the User's Manual to get
`g2km_gridseq.nc` (or similar).  Then do

    $ ./run-decoupled.sh 5 g2km_gridseq.nc              # 5 year runs

This run produces six files: `routing-decoupled.nc`, `ex_routing-decoupled.nc`, `ts_routing-decoupled.nc`, `distributed-decoupled.nc`, `ex_distributed-decoupled.nc`, `ts_distributed-decoupled.nc`.

To generate map-plane figures from the paper do:

    $ ln -s ../basemapfigs.py
    $ ./allfigs.sh g2km_gridseq

The second script `allfigs.sh` calls other figure-generating
scripts: `basemapfigs.py`, `genGreenfig.sh`, `genscatfig.sh`,
and `showPvsW.py`.  It also uses `mogrify` from the [Imagemagick](http://www.imagemagick.org/) tools.)

To show some additional hydrology time-series do

    $ ./hydro-tsshow.py ts-routing.png ts_routing-decoupled.nc
    $ ./hydro-tsshow.py ts-distributed.png ts_distributed-decoupled.nc

