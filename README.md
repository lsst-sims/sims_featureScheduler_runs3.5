
[![DOI](https://zenodo.org/badge/808327321.svg)](https://zenodo.org/doi/10.5281/zenodo.13260262)


# sims_featureScheduler_runs3.5
The further adventues of simulating Rubin observing strategies


Completed simulations can be found at: https://s3df.slac.stanford.edu/data/rubin/sim-data/sims_featureScheduler_runs3.5/

Latest analysis of runs often served at: http://astro-lsst-01.astro.washington.edu:8080/


List from fed on slack:

* baseline - with unifrom rolling, 29.x sec exposures (no-snaps) (38 in uband), ToO, new galactic footprint, Roman field at some point, what filter balance on Clouds and SCP?
* noToO - same as baseline but noToO
* 4cycleRolling - same as baseline but 4 cycles of rolling
* sameFilterBalanceonMC - same as baseline but different filter balance on Clouds and SP: if it was bluer in baseline then implement same across all low dust footrpint, if it was the same implement bluer
* snaps - same as baseline byt with 2x14.x snaps in grizy
* EarlyScienceHighestDT - as baseline but high down time (pessimistic)**
* EarlyScienceHighDT - as baseline but downtime slightly higher than current expectation**
* EarlyScienceExpectedDT -  as baseline but downtime as current expectation**
* EarlyScienceExpectedLowDT -  as baseline but slightly less downtime as current expectation**
* EarlyScienceExpectedLowestDT -  as baseline but downtime very low (this may be what we havesimulated so far… :sweat_smile:)**
* SlewTime: as as baseline but slower slew time acceleration jerk (to be discussed with RL and obs strategy team)

