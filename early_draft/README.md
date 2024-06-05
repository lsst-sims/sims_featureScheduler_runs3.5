Making an early draft of combining results from 3.4

Things we are updating:

* updating u to 38s and a few more visits (in response to throughput changes)
* adding in the early uniform rolling
* taking LMC/SMC out of rolling
* updated LMC/SMC filter balance. Trying to turn down u and g a little so dark time not oversubscribed
* making long gap blobs respect pre-scheduled observations.

--looks like the flag to turn off uniform could use some more work, but if we want uniform releases it's fine.