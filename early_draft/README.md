Making an early draft of combining results from 3.4

Updates:

* updating u to 38s and a few more visits (in response to throughput changes)
* adding in the early uniform rolling
* making long gap blobs respect pre-scheduled observations.
* minor shifts to DDF central positions
* updates to DDF season. First season now auto-scales rather than being manually set.
* taking LMC/SMC out of rolling
* updated LMC/SMC filter balance. Trying to turn down u and g a little so dark time not oversubscribed


--looks like the flag to turn off uniform could use some more work, but if we want uniform releases it's fine.

XXX--may want to add some of the early science things
XXX--Do we add ToO's into the baseline?

trying to add ToO events to early_draft_too.py
