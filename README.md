# StickyPriceModel

This repository contains the most recent dynare and matlab files for the sticky price model. 

This folder is small, so it should be easy to work with on Board computers. Just use the download and upload links:

https://github.com/pedm/StickyPriceModel/archive/master.zip

https://github.com/pedm/StickyPriceModel/upload/master


Oct 24 Update:
Estimation.m now contains everything related to estimation. Here you can choose the parameters to estimate, set the bounds, etc. Meanwhile the mod file just runs the model with one set of parameters.


Albert: I was playing around with it and reached the results below and attached, which are fairly similar to what you got earlier, and which are quite good I think. I only did a couple of minor changes to the setting (see Estimation.m attached): I added the “sigman” parameter, which I think is necessary to get the correct impact response of R&D (is that right? Or are you scaling the model IRFs somewhere?); and I also played around with the bounds a bit. What I suggest we do until tomorrow is that we both experiment a bit more (e.g. trying different starting values, maybe playing a little with the bounds here and there). I also want to try excluding some parameters to see how far it can get. Once we are satisfied with this part I think we have more than enough to start writing the draft.
