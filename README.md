# GrFNNCochlea
<img src="https://MusicDynamicsLab.github.io/Figures/tuningCurves.jpg" width="800">

Cochlear model using Gradient Frequency Neural Network (GrFNN) toolbox.

The script `cochleaUnidirectional` sets up a two-layer, unidirectionally-coupled GrFNN parameterized to model the human cochlea. The script `cochleaBidirectional` sets up a two-layer, bidirectionally-coupled GrFNN. They both have default stimuli designed to replicate the results of the tuning curve fits in `TCfitGUI`, but the user can of course make their own stimuli to pass to the cochlea model. 

After a cochlea model has finished running, the script `cochleaMovie` may be run. This script will display a time-domain representation of the basilar membrane-organ of Corti dynamics over the course of the simulation just run, and write an .avi movie file into the current Matlab directory of that animation.

`TCfitGUI` runs a graphical interface to compare model tuning curves to auditory nerve tuning curves obtained from a macaque monkey. The user can adjust a variety of parameters to observe the effects on the tuning curves, for both unidirectional and bidirectional two-layer cochlear models as well as a single-layer model. 
