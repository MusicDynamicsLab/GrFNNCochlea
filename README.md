# GrFNNCochlea
<img src="https://MusicDynamicsLab.github.io/Figures/tuningCurves.jpg" width="800">

GrFNN Cochlea is a cochlear model that uses Gradient Frequency Neural Networks (GrFNN) and requires the GrFNN Toolbox (https://github.com/MusicDynamicsLab/GrFNNToolbox).

This repository contains code to use a canonical cochlear model. This model is based on a coupled oscillator system representing the coupling between the basilar membrane and the organ of Corti of the cochlea. The basilar membrane (BM) is represented here as a linear oscillator, and the organ of Corti (OC) is represented as a highly nonlinear critical Hopf oscillator. Each coupled oscillator system represents a specific place along the cochlear spiral, and the natural frequencies of the oscillators in the system correspond to the center frequency of that place in the cochlea. The scripts in this repository set up a gradient-frequency network of many BM oscillators, each coupled to its respective OC oscillator in another network. This system as a whole represents a large portion of the human cochlea. 

The script `cochleaUnidirectional` sets up a two-layer, unidirectionally-coupled GrFNN parameterized to model the human cochlea. The script `cochleaBidirectional` sets up a two-layer, bidirectionally-coupled GrFNN. They both have default stimuli designed to replicate the results of the tuning curve fits in `TCfitGUI`, but the user can of course make their own stimuli to pass to the cochlea model. 

After a cochlea model has finished running, the script `cochleaMovie` may be run. This script will display a time-domain representation of the basilar membrane-organ of Corti dynamics over the course of the simulation just run, and write an .avi movie file into the current Matlab directory of that animation.

`TCfitGUI` runs a graphical interface to compare model tuning curves to auditory nerve tuning curves obtained from a macaque monkey. The user can adjust a variety of parameters to observe the effects on the tuning curves, for both unidirectional and bidirectional two-layer cochlear models as well as a single-layer model. 
