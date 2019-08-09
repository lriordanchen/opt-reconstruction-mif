Here are instructions for how to do tomographic reconstructions with your data collected from one
of the Mesoscopic Imaging Facility OPT set-ups.

This folder contains a Jupyter notebook for running tomographic reconstructions step-by-step.
It also contains a python script for batch running several tomographic reconstuctions.  

There is also a notebook (Data_organizer_wetlabOPT) which takes the data in the form it came out from the OPT 
(tiffs with different names in different folders) at MIF and puts it into the correct formats for reconstruction.  
Both reconstruction codes require the data to be in hdf5 files, but tiff stacks 
are also generated for users who prefer handling data in this format.

The python script can be run by entering the following style of command into your (Anaconda) command prompt 
from the directory containing the script opt-reconstruction.py:
python opt-reconstruction.py PathToSample1 PathToSample2 PathToSample3 ...

For example, if you have two octopus samples you want to reconstruct, your command might look like
C:\Users\me\Random_codes\tomopy-lec\> python opt-reconstruction.py F:\Octopus\experiment1 F:\Octopus\experiment2

By default, it will not run tilt correction or mlem.  If you want to run either of these, you can add 
tilt_correction=true and/or mlem=true (all lower-case) to the end of the command, example:
python opt-reconstruction.py F:\Octopus\experiment1 F:\Octopus\experiment2 mlem=true