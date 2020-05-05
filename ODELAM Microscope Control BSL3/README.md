ODELAM Microscope Control

This MATLAB software is utilized for collecting data on a Nikon Ti-E microscope.  
This software is given as an example only and is not supported in any way.


## Installing

This software requires:
MATLAB                                                Version 9.4         (R2018a)
Global Optimization Toolbox                           Version 3.4.4       (R2018a)
Image Processing Toolbox                              Version 10.2        (R2018a)
MATLAB Compiler                                       Version 6.6         (R2018a)
Optimization Toolbox                                  Version 8.1         (R2018a)
Parallel Computing Toolbox                            Version 6.12        (R2018a)
Statistics and Machine Learning Toolbox               Version 11.3        (R2018a) 

Also to interface with the microscope, micromanager must be installed and configured.  

For information on installing Micro-manager and configuring matlab to interface with the Micro-manager API please see:

https://micro-manager.org/wiki/Micro-Manager
https://micro-manager.org/wiki/Micro-Manager_Programming_Guide
https://micro-manager.org/wiki/Matlab_Configuration

Install micromanager on your system and use the wizard to set up hardware with your instrument.
Example micromanager config files utilized for a Ti-E are provided though they may need modifications to work on your system.  Please refer to the mircomanager website for instructions on the config files. 

In the file ActivateMicroscope.m change line 8 such that the path links to the Micro-manager config file of your choice.  
eg. mmc.loadSystemConfiguration('');

Important:  Remove any objective lens that may come into contact with the stage before initalizing the software.  This way unintentional stage movement won't accidently destroy a expensive lens.

Please check your micoscope stage and make sure that it reports position in micrometers.

Important:  set the orgin for the microscope in the file InitializeMicroscopeProperties.m line 53 mP.XYZOrigin(1,:) = [29500,-15500,8500.0]; Check that the direction of movement is not reversed.  This may require setting mP.xgrid and mP.ygrid as positve or negative to get the correct movemnt.  

Also you will need to check the camera rotation.  Make sure that the image is not reversed and if needed make changes to the image snap commands.

## Running ODELAY microscope control

To run start MATLAB and run the file ODELAY_Microscope_Control.m

You may have errors due to spcific equipment not being available. 

This software is poorly written so please use as an example but not for anything your life depends on.  

A version that utilizes python and python packages is currenly in development (as of May 2020).

## Running an experiment

Required is a microscope growth chamber, and tooling to spot cultures.  

Previous steps will require printing a slide gasket, casting agarose media, and then spoting prepared cultures onto the agarose media. 

Please see https://www.jove.com/video/55879/odelay-large-scale-method-for-multi-parameter-quantification-yeast

The software will create a set of folders and *.mat files that the images are recorded in with meta data.    






 

