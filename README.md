

Note:  This software has largely been replaced by the ODELAY-python files. 

Files for Processing ODELAY-ODELAM data on a High Performance compute cluster with the software PBSPro installed and matlab/v92.

once logged into the host excecute the following commands from the directory with ODELAM_HPC_v1.m 

    module load matlab/v92

    mcc -m ODELAM_HPC_v1.m

    export ODDIR="/Path/to/data/directory"
    export MATLIB="/Path/to/MATLAB_Runtime/v92"

next you will need to initialize the experiment with the following command:

    ./run_ODELAM_HPC_v1.sh "$MATLIB" InitExp "$ODDIR"

finally excecute the pbs run file to process the experiment in parellel.  Please make sure to change the names of paths and email addresses and such things that are required in the PBS run file and probably rename the file as something to represent the experiment.  

    qsub -q workq -J 1-80 ODELAM_Example_Experiment_Process_File.pbs

