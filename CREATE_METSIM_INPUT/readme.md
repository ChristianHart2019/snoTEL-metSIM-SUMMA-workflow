# Create input files for metSIM using snoTEL data

Once the data is downloaded and MetSim is installed, paste the files from the CREATE_METSIM_INPUTS folder into a folder with the snotel.csv files. These files include three broad categories of scripts: 
1) list files used for reading file lists and creating output file names 
2) file creation scripts SNOTEL_STATE.sh, SNOTEL_DOMAIN.sh, and SNOTEL_FORCING.sh are used for creating model STATE, DOMAIN, and FORCING input files for MetSim; SNOTEL_CONFIG is used for creating a configuration script for MetSim; SNOTEL_SUBMIT is used for creating plato-based submission scripts for each snotel file; and SNOTEL_RUN runs MetSim simulations for each snotel site. Read the READ.ME file within the folder for more information about each of these files.
3) submission scripts automate the running of file creation scripts on a high-performance computing cluster. These submission scripts are the only ones that will have to be submitted.

SUBMIT_RUN.sh has to be executed last. Use command “sbatch SUBMIT_*.sh”. Other than that, it doesn’t matter the order in which the submission scripts are executed. I’d recommend the procedure of:
sbatch submit_DOMAIN.sh
sbatch submit_FORCING.sh
sbatch submit_STATE.sh
sbatch submit_CONFIG.sh
sbatch submit_SUBMIT.sh

Then, you’ll have to wait for a few minutes until the several thousand files have been created.

Afterwards, use
sbatch submit_RUN.sh

After a few more minutes, that will finish too. And, the output files will be in the form of:

snotel_*_20170101-20180831.nc

Where * is the station number and the two dates are the start and end of the simulation. 
