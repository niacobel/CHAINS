#!/bin/bash

#########################################################################################################
###      This script will load the different modules needed by CHAINS for the different clusters      ###
#########################################################################################################

# Global requirements:
#   - Python    3.5 and later
#	- Jinja2	2 and later		(in case of required upgrade, load python 3 module then issue in command line: python -m pip install --user -U Jinja2)
#	- PyYaml	5.1 and later	(in case of required upgrade, load python 3 module then issue in command line: python -m pip install --user -U pyyaml)

# Specifically for control_launcher and results_treatment
#   - Numpy     1.14 and later  (in case of required upgrade, load python 3 module then issue in command line: python -m pip install --user -U numpy)

module --force purge

if   [ "${CLUSTER_NAME}" = "dragon1" ]; then
    module load releases/2019b
    module load Python/3.7.4-GCCcore-8.3.0
	
elif [ "${CLUSTER_NAME}" = "dragon2" ]; then
    module load releases/2019a
    module load PyYAML/5.1-GCCcore-8.2.0

elif [ "${CLUSTER_NAME}" = "lemaitre3" ]; then
    module load releases/2018b
    module load Python/3.6.6-foss-2018b

elif [ "${CLUSTER_NAME}" = "hercules2" ]; then
    module load releases/2020b
    module load Python/3.8.6-GCCcore-10.2.0

elif [ "${CLUSTER_NAME}" = "nic5" ]; then
    module load releases/2019b
    module load Python/3.7.4-GCCcore-8.3.0

elif [ "${CLUSTER_NAME}" = "lyra" ]; then
    module load releases/2020b
    module load Python/3.8.6-GCCcore-10.2.0

else
    echo "ERROR: Unknown cluster. Corresponding modules can't be loaded. Please add the cluster module information to load_modules.sh."

fi
