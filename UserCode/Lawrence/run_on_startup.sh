####
## John Gresl
## Last Edit: 5/10/2018
## Since I couldn't find an appropriate .bash file that was automatically run on startup,
## run this to create common aliases/execute other tasks.
## To get these changes to presist, use:
## $ source run_on_startup.sh
####



## Aliases
echo Creating aliases...
#alias python3="/coupp/app/home/coupp/anaconda3/bin/python3"
#alias ipython3="/coupp/app/home/coupp/anaconda3/bin/ipython3"
alias python="python"
#alias ipython2="/coupp/app/home/coupp/anaconda2/bin/ipython"
#alias SBC="cd /nashome/j/jgresl/Projects/SBCcode/"
#alias recon="cd /pnfs/coupp/persistent/grid_output/SBC-17/output"
#alias raw="cd /bluearc/storage/SBC-17-data/"


## Environment Variable Stuff
echo Setting up python Event Variables...
unset PYTHONHOME
#export PYTHONPATH=${PYTHONPATH}:/nashome/j/jgresl/Projects/
#export PATH=${PATH}:/nashome/j/jgresl/Projects/
export PYTHONPATH=${PYTHONPATH}:/Users/LawrenceLuo/Documents/College/Dark\ Matter\ Research/
export PATH=${PATH}:/Users/LawrenceLuo/Documents/College/Dark\ Matter\ Research/


## JobSub Setup
#echo Setting up jobsub suite...
#. /grid/fermiapp/products/common/etc/setups.sh
#setup jobsub_client
