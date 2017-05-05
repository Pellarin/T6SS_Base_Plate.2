
#!/bin/zsh
# ---
#$ -N emd_5739.map_256
#$ -cwd
#$ -S /bin/zsh
#$ -t 1-2
#$ -tc 10
#$ -q all.q
# ---


module ()
{
	eval `/usr/bin/modulecmd bash $*`
}

source ~/bin/riccardo_modules.sh

if ((SGE_TASK_ID==1))
then
/usr/bin/time -v /Bis/home/shanot/imp_sam-fast//setup_environment.sh python /Bis/home/shanot/recursive_gmm/job.py emd_5739.map 0.04 256
fi
