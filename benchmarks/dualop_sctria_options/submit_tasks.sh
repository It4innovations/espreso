#!/bin/bash

basedir="benchmarks/dualop_sctria_options"

if [ ! -d "${basedir}" ]
then
    echo "must be run from espreso root directory"
    exit 7
fi

rundir="${1}"

tasksdir="${rundir}/tasks"

if [ ! -d "${tasksdir}" ]
then
    echo "rundir or tasksdir '${tasksdir}' does not exist"
    exit 1
fi

hqbin="dependencies/HyperQueue/hq"

if [ ! -f "${hqbin}" ]
then
    echo "HyperQueue is not installed in dependencies"
    # just wget and extract it from github
    # or symlink to actual install dir
    exit 2
fi

if ! ${hqbin} server info > /dev/null 2> /dev/null
then
    echo "HyperQueue server is not running"
    echo "Start the server on this login node using 'hq server start' inside tmux"
    exit 3
fi

datestr="$(date +%Y%m%d_%H%M%S)"
hq_outdir="${rundir}/hq_outerr/${datestr}"
mkdir -p "${hq_outdir}"

slurm_outdir="${rundir}/slurm_outerr"
mkdir -p "${slurm_outdir}"

machine=$(cat "${rundir}/machine.txt")

num_cores_for_job="1"
if [ "${machine}" == "karolina" ]; then num_cores_for_job="16"; fi
if [ "${machine}" == "mn5" ]; then num_cores_for_job="20"; fi



taskdirsfile="${rundir}/taskdirs.txt"
find "${tasksdir}/" -mindepth 1 -maxdepth 1 -type d | sort > "${taskdirsfile}"



${hqbin} submit \
    --time-request=6m \
    --each-line "${taskdirsfile}" \
    --cpus="${num_cores_for_job} compact" \
    --stdout="${hq_outdir}/hqtask-%{JOB_ID}-%{TASK_ID}.o.txt" \
    --stderr="${hq_outdir}/hqtask-%{JOB_ID}-%{TASK_ID}.e.txt" \
    -- \
    "${basedir}/hqtask.sh"



echo
echo "run slurm job using:"
if [ "${machine}" == "karolina" ] || [ "${machine}" == "mn5" ]; then
    echo "    sbatch --array=1-10 -o ${slurm_outdir}/slurm-%j.out -e ${slurm_outdir}/slurm-%j.err ${basedir}/slurmjob_${machine}.sh"
fi
echo
echo "observe with:"
echo "    (cd ${rundir}; for d in tasks/*; do if [ -f \$d/.finished ]; then exitcode=\$(cat \$d/.finished); if [ ! \$exitcode == 0 ]; then echo "\$d"; echo exitcode \$exitcode; cat \$d/stderr.txt; echo; echo; echo; fi; fi; done > errs.txt)"
echo "    (cd ${rundir}; for d in tasks/*; do if [ -f \$d/.finished ]; then echo "\$d"; fi; done > finished.txt)"
