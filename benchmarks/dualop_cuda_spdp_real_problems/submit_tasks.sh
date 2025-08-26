#!/bin/bash

basedir="benchmarks/dualop_cuda_spdp_real_problems"

if [ ! -d "${basedir}" ]
then
    echo "must be run from espreso root directory"
    exit 7
fi

rundir="${1}"
tasksdir="${rundir}/tasks"

if [ ! -d "${tasksdir}" ]
then
    echo "rundir or taskdir '${tasksdir}' does not exist"
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
hq_outdir="${basedir}/hq_outerr/${datestr}"
mkdir -p "${hq_outdir}"

slurm_outdir="${basedir}/slurm_outerr"
mkdir -p "${slurm_outdir}"



taskdirsfile="${rundir}/taskdirs.txt"
find "${tasksdir}/" -mindepth 1 -maxdepth 1 -type d | sort > "${taskdirsfile}"



${hqbin} submit \
    --time-request=11m \
    --each-line "${taskdirsfile}" \
    --cpus="16 compact" \
    --stdout="${hq_outdir}/hqtask-%{JOB_ID}-%{TASK_ID}.o.txt" \
    --stderr="${hq_outdir}/hqtask-%{JOB_ID}-%{TASK_ID}.e.txt" \
    -- \
    "${basedir}/hqtask.sh"



rundir_name="$(echo ${rundir} | rev | cut -d/ -f1 | rev)"

echo
echo "run slurm job using:"
echo "    sbatch --array=1-10 -o ${slurm_outdir}/slurm-%j.out -e ${slurm_outdir}/slurm-%j.err ${basedir}/slurmjob.sh"
echo
echo "observe with:"
echo "    (cd ${rundir}; for d in tasks/*; do if [ -f \$d/.finished ]; then exitcode=\$(cat \$d/.finished); if [ ! \$exitcode == 0 ]; then echo "\$d"; echo exitcode \$exitcode; cat \$d/stderr.txt; echo; echo; echo; fi; fi; done > errs.txt)"
echo "    (cd ${rundir}; for d in tasks/*; do if [ -f \$d/.finished ]; then echo "\$d"; fi; done > finished.txt)"
echo
echo "later summarize with"
echo "    ./benchmarks/dualop_cuda_spdp_real_problems/summarize.py ${rundir_name}"
echo
