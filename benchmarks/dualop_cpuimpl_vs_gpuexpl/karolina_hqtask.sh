#!/bin/bash

cd "${SLURM_SUBMIT_DIR}"

machine="${1}"
tool="${2}"
ecf_file="${3}"
output_path="${4}"
uniform_clusters_domains="${5}"
element_type="${6}"
domains_x="${7}"
domains_y="${8}"
domains_z="${9}"
elements_x="${10}"
elements_y="${11}"
elements_z="${12}"
dual_operator="${13}"

shift 3

mkdir -p "${output_path}"

infofile="${output_path}/info.txt"
echo -n > "${infofile}"
echo -n "datestr " >> "${infofile}"; date +%Y%m%d_%H%M%S >> "${infofile}"
echo "machine ${machine}" >> "${infofile}"
echo "tool ${tool}" >> "${infofile}"
echo -n "node " >> "${infofile}"; hostname >> "${infofile}"
echo -n "gpu " >> "${infofile}"; nvidia-smi -L >> "${infofile}"
echo -n "cpus " >> "${infofile}"; numactl -s | grep physcpubind | cut -d' ' -f2- >> "${infofile}"
echo "ecf_file ${ecf_file}" >> "${infofile}"
echo "output_path ${output_path}" >> "${infofile}"
echo "uniform_clusters_domains ${uniform_clusters_domains}" >> "${infofile}"
echo "element_type ${element_type}" >> "${infofile}"
echo "domains_x ${domains_x}" >> "${infofile}"
echo "domains_y ${domains_y}" >> "${infofile}"
echo "domains_z ${domains_z}" >> "${infofile}"
echo "elements_x ${elements_x}" >> "${infofile}"
echo "elements_y ${elements_y}" >> "${infofile}"
echo "elements_z ${elements_z}" >> "${infofile}"
echo "dual_operator ${dual_operator}" >> "${infofile}"

cusparse_version=""
if [ "${tool}" == "cudalegacy" ]; then cusparse_version="legacy"; fi
if [ "${tool}" == "cudamodern" ]; then cusparse_version="modern"; fi
source env/karolina.gcc.openblas.openmpi.suitesparsenew.cuda.sh "${cusparse_version}"

# screw this ... I assume there are no spaces in any of the argument values inside $@
command="mpirun -n 1 --bind-to numa ./build/espreso -c \"${ecf_file}\" $@ > \"${output_path}/stdout.txt\" 2> \"${output_path}/stderr.txt\""
timeout -v 300s bash -c "${command}" 2> "${output_path}/timeout.txt"
exitcode=$?

rm -f core*

exit ${exitcode}
