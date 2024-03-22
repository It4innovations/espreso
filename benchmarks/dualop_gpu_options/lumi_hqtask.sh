#!/bin/bash

cd "${SLURM_SUBMIT_DIR}"

machine="${1}"
tool="${2}"
ecf_file="${3}"
factor_symmetry_fake="${4}"
output_path="${5}"
uniform_clusters_domains="${6}"
element_type="${7}"
domains_x="${8}"
domains_y="${9}"
domains_z="${10}"
elements_x="${11}"
elements_y="${12}"
elements_z="${13}"
concurrency_set="${14}"
concurrency_update="${15}"
concurrency_apply="${16}"
trsm1_factor_storage="${17}"
trsm2_factor_storage="${18}"
trsm1_solve_type="${19}"
trsm2_solve_type="${20}"
trsm_rhs_sol_order="${21}"
path_if_hermitian="${22}"
f_sharing_if_hermitian="${23}"
apply_scatter_gather_where="${24}"
transpose_where="${25}"
shift 4

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
echo "factor_symmetry_fake ${factor_symmetry_fake}" >> "${infofile}"
echo "output_path ${output_path}" >> "${infofile}"
echo "uniform_clusters_domains ${uniform_clusters_domains}" >> "${infofile}"
echo "element_type ${element_type}" >> "${infofile}"
echo "domains_x ${domains_x}" >> "${infofile}"
echo "domains_y ${domains_y}" >> "${infofile}"
echo "domains_z ${domains_z}" >> "${infofile}"
echo "elements_x ${elements_x}" >> "${infofile}"
echo "elements_y ${elements_y}" >> "${infofile}"
echo "elements_z ${elements_z}" >> "${infofile}"
echo "concurrency_set ${concurrency_set}" >> "${infofile}"
echo "concurrency_update ${concurrency_update}" >> "${infofile}"
echo "concurrency_apply ${concurrency_apply}" >> "${infofile}"
echo "trsm1_factor_storage ${trsm1_factor_storage}" >> "${infofile}"
echo "trsm2_factor_storage ${trsm2_factor_storage}" >> "${infofile}"
echo "trsm1_solve_type ${trsm1_solve_type}" >> "${infofile}"
echo "trsm2_solve_type ${trsm2_solve_type}" >> "${infofile}"
echo "trsm_rhs_sol_order ${trsm_rhs_sol_order}" >> "${infofile}"
echo "path_if_hermitian ${path_if_hermitian}" >> "${infofile}"
echo "f_sharing_if_hermitian ${f_sharing_if_hermitian}" >> "${infofile}"
echo "apply_scatter_gather_where ${apply_scatter_gather_where}" >> "${infofile}"
echo "transpose_where ${transpose_where}" >> "${infofile}"

source env/lumi.clang.aocl.mpich.suitesparsenew.rocm.sh

export ESPRESO_FAKE_TEMP_FACTOR_SYMMETRY="${factor_symmetry_fake}"

# screw this ... I assume there are no spaces in any of the argument values inside $@
command="srun -n 1 ./build/espreso -c \"${ecf_file}\" $@ > \"${output_path}/stdout.txt\" 2> \"${output_path}/stderr.txt\""
timeout -v 300s bash -c "${command}" 2> "${output_path}/timeout.txt"
exitcode=$?

rm -f core*

exit ${exitcode}
