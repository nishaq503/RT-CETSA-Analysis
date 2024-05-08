#!/bin/bash
version=$(<VERSION)
# container_name=polusai/rt_cetsa_analysis_test_7
container_name=polusai/rt_cetsa_analysis

inpDir=./data
outDir=./tmp
container_input_dir="/inpDir"
container_output_dir="/outDir"
# command="Rscript main.R"
command="python3 __main__.py"

# docker run -v $inpDir:/${container_input_dir} \
#             -v $outDir:/${container_output_dir} \
#             --user $(id -u):$(id -g) \
#             ${container_name}:${version} \
#             ${command} \
#             --params ${container_input_dir}/test_exp_param_full.csv \
#             --values ${container_input_dir}/test_exp_curve_all.csv \
#             --plate ${container_input_dir}/platemap.xlsx \
#             --outdir ${container_output_dir}

docker run -v $inpDir:/${container_input_dir} \
            -v $outDir:/${container_output_dir} \
            --user $(id -u):$(id -g) \
            --rm -it \
            ${container_name}:${version} \
            bash
