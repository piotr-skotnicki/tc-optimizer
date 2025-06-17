#!/bin/bash
complete -f -W "--stencil-tiling --regular-tiling --correction-tiling --correction-inv-tiling --diamond-tiling --semi-diamond-tiling --merge-tiling --split-tiling --mod-correction-tiling \
                --lex-scheduling --isl-scheduling --isl-wave-scheduling --feautrier-scheduling --sfs-tile-scheduling --sfs-single-scheduling --sfs-multiple-scheduling --free-scheduling --free-rk-scheduling --free-finite-scheduling --dynamic-free-scheduling \
                --serial-codegen --omp-for-codegen --omp-task-codegen --omp-gpu-codegen \
                --isl-map-tc --isl-union-map-tc --floyd-warshall-tc --iterative-tc --omega-map-tc --omega-union-map-tc --tarjan-tc \
                -b -R --report --cache -d --debug -D --version -v --help -h --inline --time --use-macros \
                -y --yes \
                -g --out -o \
                -m --max \
                --drop-bounds" \
    tc
