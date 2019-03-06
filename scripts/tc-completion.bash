#!/bin/bash
complete -f -W "--stencil-tiling --regular-tiling --correction-tiling --correction-inv-tiling --merge-tiling --split-tiling \
                --lex-scheduling --sfs-tile-scheduling --sfs-single-scheduling --sfs-multiple-scheduling --free-scheduling --free-rk-scheduling --free-finite-scheduling --dynamic-free-scheduling \
                --serial-codegen --omp-for-codegen --omp-task-codegen \
                --isl-map-tc --isl-union-map-tc --floyd-warshall-tc --iterative-tc --tarjan-tc \
                -b -R --report --cache -d --debug -D --version -v --help -h --inline --time --use-macros \
                -g --out -o \
                -m --max" \
    tc