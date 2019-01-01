TC Optimizing Compiler 0.2.26
=============================

Manual
------

### Usage:

    tc <input.c> <algorithm> <scheduling> <codegen> [<closure>] [<options>...]

### Algorithms:

    --stencil-tiling       Concurrent start tiling for stencils
    --regular-tiling       Tiling with regular tile shapes
    --correction-tiling    Tiling with tiles correction
    --merge-tiling         Tiling with tiles merging

### Scheduling:

    --lex-scheduling               Lexicographic order execution
    --sfs-single-scheduling        Tiling of synchronization-free slices with single sources
    --sfs-multiple-scheduling      Tiling of synchronization-free slices with multiple sources
    --sfs-tile-scheduling          Tile-wise synchronization-free slices
    --free-scheduling              Free scheduling based on R^+
    --free-rk-scheduling           Free scheduling based on R^k
    --free-finite-scheduling       Exact free scheduling for finite graphs
    --dynamic-free-scheduling      Dynamic free scheduling

### Code generators:

    --serial-codegen       Serial code generator
    --omp-for-codegen      OpenMP parallel for generator
    --omp-task-codegen     OpenMP parallel task generator

### Transitive closure:

    --isl-map-tc           ISL normalized map transitive closure (default)
    --isl-union-map-tc     ISL union map transitive closure
    --floyd-warshall-tc    Floyd-Warshall algorithm
    --iterative-tc         Iterative algorithm
    --tarjan-tc            Tarjan algorithm for finite graphs

### Options:

    -b <value>           Tile size, e.g. -b 256 -b S1:128,128 (default: 32)
    --debug   | -d       Verbose mode
    --report             Generate tile statistics report (use -R for each parameter)
    --time               Measure calculations time
    --inline             Always inline loop bounds expressions
    -D <name>=<value>    Define parameter value, e.g. -D M=2000 -D N=2600
    -R <name>=<value>    Set parameter value for report generation, e.g. --report -R M=2000 -R N=2600
    --cache <value>      Cache line length in bytes (default: 64)
    --use-macros         Use macro definitions in place of statements
    --version | -v       Print compiler info
    --help    | -h       Print help

Examples
--------

    ./src/tc ./examples/stencils/heat-1d.scop.c --stencil-tiling --omp-for-codegen -b 150,25000 --debug
    ./src/tc ./examples/polybench/bicg.scop.c --correction-tiling --sfs-single-scheduling --omp-for-codegen -b 8 --time
    ./src/tc ./examples/polybench/trisolv.scop.c --merge-tiling --free-scheduling --omp-task-codegen -b S1:16 -b S2:16,8 -b S3:16