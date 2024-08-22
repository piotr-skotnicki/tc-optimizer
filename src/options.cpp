#include "options.h"
#include "config.h"
#include "debug.h"

#include <isl/space.h>

#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <stdlib.h>
#include <stddef.h>
#include <stdarg.h>
#include <limits.h>

#define TC_STR_(A) #A
#define TC_STR(A) TC_STR_(A)

void tc_options_help()
{
    fprintf(stderr,
        "TC Optimizing Compiler " TC_CONF_VERSION "\n"
        "\n"
        " Usage:\n"
        "\n"
        "    tc <input.c> <algorithm> <scheduling> <codegen> [<closure>] [<options>...]\n"
        "\n"
        " Algorithms:\n"
        "\n"
        "    --stencil-tiling             Concurrent start tiling for stencils\n"
        "    --regular-tiling             Tiling with regular tile shapes\n"
        "    --correction-tiling          Tiling with LT tiles correction\n"
        "    --correction-inv-tiling      Tiling with GT tiles correction\n"
        "    --merge-tiling               Tiling with tiles merging\n"
        "    --split-tiling               Tiling with tiles splitting\n"
        "    --mod-correction-tiling      Tiling with LT cyclic tiles modified correction\n"
        //"    --mod-correction-inv-tiling  Tiling with GT cyclic tiles correction\n"
        "\n"
        " Schedulers:\n"
        "\n"
        "    --lex-scheduling               Lexicographic order execution\n"
        "    --isl-scheduling               Integer set library scheduler\n"
        "    --isl-wave-scheduling          Integer set library scheduler with wavefronting\n"
        "    --feautrier-scheduling         Integer set library scheduler (Feautrier scheduling)\n"
        "    --sfs-single-scheduling        Tiling of synchronization-free slices with single sources\n"
        "    --sfs-multiple-scheduling      Tiling of synchronization-free slices with multiple sources\n"
        "    --sfs-tile-scheduling          Tile-wise synchronization-free slices\n"
        "    --free-scheduling              Free scheduling based on R^+\n"
        "    --free-rk-scheduling           Free scheduling based on R^k\n"
        "    --free-finite-scheduling       Exact free scheduling for finite graphs\n"
        "    --dynamic-free-scheduling      Dynamic free scheduling\n"
        "\n"
        " Code generators:\n"
        "\n"
        "    --serial-codegen       Serial code generator\n"
        "    --omp-for-codegen      OpenMP parallel for generator\n"
        "    --omp-task-codegen     OpenMP parallel task generator\n"
        "    --omp-gpu-codegen      OpenMP offloading to GPU target\n"
        "\n"
        " Transitive closure:\n"
        "\n"
        "    --isl-map-tc           ISL normalized map transitive closure (default)\n"
        "    --isl-union-map-tc     ISL union map transitive closure\n"
        "    --floyd-warshall-tc    Floyd-Warshall algorithm\n"
        "    --iterative-tc         Iterative algorithm\n"
        "    --tarjan-tc            Tarjan algorithm for finite graphs\n"
        "\n"
        " Options:\n"
        "\n"
        "    -b <value>           Tile size, e.g. -b 256 -b S1:128,128 (default: " TC_STR(TC_CONF_DEFAULT_TILE_SIZE) ")\n"
        "    --debug   | -d       Verbose mode\n"
        "    --report             Generate tile statistics report (use -R for each parameter)\n"
        //"    --time               Measure calculations time\n"
        "    --inline             Always inline loop bounds expressions\n"
        "    -D <name>=<value>    Define parameter value, e.g. -D M=2000 -D N=2600\n"
        "    -R <name>=<value>    Set parameter value for report generation, e.g. --report -R M=2000 -R N=2600\n"
        "    --cache <value>      Cache line length in bytes (default: " TC_STR(TC_CONF_DEFAULT_CACHE_LINE) ")\n"
        "    --use-macros         Use macro definitions in place of statements\n"
        "    --version | -v       Print compiler info\n"
        "    --help    | -h       Print help\n"
        "\n"
        " e.g.:\n"
        "    ./src/tc ./examples/stencils/heat-1d.scop.c --stencil-tiling --omp-for-codegen -b 150,25000 --debug\n"
        "    ./src/tc ./examples/polybench/bicg.scop.c --correction-tiling --sfs-single-scheduling --omp-for-codegen -b 8\n" // --time
        "    ./src/tc ./examples/polybench/trisolv.scop.c --merge-tiling --free-scheduling --omp-task-codegen -b S1:16 -b S2:16,8 -b S3:16\n"
        "\n"
    );
}

void tc_options_credits()
{
    fprintf(stderr,
        "TC Optimizing Compiler\n"
        "version " TC_CONF_VERSION "\n"
        "\n"
        "Piotr Skotnicki <pskotnicki@zut.edu.pl>\n"
        "\n"
        "West Pomeranian University of Technology\n"
        "Faculty of Computer Science and Information Technology\n"
        "ul. Zolnierska 49, 71-210 Szczecin, Poland\n"
        "\n"
    );
}

struct tc_options* tc_options_alloc(int argc, char* argv[])
{
    struct tc_options* options = (struct tc_options*)malloc(sizeof(struct tc_options));
    
    options->argc = argc;
    options->argv = argv;
    
    return options;
}

void tc_options_free(struct tc_options* options)
{
    free(options);
}

char* tc_options_get_command_line(struct tc_options* options)
{
    int argc = options->argc;    
    char** argv = options->argv;
    
    int len = 0;
        
    for (int i = 0; i < argc; ++i)
    {
        len += strlen(argv[i]) + 1;
    }
    
    len += (len % sizeof(int));
    
    char* command_line = (char*)malloc(len * sizeof(char));
    command_line[0] = '\0';
        
    for (int i = 0; i < argc; ++i)
    {
        strcat(command_line, argv[i]);
        
        if (i != argc - 1)
        {
            strcat(command_line, " ");
        }
    }
    
    return command_line;
}

int tc_options_get_int(struct tc_options* options, const char* short_name, const char* long_name)
{
    int value = 0;
    
    const char* str = tc_options_get_string(options, short_name, long_name);    
    
    if (!isdigit(str[0]))
    {
        tc_error("Missing integer value for %s.", short_name);
    }
    else
    {
        value = atoi(str);    
    }
        
    return value;
}

long tc_options_get_long(struct tc_options* options, const char* short_name, const char* long_name)
{   
    long value = 0;
        
    const char* str = tc_options_get_string(options, short_name, long_name);
    
    if (!isdigit(str[0]))
    {
        tc_error("Missing integer value for %s.", short_name);
    }
    else
    {
        value = atol(str);    
    }
    
    return value;
}

const char* tc_options_get_string(struct tc_options* options, const char* short_name, const char* long_name)
{
    int argc = options->argc;
    char** argv = options->argv;
    
    const char* value = 0;
        
    const char* name = (NULL != long_name ? long_name : (NULL != short_name ? short_name : "(unknown)"));
    
    int found = 0;
    
    for (int i = 0; i < argc; ++i)
    {        
        if ((NULL != short_name && 0 == strcmp(argv[i], short_name)) || (NULL != long_name && 0 == strcmp(argv[i], long_name)))
        {   
            found = 1;
            if (i + 1 >= argc)
            {
                tc_error("Missing value for %s.", name);
            }
            else
            {
                char* str = argv[i + 1];           
                
                value = str;
            }
            break;
        }
    }
    
    if (!found)
    {
        tc_error("Missing required option %s.", name);
    }
    
    return value;
}

int tc_options_is_set(struct tc_options* options, const char* short_name, const char* long_name)
{
    int argc = options->argc;
    char** argv = options->argv;
    
    int found = 0;
    
    for (int i = 0; i < argc; ++i)
    {        
        if ((NULL != short_name && 0 == strcmp(argv[i], short_name)) || (NULL != long_name && 0 == strcmp(argv[i], long_name)))
        {   
            found = 1;            
            break;
        }
    }
        
    return found;
}

const char* tc_options_source_file(struct tc_options* options)
{
    char* file = options->argv[1];
    
    return file;
}

const char* tc_options_output_file(struct tc_options* options)
{
    return tc_options_get_string(options, "-o", "--out");
}

int tc_options_is_verbose(struct tc_options* options)
{
    return tc_options_is_set(options, "-d", "--debug");
}

int tc_options_is_report(struct tc_options* options)
{
    return tc_options_is_set(options, NULL, "--report");
}

int tc_options_cache_line(struct tc_options* options)
{
    int cache_line = TC_CONF_DEFAULT_CACHE_LINE;
    
    if (tc_options_is_set(options, NULL, "--cache"))
    {
        cache_line = tc_options_get_int(options, NULL, "--cache");
    }
    
    return cache_line;
}

enum tc_algorithm_enum tc_options_algorithm(struct tc_options* options)
{    
    enum tc_algorithm_enum value = tc_algorithm_enum_unknown;
    
    static const char* strings[] = { "--stencil-tiling", "--regular-tiling", "--correction-tiling", "--correction-inv-tiling", "--merge-tiling", "--split-tiling", "--mod-correction-tiling" };
    static enum tc_algorithm_enum values[] = { tc_algorithm_enum_stencil_tiling, tc_algorithm_enum_regular_tiling, tc_algorithm_enum_correction_tiling, tc_algorithm_enum_correction_inv_tiling, tc_algorithm_enum_merge_tiling, tc_algorithm_enum_split_tiling, tc_algorithm_enum_mod_correction_tiling };
    
    for (int i = 0; i < sizeof(strings) / sizeof(*strings); ++i)
    {        
        if (tc_options_is_set(options, NULL, strings[i]))
        {
            if (tc_algorithm_enum_unknown == value)
            {
                value = values[i];
            }
            else
            {
                tc_error("More than one algorithms specified.");
            }
        }
    }
    
    if (tc_algorithm_enum_unknown == value)
    {
        tc_error("No algorithm specified.");
    }
    
    return value;
}

enum tc_scheduling_enum tc_options_scheduling(struct tc_options* options)
{    
    enum tc_scheduling_enum value = tc_scheduling_enum_unknown;
    
    static const char* strings[] = { "--lex-scheduling", "--isl-scheduling", "--isl-wave-scheduling", "--feautrier-scheduling", "--sfs-tile-scheduling", "--sfs-single-scheduling", "--sfs-multiple-scheduling", "--free-rk-scheduling", "--free-scheduling", "--free-finite-scheduling", "--dynamic-free-scheduling" };
    static enum tc_scheduling_enum values[] = { tc_scheduling_enum_lex, tc_scheduling_enum_isl, tc_scheduling_enum_isl_wavefronting, tc_scheduling_enum_feautrier, tc_scheduling_enum_sfs_tile, tc_scheduling_enum_sfs_single, tc_scheduling_enum_sfs_multiple, tc_scheduling_enum_free_rk, tc_scheduling_enum_free_karl, tc_scheduling_enum_free_finite, tc_scheduling_enum_free_dynamic };
    
    for (int i = 0; i < sizeof(strings) / sizeof(*strings); ++i)
    {        
        if (tc_options_is_set(options, NULL, strings[i]))
        {
            if (tc_scheduling_enum_unknown == value)
            {
                value = values[i];
            }
            else
            {
                tc_error("More than one scheduling specified.");
            }
        }
    }
    
    if (tc_scheduling_enum_unknown == value)
    {
        tc_warn("No scheduling specified.");
    }
    
    return value;
}

enum tc_codegen_enum tc_options_codegen(struct tc_options* options)
{    
    enum tc_codegen_enum value = tc_codegen_enum_unknown;
    
    static const char* strings[] = { "--serial-codegen", "--omp-for-codegen", "--omp-task-codegen", "--omp-gpu-codegen" };
    static enum tc_codegen_enum values[] = { tc_codegen_enum_serial, tc_codegen_enum_omp_cpu_for, tc_codegen_enum_omp_cpu_task, tc_codegen_enum_omp_gpu };
    
    for (int i = 0; i < sizeof(strings) / sizeof(*strings); ++i)
    {        
        if (tc_options_is_set(options, NULL, strings[i]))
        {
            if (tc_codegen_enum_unknown == value)
            {
                value = values[i];
            }
            else
            {
                tc_error("More than one code generators specified.");
            }
        }
    }
    
    if (tc_codegen_enum_unknown == value)
    {
        tc_error("No code generator specified.");
    }
    
    return value;
}

enum tc_transitive_closure_enum tc_options_transitive_closure(struct tc_options* options)
{    
    enum tc_transitive_closure_enum value = tc_transitive_closure_enum_unknown;
    
    static const char* strings[] = { "--isl-map-tc", "--isl-union-map-tc", "--floyd-warshall-tc", "--iterative-tc", "--tarjan-tc" };
    static enum tc_transitive_closure_enum values[] = { tc_transitive_closure_enum_isl_map, tc_transitive_closure_enum_isl_union_map, tc_transitive_closure_enum_floyd_warshall, tc_transitive_closure_enum_iterative, tc_transitive_closure_enum_tarjan };
    
    for (int i = 0; i < sizeof(strings) / sizeof(*strings); ++i)
    {        
        if (tc_options_is_set(options, NULL, strings[i]))
        {
            if (tc_transitive_closure_enum_unknown == value)
            {
                value = values[i];
            }
            else
            {
                tc_error("More than one transitive closure algorithms specified.");
            }
        }
    }
    
    return value;
}

std::map<std::string, std::vector<int> > tc_options_blocks(struct tc_options* options)
{
    int argc = options->argc;
    char** argv = options->argv;
    
    std::map<std::string, std::vector<int> > blocks;
    
    blocks["__DEFAULT__"] = std::vector<int>(20, TC_CONF_DEFAULT_TILE_SIZE);
    
    for (int i = 0; i < argc; ++i)
    {
        if (0 == strcmp(argv[i], "-b"))
        {
            if (i + 1 >= argc)
            {
                tc_error("Missing block size.");
            }
            else
            {
                char str[256];
                
                strcpy(str, argv[i + 1]);
                
                const char* label = NULL;
                char* sizes = NULL;
                
                if (isdigit(str[0]))
                {
                    blocks["__DEFAULT__"].clear();
                    
                    label = "__DEFAULT__";                    
                    sizes = str;
                }                
                else
                {
                    label = strtok(str, ":");
                    sizes = strtok(NULL, ":");
                }
                
                if (NULL == label || NULL == sizes)
                {
                    tc_error("Invalid block size.");
                }
                else
                {
                    char* size = strtok(sizes, ",xX-./");

                    while (size != NULL)
                    {                            
                        blocks[label].push_back(atoi(size));

                        size = strtok(NULL, ",xX-./");
                    }

                    int s = blocks[label].back();
                    for (int j = blocks[label].size(); j < 20; ++j)
                    {
                        blocks[label].push_back(s);
                    }
                }
            }
        }
    }
    
    return blocks;
}

std::vector<std::vector<std::string> > tc_options_groups(struct tc_options* options)
{
    int argc = options->argc;
    char** argv = options->argv;
    
    std::vector<std::vector<std::string> > groups;
    
    for (int i = 0; i < argc; ++i)
    {
        if (0 == strcmp(argv[i], "-g"))
        {
            if (i + 1 >= argc)
            {
                tc_error("Missing group specification.");
            }
            else
            {
                char str[256];
                
                strcpy(str, argv[i + 1]);          
                                    
                char* statement = strtok(str, ",");

                if (NULL == statement)
                {
                    tc_error("Invalid group specification.");
                }
                else
                {
                    std::vector<std::string> group;
                    
                    while (statement != NULL)
                    {                            
                        group.push_back(statement);

                        statement = strtok(NULL, ",");
                    }
                    
                    groups.push_back(group);
                }
            }
        }
    }
    
    return groups;
}

__isl_give isl_set* tc_options_get_defines(struct tc_options* options, __isl_keep isl_ctx* ctx)
{
    return tc_options_collect_values(options, "-D", NULL, ctx);
}

__isl_give isl_set* tc_options_get_report_bounds(struct tc_options* options, __isl_keep isl_ctx* ctx)
{
    return tc_options_collect_values(options, "-R", NULL, ctx);
}

__isl_give isl_set* tc_options_collect_values(struct tc_options* options, const char* short_name, const char* long_name, __isl_keep isl_ctx* ctx)
{   
    int argc = options->argc;
    char** argv = options->argv;
    
    const char* name = (NULL != long_name ? long_name : (NULL != short_name ? short_name : "(unknown)"));
    
    isl_space* space = isl_space_alloc(ctx, 0, 0, 0);
    
    isl_set* bounds = isl_set_universe(space);
    
    bounds = isl_set_params(bounds);
        
    for (int i = 0; i < argc; ++i)
    {
        if ((NULL != short_name && 0 == strcmp(argv[i], short_name)) || (NULL != long_name && 0 == strcmp(argv[i], long_name)))
        {
            if (i + 1 >= argc)
            {
                tc_error("Missing value for %s.", name);
            }
            else
            {
                char str[256];
                
                strcpy(str, argv[i + 1]);
                
                char* param = strtok(str, "=");
                char* value_str = strtok(NULL, "=");

                if (NULL == param || NULL == value_str)
                {
                    tc_error("Invalid value for %s.", name);
                }
                else
                {
                    int value = atoi(value_str);

                    bounds = isl_set_insert_dims(bounds, isl_dim_param, 0, 1);
                    
                    bounds = isl_set_set_dim_name(bounds, isl_dim_param, 0, param);
                    
                    bounds = isl_set_fix_si(bounds, isl_dim_param, 0, value);
                }
            }
        }
    }
    
    return bounds;
}

static int tc_options_editorial_distance(const char* a, const char* b)
{
    const int M = strlen(a);
    const int N = strlen(b);
    
    int** tab = (int**)malloc((M + 1) * sizeof(int*));

    for (int i = 0; i <= M; ++i)
    {
        tab[i] = (int*)malloc((N + 1) * sizeof(int));
        tab[i][0] = i;
    }

    for (int j = 1; j <= N; ++j)
    {
        tab[0][j] = j;
    }

    /* ./tc ../examples/other/levenshtein.scop.c --regular-tiling --lex-scheduling --serial-codegen -b 4 */
    #define min(x,y)    ((x) < (y) ? (x) : (y))
    #define floord(n,d) (((n)<0) ? -((-(n)+(d)-1)/(d)) : (n)/(d))
    #pragma scop
    for (int ii0 = 0; ii0 <= floord(M - 1, 4); ii0 += 1)
      for (int ii1 = 0; ii1 <= floord(N - 1, 4); ii1 += 1)
        for (int c2 = 4 * ii0 + 1; c2 <= min(M, 4 * ii0 + 4); c2 += 1)
          for (int c3 = 4 * ii1 + 1; c3 <= min(N, 4 * ii1 + 4); c3 += 1)
            tab[c2][c3] = min(min(tab[c2 - 1][c3] + 1, tab[c2][c3 - 1] + 1), tab[c2 - 1][c3 - 1] + ((a[c2 - 1] == b[c3 - 1]) ? 0 : 1));
    #pragma endscop
    
    const int retval = tab[M][N];
    
    for (int i = 0; i <= M; ++i)
    {
        free(tab[i]);
    }
    free(tab);

    return retval;
}

void tc_options_check_spelling(struct tc_options* options)
{
    static const char* strings[] = {
        "--stencil-tiling", "--regular-tiling", "--correction-tiling", "--correction-inv-tiling", "--merge-tiling", "--split-tiling", "--mod-correction-tiling",
        "--lex-scheduling", "--isl-scheduling", "--isl-wave-scheduling", "--feautrier-scheduling", "--sfs-tile-scheduling", "--sfs-single-scheduling", "--sfs-multiple-scheduling", "--free-scheduling", "--free-rk-scheduling", "--free-finite-scheduling", "--dynamic-free-scheduling",
        "--serial-codegen", "--omp-for-codegen", "--omp-task-codegen", "--omp-gpu-codegen",
        "--isl-map-tc", "--isl-union-map-tc", "--floyd-warshall-tc", "--iterative-tc", "--tarjan-tc",
        "-b", "-R", "--report", "--cache", "-d", "--debug", "-D", "--version", "-v", "--help", "-h", /*"--braces", */"--inline", /*"--time", */"--use-macros",
        "-g", "--out", "-o",
        "-m", "--max"
    };
    
    int argc = options->argc;
    char** argv = options->argv;
        
    for (int i = 1; i < argc; ++i)
    {
        if ('-' == argv[i][0])
        {
            int found = 0;
            
            for (int j = 0; j < sizeof(strings) / sizeof(*strings); ++j)
            {
                if (0 == strcmp(argv[i], strings[j]))
                {
                    found = 1;
                    break;
                }
            }            
    
            if (!found)
            {
                int best_score = INT_MAX;
                
                const char* best_string = "";
                
                for (int j = 0; j < sizeof(strings) / sizeof(*strings); ++j)
                {
                    const int score = tc_options_editorial_distance(argv[i], strings[j]);
                    
                    if (score < best_score)
                    {
                        best_score = score;
                        best_string = strings[j];
                    }
                }
                
                tc_error("Unknown option: `%s'. Did you mean `%s' ?\nIf not, type `--help' for the list of available options.", argv[i], best_string);                
            }    
        }
    }
}

