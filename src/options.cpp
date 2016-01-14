#include "options.h"

#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <stdlib.h>
#include <stddef.h>
#include <stdarg.h>

void tc_options_help()
{
    fprintf(stderr,
        "TC Optimizing Compiler\n"
        "Usage:\n"
        "\n"
        "\ttc <source file> <algorithm> <options>\n"
        "\n"
        "Algorithms:\n"
        "\n"
        "\t--stencil-tiling\t- concurrent start tiling for stencils\n"
        "\t--free-schedule-tiling\t- dynamic free scheduling for tiles\n"
        "\n"
        "Options:\n"
        "\n"
        "\t-b\t\t- block size, e.g. -b 256 -b S1:128,128\n"
        "\n"
        "e.g.: \ttc ./examples/polybench/gemm.scop.c --free-schedule-tiling -b S1:256,256 -b S2:128\n"
        "\n"
    );
}

void tc_options_error(const char* msg, ...)
{
    fprintf(stderr, "Error: ");
    
    va_list varargs;
    va_start(varargs, msg);
    vfprintf(stderr, msg, varargs);
    va_end(varargs);
    
    fprintf(stderr, "\n");
    
    exit(1);
}

char* tc_options_source_file(int argc, char* argv[])
{
    char* file = argv[1];
    
    return file;
}

enum tc_algorithm_enum tc_options_algorithm(int argc, char* argv[])
{
    const char* algorithm = argv[2];
    
    enum tc_algorithm_enum value = tc_algorithm_enum_unknown;
    
    if (0 == strcmp(algorithm, "--stencil-tiling"))
    {
        value = tc_algorithm_enum_stencil_tiling;
    }
    else if (0 == strcmp(algorithm, "--free-schedule-tiling"))
    {
        value = tc_algorithm_enum_free_schedule_tiling;
    }
    
    return value;
}

std::map<std::string, std::vector<int> > tc_options_blocks(int argc, char* argv[])
{
    std::map<std::string, std::vector<int> > blocks;
    
    blocks["__DEFAULT__"] = std::vector<int>(20, 256);
    
    for (int i = 0; i < argc; ++i)
    {
        if (0 == strcmp(argv[i], "-b"))
        {
            if (i + 1 >= argc)
            {
                tc_options_error("Missing block size");
            }
            else
            {
                char* str = argv[i + 1];           
                
                if (isdigit(str[0]))
                {
                    blocks["__DEFAULT__"] = std::vector<int>(20, atoi(str));
                }                
                else
                {
                    char* label = strtok(str, ":");                    
                    char* sizes = strtok(NULL, ":");
                    
                    if (NULL == label || NULL == sizes)
                    {
                        tc_options_error("Invalid block size");
                    }
                    else
                    {
                        char* size = strtok(sizes, ",");
                        
                        while (size != NULL)
                        {                            
                            blocks[label].push_back(atoi(size));
                            
                            size = strtok(NULL, ",");
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
    }
    
    return blocks;
}
