#include "input_output.h"
#include "options.h"
#include "debug.h"

#include <stdio.h>
#include <stddef.h>

bool tc_io_confirm(struct tc_options* options, const char* question)
{
    if (tc_options_is_set(options, "-y", "--yes"))
    {
        printf("%s [y/N]: y\n", question);
        fflush(stdout);
        return true;
    }
    else
    {
        char c = '\0';

        do
        {
            printf("%s [y/N]: ", question);
            fflush(stdout);
            scanf("%c", &c);
        }
        while (c != 'Y' && c != 'y' && c != 'N' && c != 'n' && c != '\n' && c != '\r');

        return (c == 'Y' || c == 'y');
    }
}

FILE* tc_io_open_output(const char* path)
{
    FILE* file = fopen(path, "w");
    if (file == NULL)
    {
        tc_error("Failed to open output file `%s'", path);
        tc_die(tc_exit_code_input);
    }
    return file;
}

void tc_io_close(FILE* file)
{
    if (file != NULL)
    {
        fclose(file);
    }
}
