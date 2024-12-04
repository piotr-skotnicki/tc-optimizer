#include "input_output.h"
#include "options.h"

#include <stdio.h>

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
