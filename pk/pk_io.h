char *
pk_io_read_file(char const filename[]);


#ifdef PK_IO_IMPLEMENTATION

#include <stdio.h>
#include <stdlib.h>


char *
pk_io_read_file(char const filename[])
{
        FILE *fp = fopen(filename, "rb");
        if (!fp)
                return NULL;

        /* TODO: not portable */
        fseek(fp, 0, SEEK_END);
        long const size = ftell(fp);
        rewind(fp);

        char *string = malloc(size + 1);
        size_t const read = fread(string, 1, size, fp);
        fclose(fp);
        string[read] = '\0';

        return string;
}

#endif
