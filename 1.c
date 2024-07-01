#include <stdarg.h>                                                             
#include <stdio.h>         

void vout(FILE* stream, char* fmt, ...);
char fmt1[] = "%s  %s  %s\n";

int main(void)
{
    FILE* stream;

    vout(stream, fmt1, "Sat", "Sun", "Mon");
}

void vout(FILE* stream, char* fmt, ...)

{
    va_list arg_ptr;

    va_start(arg_ptr, fmt);
    vfprintf(stream, fmt, arg_ptr);
    va_end(arg_ptr);
}