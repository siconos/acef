#!/bin/sh
# Static library demo

# Create static library's object file, libhello-static.o.
# I'm using the name libhello-static to clearly
# differentiate the static library from the
# dynamic library examples, but you don't need to use
# "-static" in the names of your
# object files or static libraries.

gcc -Wall  -c -o parser.o parser.c

# Create static library.

ar rcs parser.a parser.o

gcc  -c test.c
gcc  -o test test.o -ldl
