traceformula
============
Compile with:
g++ trace.cpp -g -o trace -lpari -lmpfr -I. -O2

Run with (for example):
(echo 1000 15 4 1000) | time ./trace
(echo 7000 97 4 1000) | time ./trace
Where the parameters are # of traces to compute, level, weight, bits of precision

Output is in the newforms/ folder
