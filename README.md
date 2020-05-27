#Affinity scheduling program

To run this program, please first ensure that your system has installed intel-compilers-18.

Then, please compile this program:
```bash
icc -o test -fopenmp loops.c -O3 -lm
```

Then, please input the number of threads by using this command:
```bash
export OMP_NUM_THREADS=4
```

Finally, you can execute this program and get results:
```bash
./test
```

