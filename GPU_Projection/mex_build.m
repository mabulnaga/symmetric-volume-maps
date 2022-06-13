% Make sure to modify to point to your CUDA installation.
mexcuda -O -g tet_proj.cu -L/usr/local/cuda/lib64 -lcublas -lcusolver;
