# ppopp17_artifact

The artifact includes all of the programs that are needed to reproduce the evaluation results in our paper "Combining SIMD and Many/Multi-core Parallelism for Finite State Machines with Enumerative Speculation". The programs should be executed on an Intel Xeon Phi SE10P coprocessor with 61 cores. An Intel ICC compiler is required for compiling. Python is also require for executing the scripts that generate some of the inputs. The artifact will print out the execution time of different versions of the evaluated applications on different datasets. The execution time should be very close to those reported in the paper. The artifact will also print out the speculation success rates of different versions of applications, and the rates should be exactly the same as reported in the paper. 
