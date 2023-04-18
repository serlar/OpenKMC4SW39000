# OpenKMC4SW39000

**OpenKMC** is a high-performance parallel AKMC program developed based on [Spparks](https://spparks.github.io/).

This repository is a redesigned distribution of OpenKMC for the new generation of Sunway high-performance platform. OpenKMC is highly optimized for the SW39000 processor of the next-generation Sunway system, which can effectively utilize the heterogeneous many-core computing power and can obtain 37x acceleration.

## Enviroment
This repository needs to be compiled and run on the new-generation Sunway high-performance platform. Please make sure the repository is installed in that environment.

### Compiler
The compiler version, which ensure that OpenKMC can be compile successfully, is:

**gcc version 7.1.0 20170502 (swgcc-1307)**
## Quick Start
```bash
cd src
make swg++
cd ../SW_test/01-18MPI
./sub.sh
```
The software runs in a similar way to [Spparks](https://spparks.github.io/), see the Spparks documentation for details.

After compiling in the src directory, you will get an openkmc_swg++ binary.
Submit the job with the following command:
```bash
bsub -b -J job_name -o outfile -m 1 -q queue_name -n Number_of_MPI_processes -cgsp 64 -share_size 15000 -host_stack 256 -cross_size 32 -ldm_share_mode 4 -ldm_share_size 4 ../../src/openkmc_swg++ -in input
```
Users can set parameters on demand, such as **job_name**, **outfile**, **queue_name**, **Number_of_MPI_processes** and **input**.
The other flags(include "-cgsp 64 -share_size 15000 -host_stack 256 -cross_size 32 -ldm_share_mode 4 -ldm_share_size 4") are not recommended to be modified because there is a high probability that the program will not run after modification.

This repository contains a **\perf_log** directory. It contains a number of performance test files. Each test folder contains:
1. **input**. The input file describes the system information that OpenKMC needs to calculate, and more information is available in the [Spparks](https://spparks.github.io/) manual.
2. **sub.sh**. Job submission scripts for use on the next-generation Shineway high-performance platform.
3. **output**. The result of this test.

The **\perf_log\SW_opt_test** directory contains two tests. Each test contains two run logs **output** and **output_SSL**. These two logs are the results after and before the SW39000 adaptation optimization, respectively. From them, we can see that OpenKMC gets 37x acceleration from the many-core.

