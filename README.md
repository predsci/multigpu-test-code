# multigpu-test-code
  
This code mimics the basic MPI+OpenACC tasks of PSI's codes, for use with testing multi-GPU multi-node clusters.  It relies on the use of GPU-aware MPI.
  
Five versions are provided:  
  
`psi_multigpu_test_code_openacc.f90`:  
OpenACC for computation, device selection, and data movement.  
  
`psi_multigpu_test_code_openacc_nodata.f90`:  
OpenACC for computation and device selection only.  This should be used with managed or unified memory.  
  
`psi_multigpu_test_code_stdpar.f90`:  
Standard Fortran (no directives at all).  This should be used with managed or unified memory and a special launch script for device selection (an example script for OpenMPI is provided).  
  
`psi_multigpu_test_code_stdpar_accdata.f90`:  
Standard Fortran with OpenACC for device selection and data movement.  
  
`psi_multigpu_test_code_stdpar_ompdata.f90`:  
Standard Fortran with OpenMP Target for data movement and OpenACC+OpenMP for device selection.  
  
  
Example compile lines using the NVIDIA HPC SDK 24.11 on GPUs:  
`mpif90 psi_multigpu_test_code_stdpar.f90         -stdpar=gpu -gpu=ccnative,mem:unified  -Minfo=accel -o psi_multigpu_test_code_stdpar`  
`mpif90 psi_multigpu_test_code_stdpar_accdata.f90 -stdpar=gpu -acc=gpu -gpu=ccnative,mem:separate  -Minfo=accel -o psi_multigpu_test_code_stdpar_accdata`  
`mpif90 psi_multigpu_test_code_stdpar_ompdata.f90 -stdpar=gpu -acc=gpu -mp=gpu -gpu=ccnative,mem:separate  -Minfo=accel -o psi_multigpu_test_code_stdpar_ompdata`  
`mpif90 psi_multigpu_test_code_openacc.f90        -acc=gpu -gpu=ccnative,mem:separate -Minfo=accel -o psi_multigpu_test_code_openacc`  
`mpif90 psi_multigpu_test_code_openacc_nodata.f90 -acc=gpu -gpu=ccnative,mem:unified  -Minfo=accel -o psi_multigpu_test_code_openacc_nodata`  
  
  
Example launch commands:  
`mpiexec -npernode 4 ./psi_multigpu_test_code_openacc`  
`mpiexec -npernode 4 ./launch_multigpu_openmpi.sh psi_multigpu_test_code_stdpar`  
  



