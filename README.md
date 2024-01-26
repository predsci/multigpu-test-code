# multigpu-test-code
  
This code mimics the basic MPI+OpenACC tasks of PSI's MAS Solar MHD code, for use with testing multi-GPU multi-node clusters.  It relies on the use of CUDA-aware MPI.
  
Three versions are provided:   
  
`psi_multigpu_test_code.f`:  
OpenACC for computation, device selection, and data movement.  
  
`psi_multigpu_test_code_nodata.f`:  
OpenACC for computation and device selection only.  This should be used with managed or unified memory.  
  
`psi_multigpu_test_code_stdpar.f`:  
Standard Fortran version.  No OpenACC at all.  This should be used with managed or unified memory and a special launch script for device selection (an example script for OpenMPI is provided).  
  
Example compile lines using the NVIDIA HPC SDK 23.11 on A100 GPUs:  
`mpif90 psi_multigpu_test_code.f -acc=gpu -gpu=cc80,nomanaged,nounified -Minfo=accel -o psi_multigpu_test_code`  
`mpif90 psi_multigpu_test_code_nodata.f -acc=gpu -gpu=cc80 -Minfo=accel -o psi_multigpu_test_code_nodata`  
`mpif90 psi_multigpu_test_code_stdpar.f -stdpar=gpu -gpu=cc80 -Minfo=accel -o psi_multigpu_test_code_stdpar`  

Example launch commands:  
`mpiexec -npernode 4 ./psi_multigpu_test_code`  
`mpiexec -npernode 4 ./psi_multigpu_test_code_nodata`  
`mpiexec -npernode 4 ./launch_psi_multigpu_test_code_stdpar.sh psi_multigpu_test_code_stdpar`  
  
  



