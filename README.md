# multigpu-test-code
This code mimics the basic MPI+OpenACC tasks of PSI's MAS Solar MHD code, for use with testing multi-GPU multi-node clusters
  
Example compile line using the NVIDIA HPC SDK 22.2 on A100 GPUs:  
`mpif90 psi_multigpu_test_code.f -acc=gpu -gpu=cc80 -Minfo=accel -o psi_multigpu_test_code`
  
Example OpenMPI<4 launch command (assuming 2 GPUs per socket):  
`mpiexec -npersocket 2 ./psi_multigpu_test_code`  
