# multigpu-test-code
This code mimics the basic MPI+OpenACC tasks of PSI's MAS Solar MHD code, for use with testing multi-GPU multi-node clusters
  
Example compile line:  
`mpif90 psi_multigpu_test_code.f -acc=gpu -gpu=cc75,cuda11.6 -Minfo=accel -o psi_multigpu_test_code`
  
Example launch command:  
`mpiexec -npersocket 2 ./psi_multigpu_test_code`  
