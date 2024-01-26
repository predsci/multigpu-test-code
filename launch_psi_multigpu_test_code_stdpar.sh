#!/bin/bash
# This script assumes you are using OpenMPi and are launching with 1 GPU per MPI rank.

export CUDA_VISIBLE_DEVICES="$OMPI_COMM_WORLD_LOCAL_RANK"
echo $CUDA_VISIBLE_DEVICES

# Execute command
exec $*
