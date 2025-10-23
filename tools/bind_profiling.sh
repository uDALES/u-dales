#!/bin/bash

export nn=${OMPI_COMM_WORLD_SIZE}
export LOCAL_RANK=${OMPI_COMM_WORLD_LOCAL_RANK}
export CUDA_VISIBLE_DEVICES=${LOCAL_RANK}

echo "[LOG] local rank $LOCAL_RANK: bind to $CUDA_VISIBLE_DEVICES"
echo ""

if [[ ${LOCAL_RANK} == 0 ]]; then
  echo "Sizes $nn"
  nome=ReportUDALES_ngpu${nn}
  echo $nome
  #nsys profile --trace=cuda,nvtx,mpi,openacc --cuda-memory-usage=true --gpu-metrics-device=0 --output=${nome} $*
  nsys profile --trace=cuda,nvtx,mpi,openacc --cuda-memory-usage=true --output=${nome} $*
  #nsys profile --trace=cuda,nvtx,mpi,openacc --cuda-memory-usage=true $*
else
  $*
fi
