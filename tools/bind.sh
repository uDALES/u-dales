#!/bin/bash

export LOCAL_RANK=${OMPI_COMM_WORLD_LOCAL_RANK}
export CUDA_VISIBLE_DEVICES=${LOCAL_RANK}

echo "[LOG] local rank $LOCAL_RANK: bind to $CUDA_VISIBLE_DEVICES"
echo ""

$* 

