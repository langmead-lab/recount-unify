#!/bin/bash
set -x

IMAGE="quay.io/broadsword/recount-unify"
VER=$(cat ver.txt)

docker build $* \
    --tag ${IMAGE}:${VER} \
    --tag ${IMAGE}:latest \
    .
