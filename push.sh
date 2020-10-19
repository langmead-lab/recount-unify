#!/bin/sh

IMAGE="quay.io/broadsword/recount-unify"
VER=$(cat ver.txt)

docker push $* ${IMAGE}:${VER}
#docker push $* ${IMAGE}:latest
