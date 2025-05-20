#!/bin/bash
cmd="$*"
podman run -it --rm \
  --user $(id -u):$(id -g) \
  --userns=keep-id \
  -v "$PWD":"$PWD" \
  $(for d in /tmp /home/submit /work/submit /ceph/submit /scratch/submit /cvmfs /etc/grid-security; do echo -v $d:$d; done) \
  -w "$PWD" \
  gitlab-registry.cern.ch/cms-cloud/combine-standalone:v9.2.1 "$cmd"
