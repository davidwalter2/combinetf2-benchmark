#!/usr/bin/env bash

# Read target number of CPUs (default 8)
NUM_CPUS="${1:-8}"
LAST_CPU=$((NUM_CPUS - 1))
CPU_RANGE="0-$LAST_CPU"
MEMS="0-1"

CG=/sys/fs/cgroup/mycpulimit

echo "[+] Creating/initialising $CG for CPUs $CPU_RANGE"

sudo mkdir -p "$CG"
sudo bash -c "echo +cpuset > /sys/fs/cgroup/cgroup.subtree_control || true"
sudo bash -c "echo $CPU_RANGE > $CG/cpuset.cpus"
sudo bash -c "echo $MEMS > $CG/cpuset.mems"
sudo chown -R $USER:$USER "$CG"

echo "[+] Moving shell (PID $$) into $CG"
sudo bash -c "echo $$ > $CG/cgroup.procs"

echo "[+] Now restricted to CPUs:"
taskset -cp $$
