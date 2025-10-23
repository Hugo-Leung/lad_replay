#!/bin/bash

SCRIPT_PATH="$(dirname "$(readlink -f ${BASH_SOURCE[0]})")"
TAR_PATH_PREFIX="${SCRIPT_PATH}/../lad_replay.tar.gz"
TAR_PATH=${TAR_PATH_PREFIX}
i=0
while [ -e "${TAR_PATH}" ]; do
  i=$((i+1))
  TAR_PATH="${TAR_PATH_PREFIX}.$i"
done

cd ${SCRIPT_PATH}
tar --exclude='*.root' -czf ${TAR_PATH} .
cd -
