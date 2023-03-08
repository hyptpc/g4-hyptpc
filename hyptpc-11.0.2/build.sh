#!/bin/sh

set -e
type cmake3 >/dev/null 2>&1 && alias cmake=cmake3

main_dir=$(dirname `readlink -f $0`)
obj_dir=$main_dir/obj

mkdir -pv $obj_dir
cd $obj_dir
#cmake .. -DCMAKE_INSTALL_PREFIX=$G4WORKDIR
cmake .. -DCMAKE_INSTALL_PREFIX=$main_dir
cmake --build . -- -j4
#cmake --install .
