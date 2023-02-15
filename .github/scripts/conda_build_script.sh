#!/bin/bash

set -ex
set -o pipefail


check_if_setup_file_exists() {
    if [ ! -f setup.py ]; then
        echo "setup.py must exist in the directory that is being packaged and published."
        exit 1
    fi
}

check_if_meta_yaml_file_exists() {
    if [ ! -f meta.yaml ]; then
        echo "meta.yaml must exist in the directory that is being packaged and published."
        exit 1
    fi
}

build_package(){
    # Build for Linux
    conda build -c conda-forge --output-folder . .
}

upload_package(){
    export ANACONDA_API_TOKEN=$INPUT_ANACONDATOKEN
    anaconda upload --label main noarch/*.tar.bz2
}

check_if_setup_file_exists
cd conda
check_if_meta_yaml_file_exists
build_package
upload_package