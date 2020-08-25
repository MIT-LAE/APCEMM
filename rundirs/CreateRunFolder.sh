#!/bin/bash

for folder in "$@"
do
    mkdir "$folder"
    cp SampleRunDir/{Batch_APCEMM.sh,input.apcemm,Makefile,species.dat} $folder
    echo "Successfully copied to $folder"
done
