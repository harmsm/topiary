#!/bin/bash

cd source
conda run -n topiary python ../script/build_api_docs.py
cd ..

rm -rf build
conda run -n topiary make html
