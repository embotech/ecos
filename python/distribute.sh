#!/bin/sh
rm -rf external include src
mkdir include src external external/amd external/ldl
cp -r ../include .
cp -r ../external/amd/src external/amd/
cp -r ../external/ldl/src external/ldl/
cp -r ../src .
cp -r ../external/amd/include external/amd/
cp -r ../external/ldl/include external/ldl/
cp -r ../external/SuiteSparse_config external/
cp setup.py tmp.py
cp .setup_dist.py setup.py
sudo python setup.py sdist upload
mv tmp.py setup.py
