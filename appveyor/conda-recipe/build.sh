#!/bin/sh
cp -r $RECIPE_DIR/../.. $SRC_DIR
cd python
$PYTHON setup.py install bdist_wheel
cd ..