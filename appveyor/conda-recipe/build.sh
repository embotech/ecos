#!/bin/sh
cp -r $RECIPE_DIR/../.. $SRC_DIR
$PYTHON setup.py install bdist_wheel