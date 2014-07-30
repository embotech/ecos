xcopy /e "%RECIPE_DIR%\..\.." "%SRC_DIR%"
cd python
"%PYTHON%" setup.py install bdist_wheel
cd ..
if errorlevel 1 exit 1