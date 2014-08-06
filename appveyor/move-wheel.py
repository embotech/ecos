import os
import shutil
import glob
from conda_build import config

wheel_path = os.path.join(config.croot, "work\python\dist\*.whl")
wheel = glob.glob(wheel_path)[0]

shutil.move(wheel, '.')