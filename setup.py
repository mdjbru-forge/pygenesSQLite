# Written usingn resources from:
# https://packaging.python.org/en/latest/distributing.html#working-in-development-mode
# https://github.com/pypa/sampleproject/blob/master/setup.py

from setuptools import setup, find_packages

setup(name = "pygenes",
      version = "0.0.4",
      py_modules = ["pygenes", "pygenesScripts"],
      entry_points =  {
          "console_scripts" : [
              "pygenes=pygenesScripts:main"
          ]
      }
)
