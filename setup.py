# Written usingn resources from:
# https://packaging.python.org/en/latest/distributing.html#working-in-development-mode
# https://github.com/pypa/sampleproject/blob/master/setup.py

from setuptools import setup, find_packages

setup(name = "pygenesSQLite",
      version = "0.0.11",
      py_modules = ["pygenesSQLite", "pygenesSQLiteScripts"],
      entry_points =  {
          "console_scripts" : [
              "pygenesSQLite=pygenesSQLiteScripts:main"
          ]
      }
)
