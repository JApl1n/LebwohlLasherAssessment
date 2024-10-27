from setuptools import setup, Extension
from Cython.Build import cythonize
import numpy as np

setup(
    name='Cythonised Lebwohl Lasher Program',
    ext_modules=cythonize("cythonLebwohlLasher.pyx"),
    include_dirs=[np.get_include()], #ensures module loaded for compilation
)
