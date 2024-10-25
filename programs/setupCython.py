from setuptools import setup
from Cython.Build import cythonize

setup(
        name='Cythonised Lebwohl Lasher Program',
        ext_modules=cythonize("cythonLebwohlLasher.pyx"),
)
