from setuptools import Extension, setup
from Cython.Build import cythonize

extensions = [
    Extension ( "pyathad",
              ['pyathad.pyx', 'PythonStream.cpp'],
              language = 'c++',
              extra_compile_args=["-std=c++17"],
              libraries = ['athad'],
              include_dirs = ['../atmosphere', '../lib', '../tinyxml2'],
              library_dirs = ['..'],
              extra_link_args = ['-fopenmp'],
              )]

setup(
    name = 'pyathad',
    ext_modules = cythonize(extensions, language_level = "2"),
    )
