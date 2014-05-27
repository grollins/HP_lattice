from sys import platform, prefix, version_info
from os.path import join, expanduser
from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext

from Scientific import N
assert N.package == "NumPy"

compile_args = []
include_dirs = ['.',]

compile_args.append("-DNUMPY=1")
if platform == 'win32':
    include_dirs.append(
        join(prefix,"Lib/site-packages/numpy/core/include")
    )
else:
    include_dirs.append(
        join(prefix, "lib/python%s.%s/site-packages/numpy/core/include" % \
                      version_info [:2])
    )

extensions = [Extension("vec2coords", ["vec2coords.pyx"],
                extra_compile_args=compile_args,
                include_dirs=include_dirs)]

setup(
    cmdclass={'build_ext': build_ext},
    ext_modules=extensions
)
