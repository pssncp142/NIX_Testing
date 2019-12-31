from setuptools import setup, find_packages
from setuptools import Extension
from distutils.command.build import build
from subprocess import call
import os, numpy


include_dirs = [numpy.get_include(), '/usr/include/', '/home/ydallilar/.software/local/include/', 'hdrldemo-1.2.0/hdrl/']
library_dirs = ['/usr/lib/', '/home/ydallilar/.software/local/lib/', 'hdrl/.lib/']
libraries = ['cplcore', 'cpldrs', 'cplui', 'cpldfs', 'cext', 
        'fftw3', 'wcs', 'm', 'gsl', 'gslcblas']
extra_link_args = ['-lgomp']
extra_compile_args = ['-fopenmp']
define_macros = [('HDRL_USE_EXPERIMENTAL', None), ('PIC', None)]
extra_objects = ['hdrldemo-1.2.0/hdrl/.libs/libhdrl.a']

# define the extension module
hdrl2 = Extension('NIX_Testing.HDRL2', sources=['NIX_Testing/HDRL2.c'],  include_dirs=include_dirs, 
        libraries=libraries, library_dirs=library_dirs, 
        extra_compile_args=extra_compile_args, extra_link_args=extra_link_args,
        define_macros=define_macros, extra_objects=extra_objects)


scripts=['scripts/'+script for script in os.listdir('scripts')]
print(scripts)

class customBuild(build):

    def run(self):

        configure_cmd = ["./configure", 
                "--prefix=/home/ydallilar/.software/local"
                "--with-cpl=/home/ydallilar/.software/local",
                "--with-gsl=/home/ydallilar/.software/local"]
        build_cmd = ["make"]
        hack_cmd = ["sed", "-i", "s/hdrldemo_utils/hdrldemo_python/g", "hdrldemo/Makefile"]

        hdrl_path = "./hdrldemo-1.2.0"
        target = hdrl_path + "/hdrldemo/.libs/libhdrldemo.so"

        print(target)

        def configuref():
            call(configure_cmd, cwd=hdrl_path)
        def hack():
            call(hack_cmd, cwd=hdrl_path)
        def buildf():
            call(build_cmd, cwd=hdrl_path)

        configuref()
        buildf()

        self.copy_file(target, './NIX_Testing/lib')

        build.run(self)

setup(
    name='NIX_Testing',
    version='0.1.0',
    author='Yigit Dallilar',
    author_email='ydallilar@mpe.mpg.de',
    ext_modules = [hdrl2],
    packages = ['NIX_Testing'],
    package_data={'NIX_Testing': ['lib/libhdrldemo.so']},
    scripts=scripts, 
    url='https://github.com/pssncp142/NIX_Testing',
    description='A python module for analyzing NIX Testing data',
    long_description=open('README.md').read(),
    install_requires=[
        "numpy", "matplotlib", "sep", "astropy"
    ],
    cmdclass = {'build': customBuild}
)

