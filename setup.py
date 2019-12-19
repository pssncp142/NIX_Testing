from setuptools import setup, find_packages
from distutils.command.build import build
from subprocess import call
import os



scripts=['scripts/'+script for script in os.listdir('scripts')]
print scripts

class customBuild(build):

    def run(self):

        configure_cmd = ["./configure", 
                "--prefix=/home/ydallilar/.software/local"
                "--with-cpl=/home/ydallilar/.software/local",
                "--with-gsl=/home/ydallilar/.software/local"]
        build_cmd = ["make"]

        hdrl_path = "./hdrldemo-1.2.0"
        target = hdrl_path + "/hdrldemo/.libs/libhdrldemo.so"

        print target

        def configuref():
            call(configure_cmd, cwd=hdrl_path)
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

