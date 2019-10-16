from distutils.core import setup
import os


scripts=['scripts/'+script for script in os.listdir('scripts')]
print scripts

setup(
    name='NIX_Testing',
    version='0.1.0',
    author='Yigit Dallilar',
    author_email='ydallilar@mpe.mpg.de',
    py_modules=['NIX_Testing'],
    scripts=scripts, 
    url='https://github.com/pssncp142/NIX_Testing',
    description='A python module for analyzing NIX Testing data',
    long_description=open('README.md').read(),
    install_requires=[
        "numpy", "matplotlib", "sep", "astropy"
    ],
)

