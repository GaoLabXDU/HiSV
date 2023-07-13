
"""
Setup script for HiSV.

"""
import os, sys, glob
import setuptools

def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()

if (sys.version_info.major!=3) or (sys.version_info.minor!=7):
    print('PYTHON 3.7 IS REQUIRED. YOU ARE CURRENTLY USING PYTHON {}'.format(sys.version.split()[0]))
    sys.exit(2)

# Guarantee Unix Format
for src in glob.glob('scripts/*'):
    text = open(src, 'r').read().replace('\r\n', '\n')
    open(src, 'w').write(text)

setuptools.setup(
    name = 'hisv',
    version = "1.0.7",
    author = "Li Junping",
    author_email = 'lijunping02@qq.com',
    url = 'https://github.com/GaoLabXDU/HiSV',
    description = 'A computational pipeline for structural variation detection from Hi-C data',
    keywords = ("Hi-C", "structural variation", "control-free"),
    scripts = glob.glob('scripts/*'),
    packages = setuptools.find_packages(),
    include_package_data = True,
    platforms = "any",
    license="MIT Licence",
    install_requires = [
        "numpy==1.21.6",
        "pandas==1.3.5",
        "pyBigWig==0.3.22",
        "prox-tv==3.3.0",
        "pysam==0.21.0",
        "pytest==7.4.0",
        "pygam==0.8.0",
        "scikit-learn==1.0.2",
        "cooler==0.9.2",
        ]
    )
