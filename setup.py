from setuptools import setup, find_packages

setup(
    name='MAProt',
    version='0.1.0',
    description='python tools for proteomics data visualization',
    author='David Mast',
    author_email='david.h.mast@gmail.com',
    packages=find_packages(),  # List of packages to include
    install_requires=[
        'Bio',
        'biopython',
        'biopython',
        'dash_bio=',
        'matplotlib',
        'pyteomics',
        'Requests',
        'setuptools',
    ],       # List of dependencies
)