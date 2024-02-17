"""
MQuad - Mixture Modelling for Mitochondrial Mutation detection in single-cell 
omics data
See: https://github.com/aaronkwc/MQuad
"""

# Always prefer setuptools over distutils
from setuptools import setup, find_packages
# To use a consistent encoding
from codecs import open
from os import path

here = path.abspath(path.dirname(__file__))

# Set __version__ for the project.
exec(open("./mquad/version.py").read())

# Get the long description from the relevant file
with open(path.join(here, 'README.rst'), encoding='utf-8') as f:
    long_description = f.read()
    
reqs = ['numpy>=1.9.0', 'scipy>=1.4.0', 'matplotlib', 'BBMix>=0.2.1', 'vireoSNP', 'kneed',
        'pandas', 'seaborn', 'scikit-learn']

setup(
    name='mquad',

    # Versions should comply with PEP440.  For a discussion on single-sourcing
    # the version across setup.py and the project code, see
    # https://packaging.python.org/en/latest/single_source_version.html
    version=__version__,

    description='MQuad - Mixture Modelling for Mitochondrial Mutation detection',
    long_description=long_description,

    # The project's main homepage.
    url='https://github.com/aaronkwc/MQuad',

    # Author details
    author='Aaron Kwok',
    author_email='aaronkwc@connect.hku.hk',

    # Choose your license
    license='Apache-2.0',

    # What does your project relate to?
    keywords=['mitochondrial mutation', 'single-cell omics data', 
              'binomial mixture model', 'model selection'],

    # You can just specify the packages manually here if your project is
    # simple. Or you can use find_packages().
    packages=find_packages(),

    entry_points={
          'console_scripts': [
              'mquad = mquad.mquad_CLI:main',
              ],
          }, 

    # List run-time dependencies here.  These will be installed by pip when
    # your project is installed. For an analysis of "install_requires" vs pip's
    # requirements files see:
    # https://packaging.python.org/en/latest/requirements.html
    
    install_requires=reqs,

    extras_require={
        'docs': [
            #'sphinx == 1.8.3',
            'sphinx_bootstrap_theme']},

    py_modules = ['mquad']

    # buid the distribution: python setup.py sdist
    # upload to pypi: twine upload dist/...
)
