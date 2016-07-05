from setuptools import setup

setup(
    name = 'LE4PD',
    version = '0.2',
    author = 'Pablo Romano',
    description = 'Python API for Protein Dynamics using the Langevin Formalism (LE4PD)',
    url = 'https://github.com/GuenzaLab/LE4PD',

    packages = ['LE4PD'],
    install_requires=[
        'numpy',
        'scipy',
        'matplotlib',
        'mdtraj',
        'openmm'
    ],
    zip_safe = False
)
