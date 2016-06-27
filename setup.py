from setuptools import setup

setup(name = 'LE4PD',
    version = '0.1',
    description = 'Python API for Protein Dynamics using the Langevin Formalism (LE4PD)',
    #url = 'https://github.com/pgromano/LE4PD',
    packages = ['LE4PD'],
    install_requires=[
        'numpy',
        'scipy',
        'matplotlib',
        'mdtraj'
    ],
    zip_safe = False
)
