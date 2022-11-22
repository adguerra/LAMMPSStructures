from setuptools import setup, find_packages

VERSION = '0.0.1' 
DESCRIPTION = 'Python Wrapper for LAMMPS'
LONG_DESCRIPTION = 'Use this to create a lammps input file with all of the pieces you want to simulate'

# Setting up
setup(
        name="lammpsWithPython", 
        version=VERSION,
        author="Arman Guerra",
        author_email="<adguerra@bu.edu>",
        description=DESCRIPTION,
        long_description=LONG_DESCRIPTION,
        packages=find_packages(),
        install_requires=[], # add any additional packages that 
        # needs to be installed along with your package. Eg: 'caer'
        
        keywords=['python', 'first package'],
        classifiers= [
            "Development Status :: 3 - Alpha",
            "Intended Audience :: Education",
            "Programming Language :: Python :: 2",
            "Programming Language :: Python :: 3",
            "Operating System :: MacOS :: MacOS X",
            "Operating System :: Microsoft :: Windows",
        ]
)
