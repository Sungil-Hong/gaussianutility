import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="gaussianutility",
    version="0.1.6",
    author="Sungil Hong",
    author_email="s.hong@pitt.edu",
    description="Useful utilities to use Gaussian09/16 software",
    long_description=long_description,
    long_description_content_type="text/markdown",
    entry_points={
        'console_scripts': [
            'com2vasp=gaussianutility.com2vasp:main',
            'com2xyz=gaussianutility.com2xyz:main',
            'xyz2com=gaussianutility.xyz2com:main',
            'freezeLayer=gaussianutility.freezeLayer:main',
            'gibbsTemp=gaussianutility.gibbsTemp:main',
            'out2com=gaussianutility.out2com:main',
            'printE=gaussianutility.printE:main',
            'sortInput=gaussianutility.sortInput:main',
            'spectrum=gaussianutility.spectrum:main',
            'vasp2com=gaussianutility.vasp2com:main'
        ],
    },
    url="https://github.com/Sungil-Hong/gaussianutility",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent"],
    python_requires='>=3.6',
    install_requires=['numpy>=1.17.2',
                      'pandas>=1.1.0',
                      'periodictable>=1.6.1',
                      'argparse>=1.1',
                      'matplotlib>=3.5.1',
                      'mendeleev>=0.9.0']
)
