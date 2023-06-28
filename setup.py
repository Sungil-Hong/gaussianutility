import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="gaussianutility",
    version="0.2.0",
    author="Sungil Hong",
    author_email="suh33@pitt.edu",
    description="Useful utilities to use Gaussian09/16 software",
    long_description=long_description,
    long_description_content_type="text/markdown",
#    entry_points={'console_scripts': [
#                'xgeom = gaussianutility.xgeom:main',
#                'pathFeature = gaussianutility.pathFeature:main',
#                'genFeat = gaussianutility.genFeat:main',
#                'freezeLayer = gaussianutility.freezeLayer:main',
#                'inputSort = gaussianutility.inputSort:main',
#                ]},
    url="https://github.com/Sungil-Hong/gaussianutility",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
    install_requires=['numpy>=1.17.2',
                      'pandas>=1.5.3',
                      'periodictable>=1.6.1',
                      'ase>=3.22.1',
                      'math',
                      'sys',
                      'argparse>=1.1',
                      'matplotlib>=3.7.1',
                      'itertools',
                      'mendeleev>=0.14.0']
)
