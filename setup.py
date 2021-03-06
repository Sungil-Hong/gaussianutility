import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="gaussianutility",
    version="0.1.5",
    author="Sungil Hong",
    author_email="suh33@pitt.edu",
    description="Useful utilities to use Gaussian09 software",
    long_description=long_description,
    long_description_content_type="text/markdown",
    entry_points={'console_scripts': [
                'xgeom = gaussianutility.xgeom:main',
                'pathFeature = gaussianutility.pathFeature:main',
                'genFeat = gaussianutility.genFeat:main',
                'freezeLayer = gaussianutility.freezeLayer:main',
                'oniomSort = gaussianutility.oniomSort:main'
                ]},
    url="https://github.com/Sungil-Hong/gaussianutility",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
    install_requires=['numpy>=1.17.2']
)
