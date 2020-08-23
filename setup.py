import setuptools

import synple

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="synple", # Replace with your own username
    version="0.2",
    author="Carlos Allende Prieto",
    author_email="callende@iac.es",
    description="A simple package to compute synthetic stellar spectra",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/callendeprieto/synple",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
)
