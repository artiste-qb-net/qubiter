from setuptools import setup, find_packages
from os import path
from io import open

with open("README.md", "r") as fh:
    long_description = fh.read()

setup(
    name="qubiter",
    version="1.0.0",
    author="qubiter contributors",
    keywords = ('quantum compiler'),
    author_email="rrtucci@gmail.com",
    description="Python tools for reading, writing, compiling, simulating quantum computer circuits.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/artiste-qb-net/qubiter",
    packages=find_packages(),
    #include_package_data=True,
    package_data={
        '': ['*.txt', '*.pdf', '*.ipynb', '*.md', '*.tex', '*.sty'],
        'qubiter': ['quantum_CSD_compiler/*.md'],
    },	
    data_files=[('qubiter', ['go_to_qubiter_inner_notebooks.ipynb'])],
    install_requires=[
        'numpy',
        'scipy'
    ],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: BSD License",
        "Operating System :: OS Independent",
    ],
)