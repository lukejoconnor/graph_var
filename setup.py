from setuptools import setup, find_packages

setup(
    name="graph_var",
    version="0.1.0",
    packages=find_packages(),
    install_requires=[
        "networkx~=3.2.1",
        "numpy~=1.26.4",
        "scipy~=1.12.0",
        "setuptools~=68.2.0",
        "tqdm",
        "jupyter",
        "matplotlib"
    ],
    entry_points={
        'console_scripts': [
            'graph_var=graph_var.cli:main',
        ],
    },
    author="",
    author_email="",
    description="A tool for analyzing genetic variants in pangenome graphs",
    long_description=open("README.md").read(),
    long_description_content_type="text/markdown",
    python_requires=">=3.6",
)
