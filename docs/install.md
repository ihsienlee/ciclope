# Installation
Install using pip. The flag `[all]` will install optional dependencies needed to run full pipelines and examples.
```shell
pip install ciclope[all]
```
---

## Development installation
Create and activate a virtual environment for development:
```shell
conda env create -n ciclope
conda activate ciclope
```
Clone the git repository:
```shell
git clone https://github.com/gianthk/ciclope.git
```
Navigate to the repo folder and install the package using pip:
```shell
cd ciclope
pip install .
```
Or install the package locally for development with [Flit](https://flit.pypa.io/en/latest/index.html):
```shell
flit install --symlink
```
---
## Build instructions
Build package:
```shell
python3 -m build
```
Upload to pypi:
```shell
python3 -m twine upload --repository pypi dist/*
```
Uninstall the existing build:
```shell
python -m pip uninstall ciclope
```