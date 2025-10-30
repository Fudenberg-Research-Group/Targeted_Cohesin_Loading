import os
import io
from setuptools import setup, find_packages

VERSION = "0.1.0"
DESCRIPTION = "Dynamic barriers lattice translocators for loop extrusion simulations"


def _read(*parts, **kwargs):
    filepath = os.path.join(os.path.dirname(__file__), *parts)
    encoding = kwargs.pop("encoding", "utf-8")
    with io.open(filepath, encoding=encoding) as fh:
        text = fh.read()
    return text


def get_requirements(path):
    content = _read(path)
    return [
        req
        for req in content.split("\n")
        if req != "" and not (req.startswith("#") or req.startswith("-"))
    ]

install_requires = get_requirements("requirements.txt")

setup(
    name="Target_cohesin_loading",
    version=VERSION,
    description=DESCRIPTION,
    #url="https://github.com/Fudenberg-Research-Group/Target_cohesin_loading",
    url="https://github.com/hrahmanin/Target_cohesin_loading/",
    author="Hadi Rahmaninejad",
    author_email="rahmanin@usc.edu",
    license="MIT",
    packages=find_packages(include=["target_cohesin", "target_cohesin.*"]),
    install_requires=install_requires
)