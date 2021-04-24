
from setuptools import find_packages, setup

with open("requirements.txt", "r") as f:
    required = f.read().splitlines()

with open("README.md", "r") as f:
    long_description = f.read()

setup(
    name="joy_plotter",
    version="0.1.0",
    # packages=["joy_plotter"],
    packages=find_packages(),
    url="https://github.com/josephwkania/joy-plotter",
    install_requires=required,
    long_description=long_description,
    long_description_content_type="text/markdown",
    author="Joseph W Kania",
    author_email="60994292+josephwkania@users.noreply.github.com",
    description="Makes Joy Division like plots from your .h5 cand files",
    # entry_points={"console_scripts": ["joy_plotter=joy_plotter.joy_plotter:main", ], },
    scripts=["joy_plotter/joy-plotter.py", "joy_plotter/line_plotting.py"],
    tests_require=["pytest", "pytest-cov"],
    classifiers=[
        "Programming Language :: Python :: 3.6",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Topic :: Documentation",
        "License :: OSI Approved :: MIT License",
        "Natural Language :: English",
    ],
)
