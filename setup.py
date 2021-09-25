import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="eeisp",
    version="0.3.0",
    license="GPL3.0",
    install_requires=['numpy','pandas','scipy','math','multiprocessing','time'],
    author="Natsu Nakajima, Ryuichiro Nakato",
    author_email="rnakato@iqb.u-tokyo.ac.jp",
    description="Identify gene pairs that are codependent and mutually exclusive from single-cell RNA-seq data.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/nakatolab/EEISP",
    keywords="eeisp scRNA-seq",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GPL-3.0 License",
        "Operating System :: OS Independent",
    ],
    python_requires=">=3.6",
)
