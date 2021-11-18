import setuptools

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setuptools.setup(
    name="py_nemo_rebuild",
    version="0.0.1",
    author="Piergiuseppe Fogli",
    author_email="piergiuseppe.fogli@cmcc.it",
    description="Rebuild NEMO/XIOS multiple output/restart files in a single file.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/CMCC-Foundation/py_nemo_rebuild",
    project_urls={
        "Bug Tracker": "https://github.com/pypa/sampleproject/issues",
    },
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GPL",
        "Operating System :: OS Independent",
    ],
    package_dir={"": "src"},
    packages=setuptools.find_packages(where="src"),
    python_requires=">=3.6",
    entry_points={
        'console_scripts': [
            'nemo_rebuild_py=py_nemo_rebuild.nemo_rebuild:rebuild',
        ],
    },
)
