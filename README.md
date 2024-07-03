# Modiview

[![Create env and Test](https://github.com/jaroslav-zukov/modiview/actions/workflows/build.yml/badge.svg?branch=master)](https://github.com/jaroslav-zukov/modiview/actions/workflows/build.yml)
[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)

### Prerequisites
Conda is installed on the local machine. You are running commands in the same directory as README file. 

### Installation guide
1. Create a new conda environment
```
conda env create -f environment.yml -n modiview -y
```
2. Activate conda 
```
conda activate modiview
```
3. Run the application
```
python modiview/app.py
```

### Supported modifications
- 5-Methylcytosine
- 5-Hydroxymethylcytosine
- N(4)-methylcytosine
- 6-Methyladenine

### Usage
Try uploading the `calls.bam` file from the repo or your own file with supported modifications.

### Future development
- Add more modifications
- Add more tests
- Possibility to add multiple bam files from different modifications calling models on the same file

### For developers
For code formatting use [black](https://github.com/psf/black)