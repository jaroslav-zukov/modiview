# Modiview

[![Create env and Test](https://github.com/jaroslav-zukov/modiview/actions/workflows/build.yml/badge.svg?branch=master)](https://github.com/jaroslav-zukov/modiview/actions/workflows/build.yml)

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

### Usage
Try uploading the `calls.bam` file from the repo or your own file with C+h, C+m modifications.

### Future development
- Add more modifications
- Add more tests
- Possibility to add multiple bam files from different modifications calling models 