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

### Supported modifications
- 5-Methylcytosine
- 5-Hydroxymethylcytosine
- N(4)-methylcytosine
- 6-Methyladenine

### Usage
Try uploading the `calls.bam` file from the repo attachments or your own file with supported modifications.

### Preview of the application
Preview the whole read
![Preview](attachments/preview_whole.png)
Zoom in on specific read region
![Zoom](attachments/preview_zoomed.png)


### Future development
- Add more modifications
- Add more tests
- Possibility to add multiple bam files from different modifications calling models on the same file