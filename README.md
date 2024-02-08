# CompMutTB - A framework to identify compensatory mutations from *Mycobacterium tuberculosis* whole genome sequences
## Setup
To setup and install dependencies for CompMutTB, please use the following code. Clone the github repository and install dependencies into a conda environment using the CompMutTB.yml file and setup Rscript. 
```
conda create -n CompMutTB
git clone https://github.com/NinaMercedes/CompMutTB.git
conda activate CompMutTB
cd CompMutTB
conda env update --file conda/CompMutTB.yml
Rscript conda/setup.R
conda deactivate
```
