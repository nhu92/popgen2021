# Calculating HWE and finding outliers of given population
Nan Hu / Jan 24, 2020

---
## Aims
1. Starting to work on HPCC
2. Building up software environment
3. Getting familiar to VCF file format
4. Using vcftools and/or angsd to calculate HWE
5. Result visulization and interpretation 
---

## Environment Setup
Before starting working on course project. We need to get familiar with working with super computer clusters and Linux environment. TTU has provided HPCC (High Performance Computer Center) services. In order to get access to HPCC resources, we need to [apply](https://www.depts.ttu.edu/hpcc/accounts/studentrequest.php) for an account for research use. 

Software for this class mostly run under Linux operating system, which is the OS of HPC. On the website of HPCC, you could find some tutorials of using HPC and basics of Linux. They also offers regular [training](https://www.depts.ttu.edu/hpcc/about/training.php) on these topics. There are previous training records to help you quickly get into  bioinformatics fields.

When you are ready to use Linux for this course, this practice requires installation of software environment called [conda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html). This is an open-source package management software helps installing other analytical tools and prepare required environment. It is simple to use Conda to install packages by just typing: ```conda install <package name>```. For more simple usages, you can visit [conda cheat sheet](https://docs.conda.io/projects/conda/en/latest/user-guide/cheatsheet.html).

Now, we can start using conda to install one useful software called [vcftools](http://vcftools.sourceforge.net/). Type below command to prompt:
```bash
conda install vcftools
```
We should see the software searching for packages. When system asks you if you want to install packages, type 'y' and hit enter. Then, we just wait it automatically finishing installation.

When finished isntallation, try to call help document to test if the software works. Sometimes, you need to relogin to remote computers to let some source files reloaded.
```bash
vcftools --help
```
If this command gives you correct information, it means you successfully created software environment for this practice.

Another software which is quite often used in population genetics is called 'angsd'. It is not available to install through conda yet (I am not sure. If I am wrong please tell me). The installation notes are [here](https://github.com/ANGSD/angsd/blob/master/README.md).

> Some of you having trouble in installing 'angsd' on HPCC. The author of the software offers the local installation processes which will install both htslib and angsd. However, if you still stuck in install the package. Here are some steps I take to install it and it works. 
```bash
# git clone 'angsd' software package from github:
git clone https://github.com/angsd/angsd.git

# download newest htslib:
wget https://github.com/samtools/htslib/releases/download/1.11/htslib-1.11.tar.bz2

# decompress it
tar -xf htslib-1.11.tar.bz2

# load gcc compiler on HPCC
module load gcc/10.1.0

# install htslib manually
cd htslib-1.11
./configure
make

# back to angsd directory and install angsd
cd ../angsd
make HTSSRC=../htslib-1.11

# test if it works
angsd --help
```

## Calculating HWE and do Fisher's Exact Test
The datasets we use for calculating HWE are from three species: *Salix phlebophylla*, *Salix nivalis*, and *Salix reticulata*. We prepared the data through initial read alignment, genotyping and hard filtering. So far, these data are prepared for this practice but might not be approporate for later analysis. There would be more filtering processes applied to these data in the future.

Current data file for three species are stored in HPCC directory: ```/lustre/scratch/nhu/popgen2021/raw/```. (*It is not there right now but I promise I will finish prepare those and get it available before midnight of Jan. 29, 2021 - Nan*)(It is done now! - Nan). You can copy them to your own working directory for easier access.

Now, create a text file and type below commands into. Then, submit it.
> For beginners, use `nano HWE_vcftools.sh` to create a job submission script and copy following scripts block to it. After edit it, hit Ctrl+X and Enter to save the file.
```bash
#!/bin/bash
#SBATCH -J HWEcalc
#SBATCH -o %x.o%j
#SBATCH -e %x.e%j
#SBATCH -p nocona
#SBATCH -N 1
#SBATCH -n 32

vcftools --vcf <vcf file> --hardy --out <output file name>
```
> HPCC uses new job submission system this year. You should use `sbatch HWE_vcftools.sh` to submit it. Use `squeue --me` or `squeue -u <your eraider>` to check job status. Use `scancel <job ID>` to cancel failed jobs.

The result of this will reports a p-value for each site from a Hardy-Weinberg Equilibrium test (as defined by Wigginton, Cutler and Abecasis (2005)). It will also contain the Observed numbers of Homozygotes and Heterozygotes and the corresponding Expected numbers under HWE.

Alternatively, you can use 'angsd' to calculate the HWE, too. Submission scripts of angsd are:
```bash
#!/bin/bash
#SBATCH -J HWEcalc
#SBATCH -o %x.o%j
#SBATCH –e %x.e%j
#SBATCH -p nocona
#SBATCH -N 1
#SBATCH -n 32

angsd -vcf-gl <vcf file> -doHWE 1
```
Similar results should be generated by angsd. 
> I am receiving errors about exceeding memories for angsd running on RedRaider cluster.

## Data transferring
After having processed commands on HPCC, you will see your result file under your working directory. For example, if we use vcftools to calculate HWE, the result will have a suffix of ".HWE". Now, view the file with command:
```bash
less <HWE result>
```
The expected results will be in a TAB spaced table form, with 8 columns recording HWE information per segregating site. In order to do further analysis like filtering, statistical analysis, and visualizing, some of you may want your data to be on your local computer. There are several ways to transferring data across HPCC and our local computer. 
### Transferring with Globus (recommended)
HPCC officially suggests us to use Globus to transferring data. Here is a tutorial from HPCC. It requires to install personal Globus client.
https://www.depts.ttu.edu/hpcc/userguides/general_guides/file_transfer.php

### Transferring with file transferring software
A lot of software support data transferring across remote computer and local computer. I use [xftp](https://www.netsarang.com/en/xftp/) and it is free for researchers. You can find any software having this function to use.

Some of you suggest [filezilla](https://filezilla-project.org/). It looks the similar thing like xftp.

### Transferring with 'scp' command (not recommended)
There is a simple way to download data from HPCC to local computer. On windows system, type 'cmd' on search box to call the command line. On Mac, use terminal window. We can use following command to copy data from remote computer:
```bash
scp <eraider>@login.hpcc.ttu.edu:<Remote file> <local directory>
```
Here is an example:
```bash
scp nhu@login.hpcc.ttu.edu:/home/nhu/snivalis.HWE .
cwd
```
And it will ask you about your eraider password. The file will download to your current directory on you local computer. 'cwd' command is to see the current working directory.

This way to transferring data is quick but not recommended by HPCC since when it transferring large data, it is not quite reliable and might cause issues on login node of HPCC.

## Plotting HWE
See [HWEplot.R](https://github.com/gudusanjiao/popgen2021/blob/main/HWEplot.R).
