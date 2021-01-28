# Calculating HWE and finding outliers of given population
Nan Hu / Jan 24, 2020

---
## Aims
1. Starting to work on HPCC
2. Building up software environment
3. Getting familiar to VCF fire format
4. Using vcftools and/or angsd to calculate HWE
5. Result display and interpretation 
---

## Environment preparation
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
If this command gives you correct information, this means you successfully created software environment for this practice.

## Calculating HWE and do Fisher's Exact Test
The datasets we use for calculating HWE are from three species: *Salix phlebophylla*, *Salix nivalis*, and *Salix reticulata*. We prepared the data through initial read alignment, genotyping and hard filtering. So far, these data are prepared for Thus, the data is ready for most population
