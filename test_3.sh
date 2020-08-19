#!/bin/bash
#$ -cwd
##### Specify job name
#$ -N test_3
##### Output file
#$ -o testo_3
##### Error file
#$ -e teste_3
##### number of cores
#$ -pe smp 1
##### memory per core
#$ -l mem_free=2G
##### When email notifications
#$ -m bea
##### Email address
#$ -M llanowaremissary@gmail.com

module load matlab

matlab -nosplash -nodesktop MISS_demo_Zeisel_3.m  

