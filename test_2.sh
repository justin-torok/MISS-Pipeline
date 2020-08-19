#!/bin/bash
#$ -cwd
##### Specify job name
#$ -N test_2
##### Output file
#$ -o testo_2
##### Error file
#$ -e teste_2
##### number of cores
#$ -pe smp 1
##### memory per core
#$ -l mem_free=2G
##### When email notifications
#$ -m bea
##### Email address
#$ -M llanowaremissary@gmail.com

module load matlab

matlab -nosplash -nodesktop MISS_demo_Zeisel_2.m  
