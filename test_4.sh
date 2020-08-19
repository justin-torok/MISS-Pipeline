#!/bin/bash
#$ -cwd
##### Specify job name
#$ -N test_4
##### Output file
#$ -o testo_4
##### Error file
#$ -e teste_4
##### number of cores
#$ -pe smp 1
##### memory per core
#$ -l mem_free=2G
##### When email notifications
#$ -m bea
##### Email address
#$ -M llanowaremissary@gmail.com

module load matlab

matlab -nosplash -nodesktop MISS_demo_Zeisel_4.m  

