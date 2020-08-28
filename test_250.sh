#!/bin/bash
#$ -cwd
##### Specify job name
#$ -N test_250
##### Output file
#$ -o testo_250
##### Error file
#$ -e teste_250
##### number of cores
#$ -pe smp 1
##### memory per core
#$ -l mem_free=2G
##### When email notifications
#$ -m bea
##### Email address
#$ -M llanowaremissary@gmail.com

module load matlab

matlab -nosplash -nodesktop MISS_demo_Tasic_Wynton_250.m  

