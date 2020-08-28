#!/bin/bash
#$ -cwd
##### Specify job name
#$ -N test_50
##### Output file
#$ -o testo_50
##### Error file
#$ -e teste_50
##### number of cores
#$ -pe smp 1
##### memory per core
#$ -l mem_free=2G
##### When email notifications
#$ -m bea
##### Email address
#$ -M llanowaremissary@gmail.com

module load matlab

matlab -nosplash -nodesktop MISS_demo_Tasic_Wynton_50.m  

