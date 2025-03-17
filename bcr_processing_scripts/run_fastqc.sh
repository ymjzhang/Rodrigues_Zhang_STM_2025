#!/bin/sh
###########
echo $1
module load fastqc
fastqc $1
