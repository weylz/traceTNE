#!/bin/bash
# File    : run_traceTNE.sh
# Time    : 2024/10/01 08:08:00
# Author  : Wenyong Zhu
# Version : 1.0.0
# Desc    : ...

export PATH=/path/to/traceTNE:$PATH
chmod 777 /path/to/traceTNE/*

traceTNE.sh -H

traceTNE.sh \
    -B /path/to/blocklist.bed \
    -S /path/to/(hg38.chrom.sizes or hg19.chrom.sizes or ...) \
    -I ProjectName(FolderName)
