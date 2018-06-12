#!/usr/bin/python

#####################################
#####        HPG Lab            #####
#####    updated May 2018       #####
#####       MJJ                 #####
#####################################

# Driver and merging script for phylogenetics and imputation/correction
# Requires 'barcode' file with fastq prefix, output name and internal barcodes
# Author: Matthew Jobin, UCSC Human Paleogenomics Lab

import argparse
from argparse import RawTextHelpFormatter
import os
import datetime
import progressbar
import subprocess
from subprocess import Popen, PIPE

def bash_command(cmd):
    cmdfile.write(cmd)
    cmdfile.write("\n\n")
    subp = subprocess.Popen(['/bin/bash', '-c', cmd], stdout=PIPE, stderr=PIPE)
    stdout, stderr = subp.communicate()
    if verbose:
        print stdout
    logfile.write(stdout)
    if verbose:
        print stderr
    logfile.write(stderr)
    return stdout


if __name__ == "__main__":

    print "\n****************\nTUBEAMP\n****************\n"



    parser = argparse.ArgumentParser(description="# This script:\n"
                                                    "1. does stuff.\n\t"
                                                 "- ", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-bc_file', metavar='<bc_file>', help='location of barcode files, Must have a newline at end.', default="")
    parser.add_argument('-bam_list', metavar='<bam_list', help='List of BAM files', default="")
    parser.add_argument('-wd', metavar='<wd>', help='Working directory. Defaults to current.', default='.')
    parser.add_argument('-verbose', dest='verbose', help='Print stdout and stderr to console.',
                        action='store_true')
    parser.set_defaults(verbose=False)
    parser.add_argument('-overwrite', dest='overwrite', help='Overwrite existing files and directories.',
                        action='store_true')
    parser.set_defaults(overwrite=False)



    #HAPLOGREP
    #Generate a consensus sequence using samtools mpileup
    fbamlist = []
    regdic = {}
    print "\nCreating mpileup consensus and writing as a VCF..."
    bar = progressbar.ProgressBar()
    for i in bar(range(flength)):
        sample = flist[i]

        fbamlist.append(sample + ".bam")
        # get list of depth 1 +
        #open file to read
        regionsfile = open(sample + "-regions.txt", 'w')
        regionsfile.write("CHROM\tPOS\tPOS_TO\n")
        regline = ""

        depthinfo = bash_command("samtools depth " + sample + ".bam", False)

        depthlines = depthinfo.split("\n")
        curstart = 0
        curend = 0
        lastpos = 0
        ingap = True
        firstthrough = True
        for dline in depthlines:
            dcols = dline.split("\t")
            if len(dcols) < 3:
                continue
            curpos = int(dcols[1])
            curdepth = int(dcols[2])
            if firstthrough:
                if curdepth >= mindepth:
                    curstart = curpos
                    firstthrough = False
            else:
                if curpos <= lastpos + maxgap: #Continuous run
                    if curdepth >= mindepth:
                        pass
                    else: #depth too small, stop here
                        curend = lastpos
                        regionsfile.write(dcols[0] + "\t" + str(curstart) + "\t" + str(curend) + "\n")
                        regline = regline + str(curstart) + "-" + str(curend) + ";"
                        curstart = curpos
                        ingap = True
                else: #discontinous
                    curend = lastpos
                    regionsfile.write(dcols[0] + "\t" + str(curstart) + "\t" + str(curend) + "\n")
                    regline = regline + str(curstart) + "-" + str(curend) + ";"
                    curstart = curpos
                    ingap = True
            lastpos = curpos
        regdic[sample + ".bam"] = regline



    bash_command("java -jar /data/scripts/haplogrep-2.1.1.jar --format vcf --in " + bcname + "-E.vcf.gz --out " + bcname + "-E.hsd --phylotree 17", verbose)

    # then edit that HSD file and do it again
    hsdoutlines = []
    hsdfirstfile = open(bcname + "-E.hsd", 'r')
    hsdheadline = "ERROR in HSD conversion"
    for hsdline in hsdfirstfile:
        hsdnewline = ""
        hsdcols = hsdline.strip().split("\t")
        if hsdcols[0] == "SampleID":
            hsdheadline = hsdline.strip()
        else:
            if len(regdic[hsdcols[0]]) > 1: # LEAVE OUT indivs for whom where are no viable regions
                for i in range(len(hsdcols)):
                    if i == 0 or i == 5:
                        hsdnewline = hsdnewline + hsdcols[i] + "\t"
                    elif i == 1:
                        hsdnewline = hsdnewline + regdic[hsdcols[0]] + "\t"
                    elif i == 2:
                        hsdnewline = hsdnewline + "?\t"
                    else:
                        pass
                hsdoutlines.append(hsdnewline)
    hsdfirstfile.close()

    hsdfinalfile = open(bcname + "-E-SECOND.hsd", 'w')
    hsdfinalfile.write("SampleID\tRange\tHaplogroup\tPolymorphisms")
    hsdfinalfile.write("\n")
    for hsdoutline in hsdoutlines:
        hsdoutcols = hsdoutline.split("\t")
        if len(hsdoutcols[3]) > 0: # LEAVE OUT indivs for whom where are no viable regions
            hsdfinalfile.write(hsdoutline)
            hsdfinalfile.write("\n")
    hsdfinalfile.close()

    bash_command("java -jar /data/scripts/haplogrep-2.1.1.jar --format hsd --in " + bcname + "-E-SECOND.hsd --out " + bcname + "-E-FINAL.hsd --phylotree 17", verbose)