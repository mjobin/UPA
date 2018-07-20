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
    subp = subprocess.Popen(['/bin/bash', '-c', cmd], stdout=PIPE, stderr=PIPE)
    stdout, stderr = subp.communicate()
    return stdout

def haplogrep_gen_hsd(flist):
    flength = len(flist)
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

        depthinfo = bash_command("samtools depth " + sample + ".bam")

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

    mtmcmd = "bcftools mpileup -I -d 8000 -Ou -f " + ref + " "
    for fbam in fbamlist:
        mtmcmd = mtmcmd + fbam + " "

    mtmcmd = mtmcmd + "| bcftools call -V indels --ploidy 1 -Ou -m -v -o " + bcname + ".vcf"

    bash_command(mtmcmd)
    return bcname + ".vcf"



def haplogrep_java(invcf):

    filebase, filext = os.path.splitext(invcf)

    firstfile = filebase + "-FIRST.hsd"

    bash_command("java -jar /data/scripts/haplogrep-2.1.1.jar --format vcf --in " + invcf + " --out " + firstfile + " --phylotree 17")

    # then edit that HSD file and do it again
    hsdoutlines = []
    hsdfirstfile = open(firstfile, 'r')
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

    secondfile = filebase + "-FIRST.hsd"

    hsdsecondfile = open(secondfile, 'w')
    hsdsecondfile.write("SampleID\tRange\tHaplogroup\tPolymorphisms")
    hsdsecondfile.write("\n")
    for hsdoutline in hsdoutlines:
        hsdoutcols = hsdoutline.split("\t")
        if len(hsdoutcols[3]) > 0: # LEAVE OUT indivs for whom where are no viable regions
            hsdsecondfile.write(hsdoutline)
            hsdsecondfile.write("\n")
    hsdsecondfile.close()

    finalfile = filebase + "-FINAL.hsd"

    bash_command("java -jar /data/scripts/haplogrep-2.1.1.jar --format hsd --in " + hsdsecondfile + " --out " + finalfile + " --phylotree 17")



if __name__ == "__main__":

    print "\n****************\nMOMMY\n****************\n"
    print "So far, I'm just a module...."




    exit(1)
