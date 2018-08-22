#!/usr/bin/python

#####################################
#####        HPG Lab            #####
#####    updated July 2018      #####
#####       MJJ                 #####
#####################################

# File format conversion and utility script for use with UPA
# Author: Matthew Jobin, UCSC Human Paleogenomics Lab
import shutil
import subprocess
import os
from subprocess import PIPE

def bash_command(cmd, verbose, cmdfile, logfile):
    cmdfile.write(cmd)
    cmdfile.write("\n\n")
    subp = subprocess.Popen(['/bin/bash', '-c', cmd], stdout=PIPE, stderr=PIPE)
    stdout, stderr = subp.communicate()
    subp.wait()
    if verbose:
        print stdout
    logfile.write(stdout)
    if verbose:
        print stderr
    logfile.write(stderr)
    return stdout

def vcf_name_strip(vcffilename):
    basecols = vcffilename.split(".")
    vcfstrippedname = basecols[0]
    vcfstrippedname = vcfstrippedname + "-striptmp.vcf.gz"
    file_data = open(vcffilename, 'rb')
    outstrip = open(vcfstrippedname, 'wb')


    for file_line in file_data:
        cols = file_line.split('\t')
        if cols[0] == '#CHROM':  # Header line of VCF file

            if cols[8] == 'FORMAT':  # On header line, a FORMAT column next to the fixed columns?

                fixedgenos = cols[:9]
                orig_names = cols[9:]  # If so, remaining columns are the genotypes
                for orig_name in orig_names:

                    nodir = os.path.basename(orig_name)
                    dotsplits = nodir.split(".")
                    dotsplit = dotsplits[0]
                    nounder = dotsplit.replace("_", "-")
                    finalname = nounder[:15]

                    fixedgenos.append(finalname)
                outstrip.write('\t'.join(fixedgenos))
                outstrip.write("\n")
            else:
                print "Error. VCF file with no genotype. Cannot create sequence data."
                return
        else:
            outstrip.write(file_line)
    file_data.close()
    outstrip.close()
    shutil.move(vcfstrippedname, vcffilename)

def poplist_alter(poplistfile, namebase):
    poplist = {}
    poplistd = open(poplistfile, 'r')

    for poplistline in poplistd:
        popstrip = poplistline.rstrip()
        if len(popstrip) > 0:
            popcols = popstrip.split("\t")
            print popcols
            poplist[popcols[1]] = popcols[0]

    pedfilename = namebase + ".ped"
    pedfiled = open(pedfilename, 'r')
    pedtmp = open("pedtemp.ped", 'w')
    for pedfileline in pedfiled:
        pedstrip = pedfileline.rstrip()
        pedcols = pedstrip.split(" ")

        if pedcols[1] in poplist:
            pedcols[0] = poplist[pedcols[1]]
        newpedline = '\t'.join(pedcols)
        newpedline = newpedline + "\n"
        pedtmp.write(newpedline)
    pedfiled.close()
    pedtmp.close()
    shutil.move("pedtemp.ped", pedfilename)
    return pedfilename

def eigenstrat_convert(bcbase):
    eigenname = "par." + bcbase + "PED.EIGENSTRAT"
    eigennout = open(eigenname, 'w')
    eigennout.write("genotypename:    ")
    eigennout.write(bcbase)
    eigennout.write(".ped\n")
    eigennout.write("snpname:         ")
    eigennout.write(bcbase)
    eigennout.write(".map\n")
    eigennout.write("indivname:       ")
    eigennout.write(bcbase)
    eigennout.write(".ped\n")
    eigennout.write("genotypeoutname: ")
    eigennout.write(bcbase)
    eigennout.write(".geno\n")
    eigennout.write("snpoutname:      ")
    eigennout.write(bcbase)
    eigennout.write(".snp\n")
    eigennout.write("indivoutname:    ")
    eigennout.write(bcbase)
    eigennout.write(".ind\n")
    eigennout.write("familynames:     NO\n")
    eigennout.close()

    bash_command("convertf -p " + eigenname)