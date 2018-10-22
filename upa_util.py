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
import progressbar
from Bio import bgzf

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

def name_strip(orig_name):
    nodir = os.path.basename(orig_name)
    dotsplits = nodir.split(".")
    dotsplit = dotsplits[0]
    nounder = dotsplit.replace("_", "-")
    finalname = nounder[:15]
    return finalname

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
                    fixedgenos.append(name_strip(orig_name))
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

def eigenstrat_convert(bcbase, verbose, cmdfile, logfile, diploid, ancient):
    eigenname = "par." + bcbase + ".PED.EIGENSTRAT"
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
    eigennout.write("familynames:    YES\n")
    eigennout.close()

    bash_command("convertf -p " + eigenname, verbose, cmdfile, logfile)
    return eigenname



def eigenstrat_smartpca(bcbase, diploid, ancient, verbose, cmdfile, logfile):
    eigenname = "par.smartpca" + bcbase
    eigennout = open(eigenname, 'w')
    eigennout.write("genotypename:    ")
    eigennout.write(bcbase)
    eigennout.write(".geno\n")
    eigennout.write("snpname:         ")
    eigennout.write(bcbase)
    eigennout.write(".snp\n")
    eigennout.write("indivname:       ")
    eigennout.write(bcbase)
    eigennout.write(".ind\n")
    eigennout.write("evecoutname:    ")
    eigennout.write(bcbase)
    eigennout.write(".evec\n")
    eigennout.write("evaloutname:    ")
    eigennout.write(bcbase)
    eigennout.write(".eval\n")
    eigennout.write("familynames:     NO\n")
    if diploid:
        pass
    else:
        eigennout.write("fstonly:     YES\n")
    if ancient:
        eigennout.write("lsqproject:     YES\n")
    eigennout.close()

    bash_command("smartpca -p " + eigenname, verbose, cmdfile, logfile)
    return eigenname

def createrefallelelist(ref, other):

    #Create a dictionary based on position
    reflines = ref.split("\n")

    aflines = []

    refdict = {}
    for refline in reflines:
        refcols = refline.split()
        if len(refcols) > 4:
            refdict[refcols[1]] = refcols

    otherlines = other.split("\n")
    otherdict = {}
    for otherline in otherlines:
        othercols = otherline.split()
        if len(othercols) > 4:
            otherdict[othercols[1]] = othercols

    for key in refdict:
        if key in otherdict:
            print "At position " + key + " other ID is " + otherdict[key][2] + " and its REF is " + otherdict[key][3] + " while sample REF is " + refdict[key][3] + "\n"
            aflines.append(otherdict[key][2] + " " + refdict[key][3] + "\n")

    return aflines



def mergeref(refvcf, othervcf, diploid):
    #First read in the reference (normally, the sample) VCF, and create a line dictionary based on position

    mergevcf = refvcf[:-7]
    mergevcf += "-MERGED.vcf"


    print "\nReading " + refvcf +  "..."
    reffile = bgzf.BgzfReader(refvcf, 'r')
    ref_data = []
    for file_line in reffile:
        if len(file_line.rstrip()) > 0:  # Strip blank lines
            ref_data.append(file_line.rstrip())
    refheaderline = ""
    refheaderlist = []
    refdict = {}
    foundheader = False
    bar = progressbar.ProgressBar()
    for i in bar(range(len(ref_data))):
        file_line = ref_data[i]
        cols = file_line.split('\t')
        if foundheader: #from here on, its data
            refdict[cols[1]] = cols
        else: ##just add to header repository
            if cols[0] == '#CHROM':
                refheaderline = file_line
                refhdrcols = cols
                foundheader = True
            elif "##fileformat" not in file_line:
                refheaderlist.append(file_line)
    reffile.close()

    foundheader = False
    #Next, read in
    print "\nReading " + othervcf +  "..."
    otherfile = bgzf.BgzfReader(othervcf, 'r')
    other_data = []
    for file_line in otherfile:
        if len(file_line.rstrip()) > 0:  # Strip blank lines
            other_data.append(file_line.rstrip())
    otherheaderline = ""
    otherheaderlist = []
    otherdict = {}
    foundheader = False
    othersamplenames = []
    bar = progressbar.ProgressBar()
    for i in bar(range(len(other_data))):
        file_line = other_data[i]
        cols = file_line.split('\t')
        if foundheader: #from here on, its data
            otherdict[cols[1]] = cols
        else: ##just add to header repository
            if cols[0] == '#CHROM':
                otherheaderline = file_line
                othersamplenames = cols[9:]
                foundheader = True
            elif "##fileformat" not in file_line:
                otherheaderlist.append(file_line)
    otherfile.close()

    oslen = len(othersamplenames)

    print "Writing to " + mergevcf
    mergeout = open(mergevcf, 'w')


    #Merged header
    mergeout.write("##fileformat=VCFv4.2\n")
    mergeout.write("##UPA merged file headers for " + refvcf + "\n")
    for refhdrline in refheaderlist:
        mergeout.write(refhdrline)
        mergeout.write("\n")
    mergeout.write("##UPA merged file headers for " + othervcf + "\n")
    for otherhdrline in otherheaderlist:
        mergeout.write(otherhdrline)
        mergeout.write("\n")
    mergeout.write("##UPA merged " + refvcf + " and " + othervcf + " with REF alleles set to those of " + refvcf + " and all-missing sites ignored.\n")


    outhdr = refhdrcols
    for osn in othersamplenames:
        outhdr.append(osn)
    outhdrlen = len(outhdr)
    print "Header has " + str(outhdrlen) + " columns."
    hdrline = '\t'.join(outhdr)
    mergeout.write(hdrline)
    mergeout.write("\n")

    for key in refdict:
        siteline = refdict[key]
        if key in otherdict:
            print "Match at position " + key
            siteline.append(otherdict[key][9:])
        else:
            print "No match at position " + key
            for nom in range(oslen):
                if diploid:
                    siteline.append("./.")
                else:
                    siteline.append(".")
        if len(siteline) != len(outhdr):
            print "ERROR: Line in merged VCF has " + str(len(siteline)) + " but header line has " + str(len(outhdr))

        ##Now check if its all missing
        allmissing = True
        for i in xrange(9, len(siteline)):
            site = siteline[i]
            if site != "./." and site != ".":
                allmissing = False
        if allmissing:
            print "At " + key + " all sites missing, skipping."
        else:
            siteout = '\t'.join([str(x) for x in siteline])
            siteout += "\n"
            mergeout.write(siteout)
    mergeout.close()
    return mergevcf

