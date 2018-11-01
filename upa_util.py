#!/usr/bin/python

#####################################
#####        HPG Lab            #####
#####    updated Oct 2018       #####
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
import linecache
import re


def bash_command(cmd, verbose, cmdfile, logfile):
    """ Execute external program using bash

    :param cmd: Command line to be executed by bash.
    :param verbose: Verbose output to log.
    :param cmdfile: File storing external commands invoked.
    :param logfile: Output log.
    :return:
    """
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
    """ Strips long or extraneous file names. Useful for file conversions where name size limits apply.

    :param orig_name: Original file name.
    :return: Modified name.
    """
    nodir = os.path.basename(orig_name)
    dotsplits = nodir.split(".")
    dotsplit = dotsplits[0]
    nounder = dotsplit.replace("_", "-")
    finalname = nounder[:15]
    return finalname

def vcf_name_strip(vcffilename):
    """ Searches VCF file for lines containing smaple names and truncates them

    :param vcffilename: Original VCF file.
    :return: VCF file with shortened/altered names.
    """
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
    """ Use FAM-style flat file to search and restore population names

    :param poplistfile: File with population in first column and smaple name in second. FAM files work.
    :param namebase: Base name of PED file.
    :return: Name of PED file with altered population names.
    """
    poplist = {}
    poplistd = open(poplistfile, 'r')

    for poplistline in poplistd:
        popstrip = poplistline.rstrip()
        if len(popstrip) > 0:
            popcols = popstrip.split()
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
    """ Use convertf to convert from PED to EIGENSTRAT format


    :param bcbase: Base name of input file.
    :param verbose: Berbose output to log.
    :param cmdfile: File storing external commands invoked.
    :param logfile: Output log.
    :return: Name of Eigenstrat file.
    """
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
    """ Create parameter file and run SmartPCA


    :param bcbase: base name of input file
    :param diploid: are samples diploid?
    :param ancient: use lsqproject argument in SmartPCA
    :param verbose: verbose output to log
    :param cmdfile: file storing external commands invoked
    :param logfile: output log
    :return: Name of SmartPCA parameter file.
    """
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

def mergeref(refvcf, othervcf, diploid, mergefoundonly, annotate):
    """ Adds the read group information by using Picard

    :param refvcf: VCF file mapped to reference given by ref argument on input, normally the samples.
    :param othervcf: VCF file of external dataset.
    :param diploid: Are samples diploid? Are friends electric.
    :param mergefoundonly: Merged file will contain sites found in both file only.
    :param annotate: Annotate the ID column of the merged file from the external dataset (othervcf).
    :param verbose: Verbose output to log.
    :param cmdfile: File storing external commands invoked.
    :param logfile: Output log.
    :return: Name of merged VCF file.
    """
    #First read in the reference (normally, the sample) VCF, and create a line dictionary based on position


    mergevcf = refvcf[:-7]
    mergevcf += "-MERGED.vcf.gz"

    if refvcf[-3:] == ".gz":
        refun = refvcf[:-3]
        with bgzf.open(refvcf, 'rb') as f_in, open(refun, 'w') as f_out:
            shutil.copyfileobj(f_in, f_out)
        refvcf = refun

    if othervcf[-3:] == ".gz":
        otherun = othervcf[:-3]
        with bgzf.open(othervcf, 'rb') as f_in, open(otherun, 'w') as f_out:
            shutil.copyfileobj(f_in, f_out)
        othervcf = otherun


    print "\nReading " + refvcf +  "..."
    reffile = open(refvcf, 'r')
    ref_data = []
    for file_line in reffile:
        if len(file_line.rstrip()) > 0:  # Strip blank lines
            ref_data.append(file_line.rstrip())
    refheaderline = ""
    refheaderlist = []
    refdict = {}
    foundheader = False
    # bar = progressbar.ProgressBar()
    # for i in bar(range(len(ref_data))):
    for i in range(len(ref_data)):
        file_line = ref_data[i]
        cols = file_line.split()
        # print cols
        if foundheader: #from here on, its data
            # print cols[0]+"-"+cols[1] + " " + str(i)
            refdict[cols[0]+"-"+cols[1]] = i
        else: ##just add to header repository
            if cols[0] == '#CHROM':
                refheaderline = file_line
                refhdrcols = cols
                print " number of total columns in ref " + str(len(refhdrcols))
                foundheader = True
            elif "##fileformat" not in file_line:
                refheaderlist.append(file_line)
    reffile.close()


    foundheader = False
    #Next, read in
    print "\nReading " + othervcf +  "..."
    otherfile = open(othervcf, 'r')
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
            otherdict[cols[0]+"-"+cols[1]] = i
        else: ##just add to header repository
            if cols[0] == '#CHROM':
                otherheaderline = file_line
                othersamplenames = cols[9:]
                print " number of sample columns in other " + str(len(othersamplenames))
                foundheader = True
            elif "##fileformat" not in file_line:
                otherheaderlist.append(file_line)
    otherfile.close()

    oslen = len(othersamplenames)



    print "Writing to " + mergevcf
    mergeout = bgzf.BgzfWriter(mergevcf, 'wb')


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


    print "Merging...."
    bar = progressbar.ProgressBar()
    for key, lnum in bar(sorted(refdict.items(), key=refkeysort)):
    # for key, lnum in sorted(refdict.items(), key=refkeysort):
        foundother = False
        refline = linecache.getline(refvcf, lnum+1).strip() # Add one because linecache lines start on 1 not 0
        # print key + " " + str(lnum+1) + " " + refline
        refcols = refline.split('\t')
        if key in otherdict:
            foundother = True
            otnum = otherdict[key]
            otherline = linecache.getline(othervcf, otnum+1).strip()


            complist = []

            othertm = {}
            # print otherline
            othercols = otherline.split()

            # print "\n"
            #
            # print key + " " + str(lnum + 1) + " " + refcols[1] +  " Otherdict " + othercols[1]

            trueref = refcols[3]
            complist.append(trueref)
            truealts = refcols[4].split(",")
            for alt in truealts:
                complist.append(alt)

            # print "True REF " + trueref
            otherref = othercols[3]
            otheralts = othercols[4].split(",")

            if otherref in complist:
                pass
            else:
                complist.append(otherref)

            for k in range(len(otheralts)):
                if otheralts[k] in complist:
                    pass
                else:
                    complist.append(otheralts[k])



            # print complist

            otherrefloc = complist.index(otherref)
            othertm[0] = otherrefloc
            for k in range(len(otheralts)):
                othertm[k+1] = complist.index(otheralts[k])

            altlist = complist
            altlist.remove(trueref)


            # print "TM "
            # print othertm

            siteline = []
            for l in range (len(refcols)):
                if l == 4:
                    siteline.append(','.join(altlist))
                elif l == 2:
                    if annotate:
                        siteline.append(othercols[l])
                    else:
                        siteline.append(refcols[l])
                else:
                    siteline.append(refcols[l])

            #
            # print "final siteline"



            #construct
            for othersite in othercols[9:]:
                othersites = re.split("[/|]+", othersite)



                # print othersites
                olen = len(othersites)
                # print olen
                if olen > 1 and not diploid:
                    print "ERROR: not diploid but more than one site at " + key
                    exit(1)
                oconstruct = ""
                for i in xrange(olen):
                    osite = othersites[i]
                    if osite == ".":
                        oconstruct += "."
                        # print osite + " becomes ."
                    else:
                        # print osite + " becomes " + str(othertm[int(osite)])
                        oconstruct += str(othertm[int(osite)])
                    if i < olen-1:
                        oconstruct += "/" # FIXME this always ouputs the unphased marker



                siteline.append(oconstruct)
        else:
            # print key + " " + str(lnum+1) + " no match"
            if mergefoundonly:
                siteline = ""
            else:
                refline = linecache.getline(refvcf, lnum+1).strip()
                refcols = refline.split('\t')
                siteline = refcols
                for nom in range(oslen):
                    if diploid:
                        siteline.append("./.") # FIXME this always ouputs the unphased marker
                    else:
                        siteline.append(".")


        ##Now check if its all missing or empty
        allmissing = True
        for i in xrange(9, len(siteline)):
            site = siteline[i]
            if site != "./." and site != "." and site != ".|.":
                allmissing = False
        if allmissing:
            # print "At " + key + " all sites missing, skipping."
            pass
        else:
            siteout = '\t'.join([str(x) for x in siteline])
            # print siteout
            siteout += "\n"
            if mergefoundonly:
                if foundother:
                    if len(siteline) != len(outhdr):
                        print "ERROR: Line in merged VCF has " + str(len(siteline)) + " but header line has " + str(
                            len(outhdr))
                    mergeout.write(siteout)
            else:
                if len(siteline) != len(outhdr):
                    print "ERROR: Line in merged VCF has " + str(len(siteline)) + " but header line has " + str(
                        len(outhdr))
                mergeout.write(siteout)
    mergeout.close()
    return mergevcf



def refkeysort(string):
    """ Sort function for compounded chromosome-position keys in dicitonary

    :param string: Input string. Should be in format <chromosome>-<position>.
    :return: Numerical chromosome-first sort number.
    """
    thekey = string[0]
    halves = thekey.split("-")
    compound = 0
    if str(halves[0]).isdigit():
        compound += int(halves[0]) * 1000000000000  # Longer than any single chromosome
    else: #for X, Y MT
        compound += ord(halves[0][0]) * 1000000000000
    compound += int(halves[1])
    return compound
