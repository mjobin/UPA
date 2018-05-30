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

def bash_command(cmd, bverbose):
    cmdfile.write(cmd)
    cmdfile.write("\n\n")
    subp = subprocess.Popen(['/bin/bash', '-c', cmd], stdout=PIPE, stderr=PIPE)
    subp.wait()
    theout = subp.stdout.read()
    if bverbose:
        print theout
    logfile.write(theout)
    theerr = subp.stderr.read()
    if bverbose:
        print theerr
    logfile.write(theerr)
    return theout

def bash_command_bare(cmd, bverbose):
    cmdfile.write(cmd)
    cmdfile.write("\n\n")
    subp = subprocess.Popen(['/bin/bash', '-c', cmd])
    subp.wait()

def vcf_name_strip(vcffilename):
    basecols = vcffilename.split(".")
    vcfstrippedname = basecols[0]
    vcfstrippedname = vcfstrippedname + "-striptmp.vcf.gz"
    file_data = gzip.open(vcffilename, 'rb')
    outstrip = gzip.open(vcfstrippedname, 'wb')

    for file_line in file_data:
        cols = file_line.split('\t')
        if cols[0] == '#CHROM':  # Header line of VCF file
            if cols[8] == 'FORMAT':  # On header line, a FORMAT column next to the fixed columns?
                fixedgenos = cols[:9]
                orig_names = cols[9:]  # If so, remaining columns are the genotypes
                for orig_name in orig_names:
                    dotsplits = orig_name.split(".")
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
    parser.add_argument('-threads', metavar='<threads>',
                        help='The number of threads to assign to each task when possible',
                        default="23")
    parser.add_argument('-ref', metavar='<ref>', help='',
                        default='/data/genomes/hg19.fa')
    parser.add_argument('-q', metavar='<q>', help='BWA min quality. 20 provides a fairly low cutoff',
                        default="20")
    parser.add_argument('-samindex', dest='samindex', help='Generate indexes for BAM files.',
                        action='store_true')
    parser.set_defaults(samindex=False)
    parser.add_argument('-diploid', dest='diploid', help='Diploid data.',
                        action='store_true')
    parser.set_defaults(diploid=False)
    parser.add_argument('-regionrestrict', metavar='<regionrestrict>', help='Restrict to a region.',
                        default='')
    parser.add_argument('-vcfchromrename', metavar='<vcfchromrename>', help='Use this is you are SURE you are merging into the same region just with different names!!',
                        default='')
    parser.add_argument('-mergevcffile', metavar='<mergevcffile>', help='Larger VCF dataset to merge. MUST be indexed .bgzipped and variants ONLY!',
                        default='')
    parser.add_argument('-mommy', dest='mommy', help='Invoke mommy.py.',
                        action='store_true')
    parser.set_defaults(mommy=False)
    parser.add_argument('-daddy', dest='daddy', help='Invoke daddy.py.',
                        action='store_true')
    parser.set_defaults(daddy=False)
    parser.add_argument('-maxheight', metavar='<maxheight>',
                        help='Height of search toward root for collecting neighbors.',
                        default="3")
    parser.add_argument('-maxdepth', metavar='<maxdepth>',
                        help='Depth of search down descendent nodes for collecting neighbors.', default="3")
    parser.add_argument('-nsize', metavar='<nsize>',
                        help='Number of neighbors that must be identical in order to impute.',
                        default="3")
    parser.add_argument('-msize', metavar='<msize>',
                        help='Number of neighbors that must be identical in order to impute for missing data.',
                        default="2")
    parser.add_argument('-ncollect', metavar='<ncollect>', help='rootward, hops, distance, mono',
                        default='rootward')
    parser.add_argument('-maxhops', metavar='<maxhops>', help='Number of alternate runs.', default="5")



    #Parsing args
    args = parser.parse_args()
    wd = args.wd
    bcfile = args.bc_file
    bamlist = args.bam_list
    verbose = bool(args.verbose)
    overwrite = bool(args.overwrite)
    threads = args.threads
    ref = args.ref
    q = args.q
    samindex = bool(args.samindex)
    diploid = bool(args.diploid)
    regionrestrict = args.regionrestrict
    vcfchromrename = args.vcfchromrename
    mergevcffile = args.mergevcffile
    mommy = bool(args.mommy)
    daddy = bool(args.daddy)
    maxdepth = args.maxdepth
    maxheight = args.maxheight
    nsize = args.nsize
    msize = args.msize
    ncollect = args.ncollect
    maxhops = args.maxhops

    #Setup
    os.chdir(wd)
    cwd = os.getcwd()
    print "Working in: ", cwd
    cmdfile = open("tube_cmds", 'w')
    today = datetime.date.today()
    logfilename = wd + "/out.tubeamp." + str(today) + ".log"
    print "Logging to: ", logfilename

    logfile = open(logfilename, 'w')


    refbase = os.path.basename(ref)
    refname, fileext = os.path.splitext(refbase)

    flist = []

    bcname = ""

    print "\nChecking for input files..."
    if bcfile != "" and bamlist == "":
        print "bcfile " + bcfile
        bcbase = os.path.basename(bcfile)
        bcname, fileext = os.path.splitext(bcbase)
        bcin = open(bcfile, 'r')
        for bcline in bcin:
            print bcline
            bccols = bcline.split("\t")
            binfile = bccols[1] + "/BWA_" + refname + "/" + bccols[1] + ".M.cf." + refname + ".q" + q + ".s"
            if os.path.isfile(binfile + ".bam"):
                flist.append(binfile)
            else:
                print "ERROR: Barcode file entry " + binfile + " does not seem to have a corresponding BAM file."
                exit(1)
    elif bcfile == "" and bamlist != "":
        print "bamlist" + bamlist
        bcbase = os.path.basename(bamlist)
        bcname, fileext = os.path.splitext(bamlist)
        bamin = open(bamlist, 'r')
        for bamline in bamin:
            binfile = bamline.rstrip()
            if os.path.isfile(binfile):
                flist.append(binfile)
            else:
                print "ERROR: File " + binfile + " does not seem to have a corresponding BAM  file in " + wd
                exit(1)
    else:
        print "Either the bc_file or the bam_list args should be used, but not both."
        exit(1)

    flength = len(flist)
    print "Number of entries: ", flength

    if samindex:
        print "\nIndexing..."
        bar = progressbar.ProgressBar()
        for i in bar(range(flength)):
            sample = flist[i]
            bash_command("samtools index " + sample + ".bam", verbose)

    #CREATE MERGED VCF
    print "\nCreating mpileup consensus and writing as a VCF..."
    vcflist = []
    bar = progressbar.ProgressBar()
    for i in bar(range(flength)):
        sample = flist[i]

        depthinfo = bash_command("samtools depth " + sample + ".bam", False)

        mpileupcmd = ("bcftools mpileup -q " + q + " -d 8000 -Ou -f " + ref + " " + sample + ".bam")
        if regionrestrict:
            mpileupcmd = mpileupcmd + " -r " + regionrestrict
        mpileupcmd = mpileupcmd + " | bcftools call -Oz -m -o " + sample + ".vcf.gz - "
        if diploid:
            mpileupcmd = mpileupcmd + " --ploidy 2 "
        else:
            mpileupcmd = mpileupcmd + " --ploidy 1 "

        bash_command(mpileupcmd, verbose)
        bash_command("bcftools index -f " + sample + ".vcf.gz", verbose)

        if vcfchromrename:
            renamefile = wd + "/" + vcfchromrename
            bash_command("bcftools annotate --rename-chrs  " + renamefile + " " + sample + ".vcf.gz -Oz -o " + sample + ".a.vcf.gz", verbose)
            bash_command("bcftools index -f " + sample + ".a.vcf.gz", verbose)
            vcflist.append(sample + ".a.vcf.gz")
        else:
            vcflist.append(sample + ".vcf.gz")

    print "\nMerging sample VCF files..."
    vcfmergecmd = "bcftools merge -Oz -o " + bcname + "-MERGED.vcf.gz "
    if regionrestrict:
        vcfmergecmd = vcfmergecmd + "-r " + regionrestrict + " "
    if mergevcffile:
        print "Attempting merge with " + mergevcffile
        print "\nWARNING: This will fail in a VERY ugly way if this file and all your files were not mapped using the "
        print "same reference sequence  AND all the chromosomes are named using the same converntion!"
        vcfmergecmd = vcfmergecmd + mergevcffile + " "
    for vcff in vcflist:
        vcfmergecmd = vcfmergecmd + vcff + " "

    bash_command(vcfmergecmd, True)


    print "Stripping long names from VCF genotypes."
    vcf_name_strip(bcname + "-MERGED.vcf.gz")


    #IMPUTOR
    if imputor:
        if diploid:
            print "IMPUTOR works on haploid data only! Exiting."
            exit(1)
        print "Running IMPUTOR... should only be run for haploid data!"
        impcmd = "imputor.py -file " + bcname + "-MERGED.vcf.gz -out vcf -maxthreads " + threads + " -ncollect " + ncollect + " -maxheight " + maxheight + " -maxdepth " + maxdepth + " -passes 1 -msize " + msize + " -nsize " + nsize + " -maxhops " + maxhops
        if imptree:
            impcmd = impcmd + " -tree " + imptree
        bash_command_bare(impcmd, True)
        with open(bcname + "-MERGED-out.vcf", 'r') as f_in, gzip.open(bcname + "-E.vcf.gz", 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)
    else:
        shutil.move(bcname + "-MERGED.vcf.gz",
                    bcname + "-E.vcf.gz")  # Overwrite with imputed sequence so pipeline knows which to use
    ename = bcname + "-E"



    logfile.close()
    cmdfile.close()
    exit(0)

else:
    print "Not yet configured as a module"
    exit(1)