#!/usr/bin/python

#####################################
#####        HPG Lab            #####
#####    updated July 2018      #####
#####       MJJ                 #####
#####################################

# Driver and merging script for phylogenetics and imputation/correction
# Requires 'barcode' file with fastq prefix, output name and internal barcodes
# Author: Matthew Jobin, UCSC Human Paleogenomics Lab

import argparse
import os
import datetime
import progressbar
import gzip
import shutil
import upa_util
import upa_mito


if __name__ == "__main__":

    print "\n****************\nUPA\n****************\n"

    parser = argparse.ArgumentParser(description="# This script:\n"
                                                 "1. does stuff.\n\t"
                                                 "- ", formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('-bc_file', metavar='<bc_file>', help='location of barcode files, Must have a newline at end.',
                        default="")
    parser.add_argument('-bc_leftspec', metavar='<bc_leftspec>', help='extensions or name to the left of reference.',
                        default=".M.cf.")
    parser.add_argument('-bc_rightspec', metavar='<bc_rightspec>', help='extensions or name to the right of quality.',
                        default=".s")
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
    parser.add_argument('-vcfchromrename', metavar='<vcfchromrename>',
                        help='Use this is you are SURE you are merging into the same region just with different names!!',
                        default='')
    parser.add_argument('-mergevcffile', metavar='<mergevcffile>',
                        help='Larger VCF dataset to merge. MUST be indexed .bgzipped and variants ONLY!',
                        default='')
    parser.add_argument('-mergebamfile', metavar='<mergebamfile>',
                        help='Larger BAM dataset to merge.',
                        default='')
    parser.add_argument('-mito', dest='mito', help='Utilize upa_mito.py.',
                        action='store_true')
    parser.set_defaults(mito=False)
    parser.add_argument('-ychr', dest='ychr', help='Invoke Y chromosome functions.',
                        action='store_true')
    parser.set_defaults(ychr=False)
    parser.add_argument('-imputor', dest='imputor', help='Impute and tree construction',
                        action='store_true')
    parser.set_defaults(imputor=False)
    parser.add_argument('-vcfnamestrip', dest='vcfnamestrip', help='Strip long names from VCF entries. Can avoid some errors.',
                        action='store_true')
    parser.set_defaults(vcfnamestrip=False)
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
    parser.add_argument('-yhaplo', metavar='<yhaplo>',
                        help='Location of yhaplo. Leave blank to prevent it from running.',
                        default='')
    parser.add_argument('-poplistfile', metavar='<poplistfile>',
                        help='Text file where the FIRST column is the indiviodual and the THIRD coumn is the population. Can be the same file as your plink keeplist',
                        default='')
    parser.add_argument('-lowk', metavar='<lowk>', help='Lowest K value for Admixture run',
                        default=1)
    parser.add_argument('-hik', metavar='<hik>', help='Highest K value for Admixture run',
                        default=5)
    parser.add_argument('-reps', metavar='<reps>', help='Number of replicates to be run at each K',
                        default=10)
    parser.add_argument('-tohaploid', dest='tohaploid', help='Haploid conversion.',
                        action='store_true')
    parser.set_defaults(verbose=False)
    parser.add_argument('-tvonly', dest='tvonly', help='Transversions only.',
                        action='store_true')
    parser.set_defaults(verbose=False)
    parser.add_argument('-termcrit', metavar='<termcrit>', help='A termination criterion.',
                        default=0.0001)
    parser.add_argument('-optmethod', metavar='<optmethod>', help='Optimization method: em or block.',
                        default='block')
    parser.add_argument('-haplogrepjava', dest='haplogrepjava', help='haplogrepjava.',
                        action='store_true')
    parser.set_defaults(haplogrepjava=False)
    parser.add_argument('-maxgap', metavar='<maxgap>', help='Maximum gap in read before it is counted as a new region',
                        default=1)
    parser.add_argument('-mindepth', metavar='<mindepth>', help='Minimum depth to be in a region',
                        default=1)
    parser.add_argument('-admixture', dest='admixture', help='Run admixture.',
                        action='store_true')
    parser.set_defaults(admixture=False)


    # Parsing args
    args = parser.parse_args()
    wd = args.wd
    bcfile = args.bc_file
    bcleftspec = args.bc_leftspec
    bcrightspec = args.bc_rightspec
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
    mergebamfile = args.mergebamfile
    mito = bool(args.mito)
    ychr = bool(args.ychr)
    maxdepth = args.maxdepth
    maxheight = args.maxheight
    nsize = args.nsize
    msize = args.msize
    ncollect = args.ncollect
    maxhops = args.maxhops
    poplistfile = args.poplistfile
    imputor = bool(args.imputor)
    vcfnamestrip = bool(args.vcfnamestrip)
    admixture = bool(args.admixture)

    # adpipe
    lowk = int(args.lowk)
    hik = int(args.hik)
    reps = int(args.reps)
    tohaploid = bool(args.tohaploid)
    tvonly = bool(args.tvonly)
    termcrit = float(args.termcrit)
    optmethod = args.optmethod

    maxgap = int(args.maxgap)
    mindepth = int(args.mindepth)
    haplogrepjava = bool(args.haplogrepjava)

    # Setup
    os.chdir(wd)
    cwd = os.getcwd()
    print "Working in: ", cwd
    cmdfile = open("upa_cmds", 'w')
    today = datetime.date.today()
    logfilename = wd + "/out.upa." + str(today) + ".log"
    print "Logging to: ", logfilename

    logfile = open(logfilename, 'w')

    refbase = os.path.basename(ref)
    refname, fileext = os.path.splitext(refbase)

    flist = []

    # cmdfile.write("Arguments used:")
    # for arg in vars(args):
    #     cmdfile.write(arg)
    #     cmdfile.write("\t")
    #     # cmdfile.write(getattr(args,arg))
    #     cmdfile.write("\n")


    bcname = ""

    print "\nChecking for input files..."
    if bcfile != "" and bamlist == "":
        bcbase = os.path.basename(bcfile)
        bcname, fileext = os.path.splitext(bcbase)
        bcin = open(bcfile, 'r')
        for bcline in bcin:
            bccols = bcline.split()

            binfile = wd + "/" + bccols[1] + "/BWA_" + refname + "/" + bccols[1] + bcleftspec + refname + ".q" + q + bcrightspec
            if os.path.isfile(binfile + ".bam"):
                flist.append(binfile)
            else:
                print "ERROR: File " + binfile + " does not seem to have a corresponding BAM  file in " + wd
                exit(1)

    elif bcfile == "" and bamlist != "":
        bcbase = os.path.basename(bamlist)
        bcname, fileext = os.path.splitext(bamlist)
        bamin = open(bamlist, 'r')
        for bamline in bamin:
            binfile = bamline.rstrip()
            if os.path.isfile(binfile):
                flist.append(binfile.rsplit(".", 1)[0])
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
            upa_util.bash_command("samtools index " + sample + ".bam", verbose, cmdfile, logfile)

    # CREATE MERGED VCF
    print "\nCreating mpileup consensus and writing as a VCF..."
    vcflist = []
    bar = progressbar.ProgressBar()
    for i in bar(range(flength)):
        sample = flist[i]

        nodir = os.path.basename(sample)
        nodircols = nodir.split(".")
        finalsampname = nodircols[0]


        depthinfo = upa_util.bash_command("samtools depth " + sample + ".bam", verbose, cmdfile, logfile)

        mpileupcmd = ("bcftools mpileup -q " + q + " -d 8000 -Ou -f " + ref + " " + sample + ".bam")
        if regionrestrict:
            mpileupcmd = mpileupcmd + " -r " + regionrestrict
        mpileupcmd = mpileupcmd + " | bcftools call -Oz -m -o " + sample + ".vcf.gz - "

        if diploid:
            mpileupcmd = mpileupcmd + " --ploidy 2 "
        else:
            mpileupcmd = mpileupcmd + " --ploidy 1 "

        upa_util.bash_command(mpileupcmd, verbose, cmdfile, logfile)


        upa_util.bash_command("bcftools index -f " + sample + ".vcf.gz", verbose, cmdfile, logfile)

        if vcfchromrename:
            renamefile = wd + "/" + vcfchromrename
            upa_util.bash_command(
                "bcftools annotate --rename-chrs  " + renamefile + " " + sample + ".vcf.gz -Oz -o " + sample + ".a.vcf.gz", verbose, cmdfile, logfile)
            upa_util.bash_command("bcftools index -f " + sample + ".a.vcf.gz", verbose, cmdfile, logfile)
            vcflist.append(sample + ".a.vcf.gz")
        else:
            vcflist.append(sample + ".vcf.gz")


    print "\nMerging sample VCF files..."
    vcfmergecmd = "bcftools merge -Ov -o " + bcname + "-MERGED.vcf "
    if regionrestrict:
        vcfmergecmd = vcfmergecmd + "-r " + regionrestrict + " "
    if mergebamfile:
        mergebambase = mergebamfile.rsplit(".", 1)[0]
        mergevcffile = mergebambase + ".vcf.gz"
        depthinfo = upa_util.bash_command("samtools depth " + mergebamfile, verbose, cmdfile, logfile)
        mpileupcmd = ("bcftools mpileup -q " + q + " -d 8000 -Ou -f " + ref + " " + mergebamfile)
        if regionrestrict:
            mpileupcmd = mpileupcmd + " -r " + regionrestrict
        mpileupcmd = mpileupcmd + " | bcftools call -Oz -m -o " + mergevcffile + " - "

        if diploid:
            mpileupcmd = mpileupcmd + " --ploidy 2 "
        else:
            mpileupcmd = mpileupcmd + " --ploidy 1 "

        upa_util.bash_command(mpileupcmd, verbose, cmdfile, logfile)
        upa_util.bash_command("bcftools index -f " + mergevcffile, verbose, cmdfile, logfile)
        vcfmergecmd = vcfmergecmd + mergevcffile + " "
    elif mergevcffile:
        print "Attempting merge with " + mergevcffile
        print "\nWARNING: This will fail in a VERY ugly way if this file and all your files were not mapped using the "
        print "same reference sequence  AND all the chromosomes are named using the same convention!"
        vcfmergecmd = vcfmergecmd + mergevcffile + " "
    for vcff in vcflist:
        vcfmergecmd = vcfmergecmd + vcff + " "

    upa_util.bash_command(vcfmergecmd, True, cmdfile, logfile)


    print "Stripping long names from VCF genotypes. This should also match sample names to either the second coolumn name of a barcode file or the first element of the file name on the a BAM file list."
    upa_util.vcf_name_strip(bcname + "-MERGED.vcf")

    upa_util.bash_command("bcftools index -f " + bcname + "-MERGED.vcf", verbose, cmdfile, logfile)


    # IMPUTOR
    if imputor:
        if diploid:
            print "IMPUTOR works on haploid data only! Exiting."
            exit(1)
        print "Running IMPUTOR... should only be run for haploid data!"
        impcmd = "imputor.py -file " + bcname + "-MERGED.vcf -out vcf -maxthreads " + threads + " -ncollect " + ncollect + " -maxheight " + maxheight + " -maxdepth " + maxdepth + " -passes 1 -msize " + msize + " -nsize " + nsize + " -maxhops " + maxhops
        if imptree:
            impcmd = impcmd + " -tree " + imptree
        upa_util.bash_command_bare(impcmd, True, cmdfile, logfile)
        with open(bcname + "-MERGED-out.vcf", 'r') as f_in, gzip.open(bcname + "-E.vcf.gz", 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)
        shutil.move(bcname + "-E.vcf.gz",
                    bcname + "-MERGED.vcf.gz")  # Overwrite with imputed sequence so pipeline knows which to use

    print "Converting " + bcname + ".vcf to PED format"
    upa_util.bash_command(
        "plink --vcf " + bcname + "-MERGED.vcf --double-id --allow-extra-chr --missing-phenotype 2 --recode12 --out " + bcname, verbose, cmdfile, logfile)

    # Alter for pops
    if poplistfile:
        pedfilename = upa_util.poplist_alter(poplistfile, bcname)

    print "\nRunning SnpRelatePCA..."
    upa_util.bash_command("Rscript /data/scripts/snprelatepca.R " + bcname + " " + threads, verbose, cmdfile, logfile)

    # # Convert to EIGENSTRAT
    # upa_convert.eigenstrat_convert(bcbase)

    if ychr:
        upa_util.bash_command(yhaplo + "/callHaplogroups.py -i " + bcname + "-E-seqout.vcf", verbose, cmdfile, logfile)

    if mito:
        regdic = upa_mito.haplogrep_gen_hsd(flist, mindepth, maxgap, ref, bcname, cmdfile, logfile)
        if haplogrepjava:
            upa_mito.haplogrep_java(bcname + ".vcf", regdic, cmdfile, logfile)

    if admixture:
        print "Running Admixture..."
        upa_util.bash_command("plink --file " + bcname + " --make-bed --allow-extra-chr --out " + bcname, verbose, cmdfile, logfile)
        adpipeline = "adpipe.py -wd " + wd + " -file " + bcname
        if overwrite:
            adpipeline += " -overwrite"
        if verbose:
            adpipeline += " -verbose"
        adpipeline += " -threads " + str(threads) + " -lowk " + str(lowk) + " -hik " + str(hik) + " -reps " + str(
            reps) + " -termcrit " + str(termcrit) + " -optmethod " + str(optmethod)
        if tohaploid:
            adpipeline += " -tohaploid"
        if tvonly:
            adpipeline += " -tvonly"

        upa_util.bash_command(adpipeline, verbose, cmdfile, logfile)

    logfile.close()
    cmdfile.close()
    exit(0)

else:
    print "Not yet configured as a module"
    exit(1)

