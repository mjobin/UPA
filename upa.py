#!/usr/bin/python

#####################################
#####        HPG Lab            #####
#####    updated Sept 2018      #####
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
import upa_input





if __name__ == "__main__":

    print "\n****************\nUPA\n****************\n"

    parser = argparse.ArgumentParser(description="# This script:\n"
                                                 "1. does stuff.\n\t"
                                                 "- ", formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('-bc_file', metavar='<bc_file>', help='Location of barcode files, Must have a newline at end.',
                        default="")
    parser.add_argument('-bc_leftspec', metavar='<bc_leftspec>', help='Extensions or name to the left of reference.',
                        default=".M.cf.")
    parser.add_argument('-bc_rightspec', metavar='<bc_rightspec>', help='Extensions or name to the right of quality.',
                        default=".s")
    parser.add_argument('-bam_list', metavar='<bam_list>', help='List of BAM files', default="")
    parser.add_argument('-vcf_file', metavar='<vcf_file>', help='User-processed VCF file.', default="")
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
    parser.add_argument('-ref', metavar='<ref>', help='Reference FASTA file.',
                        default='/data/genomes/onlynumber_nochr_MT_hg19.fa')
    parser.add_argument('-q', metavar='<q>', help='BWA min quality. 20 provides a fairly low cutoff',
                        default="20")
    parser.add_argument('-samindex', dest='samindex', help='Generate indexes for BAM files.',
                        action='store_true')
    parser.set_defaults(samindex=False)
    parser.add_argument('-diploid', dest='diploid', help='Diploid data.',
                        action='store_true')
    parser.set_defaults(diploid=True)
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
    parser.add_argument('-imptree', metavar='<imptree>',
                        help='Tree for Imputor.',
                        default="")
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
    parser.add_argument('-maxhops', metavar='<maxhops>', help='Number of hops to search in hops method.', default="5")
    parser.add_argument('-yhaplo', metavar='<yhaplo>',
                        help='Location of yhaplo. Leave blank to prevent it from running.',
                        default='')
    parser.add_argument('-poplistfile', metavar='<poplistfile>',
                        help='Text file where the FIRST column is the population and the Second coumn is the individual. Can be the same file as your plink keeplist',
                        default='')
    parser.add_argument('-lowk', metavar='<lowk>', help='Lowest K value for Admixture run',
                        default=2)
    parser.add_argument('-hik', metavar='<hik>', help='Highest K value for Admixture run',
                        default=10)
    parser.add_argument('-reps', metavar='<reps>', help='Number of replicates to be run at each K',
                        default=10)
    parser.add_argument('-tohaploid', dest='tohaploid', help='Haploid conversion.',
                        action='store_true')
    parser.set_defaults(tohaploid=False)
    parser.add_argument('-tvonly', dest='tvonly', help='Transversions only.',
                        action='store_true')
    parser.set_defaults(tvonly=False)
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
    parser.add_argument('-snprelatepca', dest='snprelatepca', help='snprelatepca.',
                        action='store_true')
    parser.set_defaults(snprelatepca=False)
    parser.add_argument('-smartpca', dest='smartpca', help='smartpca.',
                        action='store_true')
    parser.set_defaults(smartpca=False)
    parser.add_argument('-ancient', dest='ancient', help='Turn on aDNA switches.',
                        action='store_true')
    parser.set_defaults(ancient=False)
    parser.add_argument('-scriptsloc', metavar='<scriptsloc>', help='Location of external scripts.',
                        default='/data/scripts/')
    parser.add_argument('-binloc', metavar='<binloc>', help='Location of binary executables.',
                        default='/usr/local/bin/')
    parser.add_argument('-stripchr', dest='stripchr', help='Strip chr from your samples chromosome names.',
                        action='store_true')
    parser.set_defaults(stripchr=False)
    parser.add_argument('-addreadgroup', dest='addreadgroup', help='Add read group (RG) back to your sample BAMs.',
                        action='store_true')
    parser.set_defaults(addreadgroup=False)
    parser.add_argument('-callmethod', metavar='<callmethod>', help='Options: bcf, genocaller.',
                        default='bcf')
    parser.add_argument('-gcbedfile', metavar='<gcbedfile>', help='UCSC BED file for use with GenoCaller.',
                        default='')
    parser.add_argument('-gcindent', metavar='<gcindent>', help='Indent depth for use with GenoCaller.',
                        default='2')
    parser.add_argument('-plinkgeno', metavar='<plinkgeno>', help='Value for PLINK geno argument.',
                        default='0.99')




    # Parsing args
    args = parser.parse_args()
    wd = args.wd
    bcfile = args.bc_file
    bcleftspec = args.bc_leftspec
    bcrightspec = args.bc_rightspec
    bamlist = args.bam_list
    vcf_file = args.vcf_file
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
    yhaplo = args.yhaplo
    maxdepth = args.maxdepth
    maxheight = args.maxheight
    nsize = args.nsize
    msize = args.msize
    ncollect = args.ncollect
    maxhops = args.maxhops
    poplistfile = args.poplistfile
    imputor = bool(args.imputor)
    imptree = args.imptree
    admixture = bool(args.admixture)
    snprelatepca = bool(args.snprelatepca)
    smartpca = bool(args.smartpca)
    ancient = bool(args.ancient)
    scriptsloc = args.scriptsloc
    binloc = args.binloc
    stripchr = bool(args.stripchr)
    addreadgroup = bool(args.addreadgroup)
    maxgap = int(args.maxgap)
    mindepth = int(args.mindepth)
    haplogrepjava = bool(args.haplogrepjava)
    callmethod = args.callmethod
    gcbedfile = args.gcbedfile
    gcindent = args.gcindent
    plinkgeno = args.plinkgeno

    # adpipe
    lowk = int(args.lowk)
    hik = int(args.hik)
    reps = int(args.reps)
    tohaploid = bool(args.tohaploid)
    tvonly = bool(args.tvonly)
    termcrit = float(args.termcrit)
    optmethod = args.optmethod




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

    logfile.write("Arguments used:\n")
    logfile.write("__________________________________________:\n")
    for arg in vars(args):
        logfile.write(arg)
        logfile.write("\t")
        logfile.write(str(getattr(args, arg)))
        logfile.write("\n")

    if mito or ychr:
        diploid = False

    bcname = ""
    samplevcffile = ""
    bampreprocess = True


    print "\nChecking for input files..."
    if vcf_file:
        bampreprocess = False
        bcbase = os.path.basename(vcf_file)
        bcname, fileext = os.path.splitext(bcbase)
    elif bcfile != "" and bamlist == "":
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

    print "\nProcessing input files..."

    if bampreprocess:
        if samindex:
            print "\nIndexing..."
            bar = progressbar.ProgressBar()
            for i in bar(range(flength)):
                sample = flist[i]
                upa_util.bash_command("samtools index " + sample + ".bam", verbose, cmdfile, logfile)
        if stripchr:
            upa_input.stripchr(flist, verbose, cmdfile, logfile)
        if addreadgroup:
            upa_input.addreadgroup(flist, binloc, verbose, cmdfile, logfile)
        print "\nIndexing..."
        bar = progressbar.ProgressBar()
        for i in bar(range(flength)):
            sample = flist[i]
            upa_util.bash_command("samtools index " + sample + ".bam", verbose, cmdfile, logfile)


    if vcf_file:
        samplevcffile = vcf_file #User submitting a VCF file
    elif callmethod == 'bcf':
        samplevcffile = upa_input.bcfmpileup(flist, ref, bcname, regionrestrict, diploid, q, cmdfile, logfile)
    elif callmethod == 'genocaller':
        samplevcffile = upa_input.genocaller(flist, gcbedfile, bcname, gcindent, ref, regionrestrict, verbose, cmdfile, logfile)
    else:
        print "EROR: Unknown calling method " + callmethod
        exit(1)


    if samplevcffile == "":
        print "ERROR. Either specify a calling method for your BAM files or submit a pre-processed VCF file."
        exit(1)


    mergedvcfname = bcname + "-MERGED.vcf"
    mergedvcfgzipname = bcname + "-MERGED.vcf.gz"

    if vcfchromrename:
        renamefile = wd + "/" + vcfchromrename
        upa_util.bash_command("bcftools annotate --rename-chrs  " + renamefile + " " + samplevcffile + " -Oz -o " + bcname + "upatmp.vcf.gz", verbose, cmdfile, logfile)
        shutil.move(bcname + "upatmp.vcf.gz",  samplevcffile )
        upa_util.bash_command("bcftools index -f " + samplevcffile, verbose, cmdfile, logfile)


    print "Stripping long names from VCF genotypes. This should also match sample names to either the second column name of a barcode file or the first element of the file name on the a BAM file list."
    upa_util.vcf_name_strip(samplevcffile)

    vcfmergecmd = "bcftools merge -Ov -o " + mergedvcfname + " "
    if regionrestrict:
        vcfmergecmd = vcfmergecmd + "-r " + regionrestrict + " "

    if mergebamfile:
        print "\nMerging sample VCF file with external BAM reference" + mergebamfile   + "..."
        mergebambase = mergebamfile.rsplit(".", 1)[0]
        mergevcffile = mergebambase + ".vcf.gz"
        mpileupcmd = ("bcftools mpileup -q " + q + " -d 8000 -Ou -f " + ref + " " + mergebamfile)
        if regionrestrict:
            mpileupcmd = mpileupcmd + " -r " + regionrestrict
        mpileupcmd = mpileupcmd + " | bcftools call -Oz -m -o " + mergevcffile + " - "

        if diploid:
            pass
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

    #Default will not merge an external, but still region restrict
    vcfmergecmd = vcfmergecmd + samplevcffile

    if mergebamfile or mergevcffile:
        upa_util.bash_command(vcfmergecmd, True, cmdfile, logfile)
    else: # If still gzipped here, unzip
        if samplevcffile.endswith(".gz"):
            with gzip.open(samplevcffile, 'rb') as f_in, open(bcname + "-MERGED.vcf", 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)
        else:
            shutil.copy(samplevcffile, mergedvcfname)

    print "Stripping long names from VCF genotypes. This should also match sample names to either the second column name of a barcode file or the first element of the file name on the a BAM file list."
    upa_util.vcf_name_strip(mergedvcfname)

    upa_util.bash_command("bcftools index -f " + mergedvcfname, verbose, cmdfile, logfile)

    # IMPUTOR
    if imputor:
        if diploid:
            print "IMPUTOR works on haploid data only! Exiting."
            exit(1)
        print "Running IMPUTOR..."
        impcmd = "imputor.py -file " + mergedvcfname + " -out vcf -maxthreads " + threads + " -ncollect " + ncollect + " -maxheight " + maxheight + " -maxdepth " + maxdepth + " -passes 1 -msize " + msize + " -nsize " + nsize + " -maxhops " + maxhops
        if imptree:
            impcmd = impcmd + " -tree " + imptree
        upa_util.bash_command(impcmd, True, cmdfile, logfile)
        shutil.move(bcname + "-MERGED-out.vcf", mergedvcfgzipname)  # Overwrite with imputed sequence so pipeline knows which to use

    print "Converting " + mergedvcfname + " to PED format"
    upa_util.bash_command("plink --vcf " + mergedvcfname + " --double-id --allow-extra-chr --missing-phenotype 2 --recode 12 --out " + bcname, verbose, cmdfile, logfile)

    # Alter for pops
    if poplistfile:
        pedfilename = upa_util.poplist_alter(poplistfile, bcname)

    if snprelatepca:
        print "\nRunning SnpRelatePCA..."
        upa_util.bash_command("Rscript " + scriptsloc + "snprelatepca.R " + bcname + " " + threads, verbose, cmdfile, logfile)


    print "\nConvert to EIGENSTRAT format..."
    eigenparname = upa_util.eigenstrat_convert(bcname, verbose, cmdfile, logfile, diploid, ancient)

    if smartpca:
        smartpcacmd = "smartpca -p " + eigenparname
        upa_util.bash_command(smartpcacmd, verbose, cmdfile, logfile)

    if ychr:
        upa_util.bash_command(yhaplo + "/callHaplogroups.py -i " + mergedvcfname, verbose, cmdfile, logfile)

    if mito:
        print "\nCreating Region dictionary..."
        regdic = {}
        vcflist = []
        bar = progressbar.ProgressBar()
        for i in bar(range(flength)):
            sample = flist[i]
            stripname = upa_util.name_strip(sample)
            regdic[stripname] = upa_mito.gen_reg_line(sample+".bam", mindepth, maxgap, cmdfile, logfile)
        upa_mito.haplogrep_gen_hsd(flist, ref, bcname, regdic, cmdfile, logfile)

        print "Stripping long names from VCF genotypes. This should also match sample names to either the second column name of a barcode file or the first element of the file name on the a BAM file list."
        upa_util.vcf_name_strip(bcname + "-4hgrp.vcf")
        if haplogrepjava:
            upa_mito.haplogrep_java(bcname + "-4hgrp.vcf", regdic, scriptsloc, cmdfile, logfile)

    if admixture:
        print "Running Admixture..."
        upa_util.bash_command("plink --file " + bcname + " --make-bed --geno " + plinkgeno + " --allow-extra-chr --out " + bcname, verbose, cmdfile, logfile)
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
    print "Not configured as a module"
    exit(1)

