#!/usr/bin/python

#####################################
#####        HPG Lab            #####
#####    updated July 2018      #####
#####       MJJ                 #####
#####################################

# Module for processing mitochondrial data from UPA
# Author: Matthew Jobin, UCSC Human Paleogenomics Lab

import os
import progressbar
import upa_util



def gen_reg_line(sample, mindepth, maxgap, cmdfile, logfile):
    regline = ""
    depthinfo = upa_util.bash_command("samtools depth " + sample, False, cmdfile, logfile)

    depthlines = depthinfo.split("\n")
    curstart = 0
    lastpos = 0
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
            if curpos <= lastpos + maxgap:  # Continuous run
                if curdepth >= mindepth:
                    pass
                else:  # depth too small, stop here
                    curend = lastpos
                    regline = regline + str(curstart) + "-" + str(curend) + ";"
                    curstart = curpos
            else:  # discontinous
                curend = lastpos
                regline = regline + str(curstart) + "-" + str(curend) + ";"
                curstart = curpos
        lastpos = curpos
    return regline


def haplogrep_gen_hsd(flist, ref, bcname, regdic, cmdfile, logfile):
    flength = len(flist)
    fbamlist = []

    print "\nCreating mpileup consensus and writing as a VCF..."
    bar = progressbar.ProgressBar()
    for i in bar(range(flength)):
        sample = flist[i]
        fbamlist.append(sample + ".bam")

    mtmcmd = "bcftools mpileup -I -d 8000 -Ov -f " + ref + " "
    for fbam in fbamlist:
        mtmcmd = mtmcmd + fbam + " "

    mtmcmd = mtmcmd + "| bcftools call -V indels --ploidy 1 -Ov -m -v -o " + bcname + "-4hgrp.vcf"

    upa_util.bash_command(mtmcmd, False, cmdfile, logfile)
    return regdic

def haplogrep_java(invcf, regdic, scriptsloc, cmdfile, logfile):
    filebase, filext = os.path.splitext(invcf)

    upa_util.bash_command("java -jar " + scriptsloc + "haplogrep-2.1.1.jar --format vcf --in " + invcf + " --out " + filebase + ".hsd --phylotree 17", False, cmdfile, logfile)


    firstfile = filebase + "-FIRST.hsd"
    # upa_util.bash_command("java -jar " + scriptsloc + "haplogrep-2.1.1.jar --format vcf --in " + invcf + " --out " + firstfile + " --phylotree 17", False, cmdfile, logfile)

    # # then edit that HSD file and do it again
    # hsdoutlines = []
    # hsdfirstfile = open(firstfile, 'r')
    # for hsdline in hsdfirstfile:
    #
    #     hsdnewline = ""
    #     hsdcols = hsdline.strip().split("\t")
    #     if hsdcols[0] == "SampleID":
    #         hsdheadline = hsdline.strip()
    #     else:
    #         if len(regdic[hsdcols[0]]) > 1:  # LEAVE OUT indivs for whom where are no viable regions
    #
    #
    #
    #             for i in range(len(hsdcols)):
    #                 if i == 0 or i == 5:
    #
    #                     variantcols = hsdcols[i].split(";")
    #                     for variant in variantcols:
    #                         varvar = variant.strip()
    #                         print varvar
    #
    #                     hsdnewline = hsdnewline + hsdcols[i] + "\t"
    #                 elif i == 1:
    #                     hsdnewline = hsdnewline + regdic[hsdcols[0]] + "\t"
    #                 elif i == 2:
    #                     hsdnewline = hsdnewline + "?\t"
    #                 else:
    #                     pass
    #             hsdoutlines.append(hsdnewline)
    # hsdfirstfile.close()
    #
    # secondfile = filebase + "-SECOND.hsd"
    #
    # hsdsecondfile = open(secondfile, 'w')
    # hsdsecondfile.write("SampleID\tRange\tHaplogroup\tPolymorphisms")
    # hsdsecondfile.write("\n")
    # for hsdoutline in hsdoutlines:
    #     hsdoutcols = hsdoutline.split("\t")
    #     if len(hsdoutcols[3]) > 0:  # LEAVE OUT indivs for whom where are no viable regions
    #         hsdsecondfile.write(hsdoutline)
    #         hsdsecondfile.write("\n")
    # hsdsecondfile.close()
    #
    # finalfile = filebase + "-FINAL.hsd"
    #
    # upa_util.bash_command("java -jar /data/scripts/haplogrep-2.1.1.jar --format hsd --in " + secondfile + " --out " + finalfile + " --phylotree 17", False, cmdfile, logfile)