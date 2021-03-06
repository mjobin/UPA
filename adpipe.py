#!/usr/bin/python

#####################################
#####        HPG Lab            #####
#####    updated July 2018      #####
#####       MJJ                 #####
#####################################

# Written and maintained by Matthew Jobin, UCSC Human Paleogenomics Lab
# Script for running Admixture etc on BAM files
# Requires 'barcode' file with fastq prefix, output name and internal barcodes

import argparse
from argparse import RawTextHelpFormatter
import os
import progressbar
import datetime
import shutil
import random
import time
import upa_util

def converttohaploid():
    """ Convert a tped file to haploid

    :return:
    """

    print "\nConverting to tped format..."
    upa_util.bash_command("plink --bed " + file + ".bed --bim " + file + ".bim --fam " + file + ".fam  --alleleACGT --recode transpose --out " + file, verbose, cmdfile, logfile)

    print "\nConverting to haploid..."
    tpedfile = open(file + ".tped", 'r')
    oldtped = []
    tlc = 0
    for tline in tpedfile:
        oldtped.append(tline)
        tlc = tlc + 1
    tpedfile.close()

    newtped = []

    bar = progressbar.ProgressBar()

    for i in bar(range(tlc)):
        tline = oldtped[i]

        cols = tline.rstrip().split()

        newtlist = []

        genotypes = cols[4:]
        nalleles = list(set(genotypes))
        alleles = []

        for g in nalleles:
            if g != '0':
                alleles.append(g)
        if tvonly:
            if 'A' in alleles and 'G' in alleles:
                continue
            if 'C' in alleles and 'T' in alleles:
                continue
        for icol in cols[0:4]:
            newtlist.append(icol)
        for j in xrange(0, len(genotypes), 2):
            thisg = random.choice([genotypes[j], genotypes[(j + 1)]])
            newtlist.append(thisg)
            newtlist.append(thisg)
        newtlist.append("\n")
        newtped.append(' '.join(newtlist))
    tpedfile.close()

    toutfile = open(file + ".h.tped", 'w')
    for toutline in newtped:
        toutfile.write(toutline)
    toutfile.close()

    shutil.copy(file + ".fam", file + ".h.tfam")

    print "\nConverting to bed format..."
    upa_util.bash_command("plink --tfile " + file + ".h --make-bed --out " + file + ".h", True, cmdfile, logfile)


def run_ad():
    """ Invokes ADMIXTURE software for the user-specified number of reps, between chosen low and high K.

    :return:
    """
    print "\nRunning Admixture..."

    famfilename = file + ".h.fam"
    ind2popname = filebase + "/" + filebase + ".ind2pop"
    famfile = open(famfilename, 'r')
    ind2pop = open(ind2popname, 'w')
    for famline in famfile:
        famcols = famline.split()
        ind2pop.write(famcols[0])
        ind2pop.write("\n")
    famfile.close()
    ind2pop.close()

    kcvss = []
    klogls = []
    kbests = {}

    pongfile = open(rdir + "pong_filemap", 'w')
    for k in xrange(lowk, (hik + 1)):
        print "\n*********** K:" + str(k)
        kreplist = []
        kcvs = []
        logls = []
        bar = progressbar.ProgressBar()
        for j in bar(xrange(0, reps + 1)):
            jreplist = []
            stdoutfilename = rdir + filebase + ".h." + str(k) + ".r" + str(j) + ".log"
            stdoutfile = open(stdoutfilename, 'w')
            stdoutfile.write(upa_util.bash_command("admixture --cv " + file + ".h.bed " + str(k) + " -j" + threads + " -s " + str(rng.getrandbits(32)) + " -C " + str(termcrit) + " -m " + optmethod, verbose, cmdfile, logfile))
            stdoutfile.close()
            pfile = filebase + ".h." + str(k) + ".P"
            qfile = filebase + ".h." + str(k) + ".Q"
            shutil.move(pfile, rdir + filebase + ".h." + str(k) + ".r" + str(j) + ".P")
            shutil.move(qfile, rdir + filebase + ".h." + str(k) + ".r" + str(j) + ".Q")
            grepcmd = "grep -h CV " + stdoutfilename
            grepline = upa_util.bash_command(grepcmd, True, cmdfile, logfile)
            grepcols = grepline.split()
            kcvs.append(float(grepcols[3]))
            jreplist.append(float(grepcols[3]))
            jreplist.append(filebase + ".h." + str(k) + ".r" + str(j))

            stdreadout = open(stdoutfilename, 'r')
            for stdoutline in stdreadout:
                if stdoutline.startswith("Loglikelihood:"):
                    loggrepcols = stdoutline.split()
                    logls.append(loggrepcols[1])
                    jreplist.append(loggrepcols[1])
                    continue
            jreplist.append(str(j))
            kreplist.append(jreplist)

            runid = "run"
            runid += str(j)
            runid += "_K"
            runid += str(k)
            pongfile.write(runid)
            pongfile.write("\t")
            pongfile.write(str(k))
            pongfile.write("\t")
            pongfile.write(filebase + ".h." + str(k) + ".r" + str(j) + ".Q")
            pongfile.write("\n")

        llbest = 0
        for i in xrange(reps):
            if logls[i] < logls[llbest]:
                llbest = i

        print ("For K = " + str(k) + " best index: " + str(llbest))
        pbestname = rdir + filebase + ".h." + str(k) + ".r" + str(llbest) + ".P"
        shutil.copy(pbestname, bestdir)
        qbestname = rdir + filebase + ".h." + str(k) + ".r" + str(llbest) + ".Q"
        shutil.copy(qbestname, bestdir)
        newqbestname = bestdir + "/" + filebase + "." + str(k) + ".r" + str(llbest) + ".Q"

        logbestname = rdir + filebase + ".h." + str(k) + ".r" + str(llbest) + ".log"
        shutil.move(logbestname, bestdir)

        kbests[k] = kreplist[llbest]

        kcvss.append(kcvs)
        klogls.append(logls)

    cvfileoutname = rdir + filebase + "-cvout.csv"
    cvfileout = open(cvfileoutname, 'w')

    for k in xrange(lowk, (hik + 1)):
        cvfileout.write(str(k))
        cvfileout.write(",")
    cvfileout.write("\n")
    for j in xrange(reps):
        for i in xrange(len(kcvss)):
            cvfileout.write(str(kcvss[i][j]))
            cvfileout.write(",")
        cvfileout.write("\n")
    cvfileout.close()
    upa_util.bash_command("Rscript cvsplot.R " + cvfileoutname, verbose, cmdfile, logfile)

    bestname = bestdir + "/bests.csv"
    bestout = open(bestname, 'w')
    bestout.write("K,CV,filebase,logL, Rep\n")
    for k, l in kbests.iteritems():
        bestout.write(str(k))
        bestout.write(",")
        for i in l:
            bestout.write(str(i))
            bestout.write(",")
        bestout.write("\n")
    bestout.close()
    pongfile.close()

    pongcmd = "pong -f -m " + rdir + "pong_filemap -i " + ind2popname
    print "\n\n\n\n Copy and paste the following to run pong:\n"
    print pongcmd
    print "\nYou will need a separate terminal logged in with ssh -X to visualize."
    print "On that terminal type :"
    print "\nfirefox &\n"
    print "And go to http://localhost:4000"
    print "\nAlternatively, you can simply download the folder " + filebase + "to your own machine"
    print "and run pong there. See http://brown.edu/Research/Ramachandran_Lab/projects/"

    elapsed_time = time.time() - start_time

    print "Number threads: " + str(threads) + " Elapsed time: " + str(elapsed_time)


if __name__ == "__main__":

    print "\n****************\nADPIPE\n****************\n"

    parser = argparse.ArgumentParser(description="# This script:\n"
                                                 "1. Prunes transverions.\n\t"
                                                 "2. Runs Admixture Reps times per each K.\n\t"
                                                 "3. Sets up visualization in PONG.\n\t"
                                                 "- ", formatter_class=RawTextHelpFormatter)

    parser.add_argument('-file', metavar='<file>', help='Name of .bed .bim .fam files WITH NO EXTENSION.',
                        required=True)
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

    args = parser.parse_args()
    wd = args.wd
    file = args.file
    verbose = bool(args.verbose)
    overwrite = bool(args.overwrite)
    threads = args.threads
    lowk = int(args.lowk)
    hik = int(args.hik)
    reps = int(args.reps)
    tohaploid = bool(args.tohaploid)
    tvonly = bool(args.tvonly)
    termcrit = float(args.termcrit)
    optmethod = args.optmethod

    os.chdir(wd)
    cwd = os.getcwd()
    print "Working in: ", cwd

    filecols = file.split(".")
    filebase = filecols[0]

    rdir = wd

    # try to create result folder, but exit if it is already there
    if os.path.exists(filebase):
        if overwrite:
            shutil.rmtree(filebase)
            os.mkdir(filebase)
            rdir = wd + "/" + filebase + "/"
        else:
            print "Error! Folder " + filebase + " already exists. Exiting."
            exit()
    else:
        os.mkdir(filebase)
        rdir = wd + "/" + filebase + "/"

    bestdir = rdir + "BEST"
    if os.path.exists(bestdir):
        print "Error! Folder " + bestdir + " already exists. Exiting."
        exit()
    else:
        os.mkdir(bestdir)

    cmdfile = open("ap_cmds", 'w')

    today = datetime.date.today()
    now = datetime.datetime.now()
    print now
    logfilename = wd + "/out.ap." + str(today) + ".log"
    print "Logging to: ", logfilename

    start_time = time.time()

    logfile = open(logfilename, 'w')

    rng = random.SystemRandom()  # Uses /dev/urandom

    print "\nChecking for input files..."

    if not os.path.isfile(file + ".bed"):
        print "ERROR: Cannot find " + wd + "/" + file + ".bed"
        exit()

    if not os.path.isfile(file + ".bim"):
        print "ERROR: Cannot find " + wd + "/" + file + ".bim"
        exit()

    if not os.path.isfile(file + ".fam"):
        print "ERROR: Cannot find " + wd + "/" + file + ".fam"
        exit()

    # Prune transversions
    if tohaploid:
        converttohaploid()
    else:
        shutil.copy(file + ".bed", file + ".h.bed")
        shutil.copy(file + ".bim", file + ".h.bim")
        shutil.copy(file + ".fam", file + ".h.fam")

    # Admixture
    run_ad()

    cmdfile.close()
    logfile.close()
