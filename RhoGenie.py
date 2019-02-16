#!/usr/bin/env python3
from collections import defaultdict
import re
import os
import textwrap
import argparse
import sys

# TODO: ADD BITSCORE CUTOFF FOR EACH HMM HIT


def deExt(string):
    ls = string.split(".")
    ls2 = ls[0:len(ls)-1]
    outstring = "".join(ls2)
    return outstring


def unique(ls, ls2):
    unqlist = []
    for i in ls:
        if i not in unqlist and i in ls2:
            unqlist.append(i)
    return len(unqlist)


def derep(ls):
    outLS = []
    for i in ls:
        if i not in outLS:
            outLS.append(i)
    return outLS


def cluster(data, maxgap):
    '''Arrange data into groups where successive elements
       differ by no more than *maxgap*

        #->>> cluster([1, 6, 9, 100, 102, 105, 109, 134, 139], maxgap=10)
        [[1, 6, 9], [100, 102, 105, 109], [134, 139]]

        #->>> cluster([1, 6, 9, 99, 100, 102, 105, 134, 139, 141], maxgap=10)
        [[1, 6, 9], [99, 100, 102, 105], [134, 139, 141]]

    '''
    # data = sorted(data)
    data.sort(key=int)
    groups = [[data[0]]]
    for x in data[1:]:
        if abs(x - groups[-1][-1]) <= maxgap:
            groups[-1].append(x)
        else:
            groups.append([x])
    return groups


def lastItem(ls):
    x = ''
    for i in ls:
        x = i
    return x


def RemoveDuplicates(ls):
    empLS = []
    for i in ls:
        if i not in empLS:
            empLS.append(i)
        else:
            pass
    return empLS


def allButTheLast(iterable, delim):
    x = ''
    length = len(iterable.split(delim))
    for i in range(0, length-1):
        x += iterable.split(delim)[i]
        x += delim
    return x[0:len(x)-1]


def secondToLastItem(ls):
    x = ''
    for i in ls[0:len(ls)-1]:
        x = i
    return x


def pull(item, one, two):
    ls = []
    counter = 0
    for i in item:
        if counter == 0:
            if i != one:
                pass
            else:
                counter += 1
                ls.append(i)
        else:
            if i != two:
                ls.append(i)
            else:
                ls.append(i)
                counter = 0
    outstr = "".join(ls)
    return outstr


def stabilityCounter(int):
    if len(str(int)) == 1:
        string = (str(0) + str(0) + str(0) + str(0) + str(int))
        return (string)
    if len(str(int)) == 2:
        string = (str(0) + str(0) + str(0) + str(int))
        return (string)
    if len(str(int)) == 3:
        string = (str(0) + str(0) + str(int))
        return (string)
    if len(str(int)) == 4:
        string = (str(0) + str(int))
        return (string)


def replace(stringOrlist, list, item):
    emptyList = []
    for i in stringOrlist:
        if i not in list:
            emptyList.append(i)
        else:
            emptyList.append(item)
    outString = "".join(emptyList)
    return outString


def remove(stringOrlist, list):
    emptyList = []
    for i in stringOrlist:
        if i not in list:
            emptyList.append(i)
        else:
            pass
    outString = "".join(emptyList)
    return outString


def remove2(stringOrlist, list):
    emptyList = []
    for i in stringOrlist:
        if i not in list:
            emptyList.append(i)
        else:
            pass
    # outString = "".join(emptyList)
    return emptyList


def removeLS(stringOrlist, list):
    emptyList = []
    for i in stringOrlist:
        if i not in list:
            emptyList.append(i)
        else:
            pass
    return emptyList


def fasta(fasta_file):
    seq = ''
    header = ''
    Dict = defaultdict(lambda: defaultdict(lambda: 'EMPTY'))
    for i in fasta_file:
        i = i.rstrip()
        if re.match(r'^>', i):
            if len(seq) > 0:
                Dict[header] = seq
                header = i[1:]
                header = header.split(" # ")[0]
                seq = ''
            else:
                header = i[1:]
                header = header.split(" # ")[0]
                seq = ''
        else:
            seq += i
    Dict[header] = seq
    # print(count)
    return Dict


def filter(list, items):
    outLS = []
    for i in list:
        if i not in items:
            outLS.append(i)
    return outLS


def delim(line):
    ls = []
    string = ''
    for i in line:
        if i != " ":
            string += i
        else:
            ls.append(string)
            string = ''
    ls = filter(ls, [""])
    return ls


parser = argparse.ArgumentParser(
    prog="RhoGenie.py",
    formatter_class=argparse.RawDescriptionHelpFormatter,
    description=textwrap.dedent('''
    *******************************************************

    Developed by Arkadiy Garber and Gustavo Ramírez;
    University of Southern California, Earth Sciences
    Please send comments and inquiries to arkadiyg@usc.edu

    *******************************************************
    '''))

parser.add_argument('-bin_dir', type=str, help="directory of ORFs in FASTA format")
parser.add_argument('-bin_ext', type=str, help="extension for files with ORFs")
parser.add_argument('-outdir', type=str, help="output directory (will be created if does not exist)", default="genie_out")
parser.add_argument('-out', type=str, help="basename of output file (default = out)", default="out")
parser.add_argument('--makeplots', type=str, help="Would you like GeoGenie to make some figures from your data? y = yes, n = no (default = n). "
                                                 "If so, you will need to have Rscipt installed. It is a way for R to be called directly from the command line. "
                                                 "Warning: this part of the program is currently under beta-testing, and if there are any problems running Rscript, "
                                                 "or installing any of the required packages, you may get a bunch of error messages at the end. "
                                                 "This will not crash the program, however, and you can still expect to "
                                                 "see the main output (CSV files) from Genie.", default="n")


# CHECKING FOR CONDA INSTALL
os.system("echo ${HMM_dir}/bitscores.txt > HMMlib.txt")
file = open("HMMlib.txt")
for i in file:
    location = i.rstrip()

os.system("rm HMMlib.txt")
try:
    bitscores = open(location)
    conda = 1
except FileNotFoundError:
    conda = 0

if conda == 0:
    parser.add_argument('-hmm_dir', type=str, help='directory of HMMs', default="NA")

    parser.add_argument('--R', type=str,
                        help="location of R scripts directory (note: this optional argument requires Rscript to be "
                             "installed on your system). The R scripts directory is in the same directory as the "
                             "Genie.py code", default="NA")


args = parser.parse_args()

# CHECKING ARGUMENTS AND PATHS
cwd = os.getcwd()
print("checking arguments")
if conda == 0:
    if args.hmm_dir != "NA":
        print(".")
    else:
        print("You have not provided the location of the HMM library via the -hmm_dir argument. Please do so, and try "
              "again. The HMM library is found within the same directory as the FeGenie executable.")
        print("Exiting")
        raise SystemExit


if args.bin_dir != "NA":
    binDir = args.bin_dir + "/"
    binDirLS = os.listdir(args.bin_dir)
    print(".")
else:
    print("Looks like you did not provide a directory of ORFs in FASTA format")
    print("Exiting")
    raise SystemExit


if args.bin_ext != "NA":
    print(".")
else:
    print('Looks like you did not provide an extension for your genomes/bins or assemblies, so FeGenie does not know'
          ' which files in the provided directory are fasta files that you would like analyzed.')
    print("Exiting")
    raise SystemExit


if args.makeplots == 'y':
    if conda == 0:
        if args.R != "NA":
            print(".")
        else:
            print('Looks like you told GeoGenie to automatically generate R plots. '
                  'However, you have not provided the location of the directory that contains the R scripts '
                  '(as required of you because you did not go through the conda-based installation.')
            print("Exiting")
            raise SystemExit


if conda == 1:
    os.system("echo ${HMM_dir} > HMMlib.txt")
    file = open("HMMlib.txt")
    for i in file:
        HMMdir = (i.rstrip())
        HMMdirLS = os.listdir(HMMdir)
    os.system("rm HMMlib.txt")
else:
    HMMdir = args.hmm_dir
    HMMdirLS = os.listdir(args.hmm_dir)

try:
    os.listdir(args.outdir)
    print("Looks like you already have a directory with the name: " + args.outdir)
    print("To avoid overwriting potentially valuable files, FeGenie will now exit. "
          "Please delete or rename said directory and try running again.")
    print("Exiting")
    raise SystemExit
except FileNotFoundError:
    print(".")
    os.system("mkdir %s" % args.outdir)
    outDirectory = "%s" % args.outdir
    outDirectoryLS = os.listdir("%s" % args.outdir)

numHMMs = 0
for i in HMMdirLS:
    if lastItem(i.split(".")) == "hmm":
        numHMMs += 1

# STARTING MAIN ALGORITHM
bits = open(HMMdir + "/bitscores.txt", "r")
bitDict = defaultdict(lambda: defaultdict(lambda: 'EMPTY'))
print("\nreading in HMM bitscore cut-offs...")
for i in bits:
    ls = i.rstrip().split("\t")
    bitDict[ls[0]]["bit"] = ls[1]
print("...")

count = 0
BinDict = defaultdict(lambda: defaultdict(lambda: 'EMPTY'))
out = open("%s/%s.csv" % (args.outdir, args.out), "w")
out.write("bin" + "," + "gene" + "," + "ORF" + "," + "evalue" + "," + "bitscore" + "," + "bitscore_cutoff" + "," + "seq" + "\n")

outTSV = open("%s/%s.tsv" % (args.outdir, args.out), "w")
outTSV.write("bin" + "\t" + "gene" + "\t" + "ORF" + "\t" + "evalue" + "\t" + "bitscore" + "\t" + "bitscore_cutoff" + "\t" + "seq" + "\n")

for i in binDirLS:
    HMMdict = defaultdict(lambda: defaultdict(lambda: 'EMPTY'))
    if not re.match(r'^\.', i) and i != (args.out + ".csv") and lastItem(i.split(".")) == args.bin_ext:
        print("")
        fastaFile = open(args.bin_dir + "/" + i, "r")
        fastaFile = fasta(fastaFile)
        os.system("mkdir " + args.bin_dir + "/" + i + "-HMM")
        count = 0
        for hmm in HMMdirLS:
            if lastItem(hmm.split(".")) == "hmm":
                count += 1
                perc = (count/numHMMs)*100-1
                sys.stdout.write("analyzing " + i + ": %d%%   \r" % (perc+1))
                sys.stdout.flush()
                if not re.match(r'^\.', hmm):
                    os.system("hmmsearch "
                              "--tblout " + binDir + "/" + i + "-HMM/" + i + "__" + hmm +
                              " -o " + binDir + "/" + i + "-HMM/" + i + "__" + hmm + ".txt " +
                              HMMdir + "/" + hmm + " " +
                              binDir + "/" + i)
                    os.system("rm " + binDir + "/" + i + "-HMM/" + i + "__" + hmm + ".txt")

        HMMresults = os.listdir(args.bin_dir + "/" + i + "-HMM")
        for result in HMMresults:
            if not re.match(r'^\.', result):
                file = open(args.bin_dir + "/" + i + "-HMM/" + result)
                for line in file:
                    line = line.rstrip()
                    if not re.match(r'^#', line):
                        ls = (delim(line))
                        orf = ls[0]
                        gene = ls[2]
                        evalue = ls[4]
                        bitscore = ls[5]
                        hmmFileName = (result.split("__")[1])
                        hmmFileName = hmmFileName.split(".hm")[0]
                        hmmName = deExt(result.split("__")[1])
                        threshold = bitDict[hmmFileName]["bit"]
                        process = bitDict[hmmFileName]["process"]
                        element = bitDict[hmmFileName]["element"]
                        if float(threshold) > 0:
                            if float(bitscore) > float(threshold):
                                if orf not in HMMdict.keys():
                                    HMMdict[orf]["hmm"] = hmmName
                                    HMMdict[orf]["evalue"] = evalue
                                    HMMdict[orf]["bitscore"] = bitscore
                                    HMMdict[orf]["process"] = process
                                    HMMdict[orf]["element"] = element
                                else:
                                    if float(bitscore) > float(HMMdict[orf]["bitscore"]):
                                        HMMdict[orf]["hmm"] = hmmName
                                        HMMdict[orf]["evalue"] = evalue
                                        HMMdict[orf]["bitscore"] = bitscore
                                        HMMdict[orf]["process"] = process
                                        HMMdict[orf]["element"] = element
                                    else:
                                        pass

        for key in HMMdict.keys():
            out.write(i + "," + HMMdict[key]["hmm"] + "," + key + "," + HMMdict[key]["evalue"] + "," + HMMdict[key]["bitscore"] + "," + bitDict[HMMdict[key]["hmm"]]["bit"] + "," + fastaFile[key] + "\n")

            outTSV.write(i + "\t" + HMMdict[key]["hmm"] + "\t" + key + "\t" + HMMdict[key]["evalue"] + "\t" + HMMdict[key][
                "bitscore"] + "\t" + bitDict[HMMdict[key]["hmm"]]["bit"] + "\t" + fastaFile[key] + "\n")

        os.system("rm -r " + args.bin_dir + "/" + i + "-HMM")

out.close()
print("")

#
# summaryDict = defaultdict(lambda: defaultdict(lambda: defaultdict(lambda: 'EMPTY')))
# summary = open("%s/%s.csv" % (args.outdir, args.out), "r")
# for i in summary:
#     ls = i.rstrip().split(",")
#     if ls != "bin":
#         genome = ls[0]
#         gene = ls[1]
#         process = ls[2]
#         substrate = ls[3]
#         orf = ls[4]
#         evalue = ls[5]
#         bitscore = ls[6]
#         sequence = ls[7]
#         bitcut = bitDict[ls[1]]["bit"]
#         summaryDict[genome][orf]["gene"] = gene
#         summaryDict[genome][orf]["process"] = process
#         summaryDict[genome][orf]["substrate"] = substrate
#         summaryDict[genome][orf]["evalue"] = evalue
#         summaryDict[genome][orf]["bitscore"] = bitscore
#         summaryDict[genome][orf]["sequence"] = sequence
#         summaryDict[genome][orf]["bitcut"] = bitcut
#
#
# # ****************************** CLUSTERING OF ORFS BASED ON GENOMIC PROXIMITY *************************************
#
# os.system("rm %s/%s.csv" % (args.outdir, args.out))
# os.system("rm %s/%s-2.csv" % (args.outdir, args.out))
# os.system("mv %s/%s-3.csv %s/%s-summary.csv" % (args.outdir, args.out, args.outdir, args.out))
# # ****************************** CREATING A HEATMAP-COMPATIBLE CSV FILE *************************************
# print("Writing heatmap-compatible CSV")
cats = ["Proteorhodopsin", "Xantharhodopsin", "Heliorhodopsin", "Heliorhodopsin_Pfam", "Bac_rhodopsin",
        "GpcrRhopsn4", "Htr2", "7tm_1", "ASRT.hmm"]

Dict = defaultdict(lambda: defaultdict(list))
final = open("%s/%s.csv" % (args.outdir, args.out), "r")
for i in final:
    ls = (i.rstrip().split(","))
    if ls[0] != "" and ls[1] != "assembly" and ls[1] != "genome":
        if not re.match(r'#', i):
            cell = ls[0]
            gene = ls[1]
            Dict[cell][gene].append(gene)

normDict = defaultdict(lambda: defaultdict(lambda: 'EMPTY'))
for i in os.listdir(args.bin_dir):
    if lastItem(i.split(".")) == args.bin_ext:
        file = open("%s/%s" % (args.bin_dir, i), "r")
        file = fasta(file)
        normDict[i] = len(file.keys())

outHeat = open("%s/%s.heatmap.csv" % (args.outdir, args.out), "w")
outHeat.write("X" + ',')
for i in sorted(Dict.keys()):
    if (not re.match(r'#', i) and i != "bin"):
        outHeat.write(i + ",")
outHeat.write("\n")

for i in cats:
    outHeat.write(i + ",")
    for j in sorted(Dict.keys()):
        if (not re.match(r'#', j) and j != "bin"):
            outHeat.write(str((len(Dict[j][i]) / int(normDict[j])) * float(100)) + ",")
    outHeat.write("\n")

outHeat.close()
print('......')
print(".......")
print("Finished!")

# ******** RUNNING RSCRIPT TO GENERATE PLOTS **************
if args.makeplots == "y":
    if conda == 0:
        Rdir = args.R
    else:
        os.system("echo ${rscripts} > r.txt")
        file = open("r.txt")
        for i in file:
            Rdir = (i.rstrip())
        os.system("rm r.txt")

    os.system("Rscript -e 'install.packages(\"ggplot2\", repos = \"http://cran.us.r-project.org\")\'")
    os.system("Rscript -e 'install.packages(\"reshape\", repos = \"http://cran.us.r-project.org\")\'")
    os.system("Rscript -e 'install.packages(\"reshape2\", repos = \"http://cran.us.r-project.org\")\'")
    os.system("Rscript -e 'install.packages(\"tidyverse\", repos = \"http://cran.us.r-project.org\")\'")
    os.system("Rscript -e 'install.packages(\"argparse\", repos = \"http://cran.us.r-project.org\")\'")
    os.system("Rscript -e 'install.packages(\"ggdendro\", repos = \"http://cran.us.r-project.org\")\'")
    os.system("Rscript -e 'install.packages(\"ggpubr\", repos = \"http://cran.us.r-project.org\")\'")
    os.system("Rscript -e 'install.packages(\"grid\", repos = \"http://cran.us.r-project.org\")\'")

    os.system("Rscript --vanilla %s/DotPlot.R %s/%s.%s.heatmap.csv %s" % (Rdir, args.outdir, args.out, args.element, args.outdir))
    os.system("Rscript --vanilla %s/dendro-heatmap.R %s/%s.%s.heatmap.csv %s" % (Rdir, args.outdir, args.out, args.element, args.outdir))


