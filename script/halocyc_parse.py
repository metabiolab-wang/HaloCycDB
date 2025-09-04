#!/usr/bin/env python
# -*- coding: utf-8 -*-

#import argparse
from optparse import OptionParser
#from symbol import argument
#from Bio import SeqIO
import sys
#import re
import os

args = sys.argv

def halocyc_gen(halocycgene):
    halocyc_dic = {}
    with open(halocycgene) as f:
        for line in f:
            line = line.strip()
            line_split = line.split("\t")
            new_line = line_split[6] + "\t" +line_split[7] + "\t" + line_split[4] + "\t" + line_split[5]
            halocyc_dic[line_split[0]] = new_line
    f.close()
    return halocyc_dic

def abund(abundance):
    abund_dic = {}
    with open(abundance) as bin_f:
        for line in bin_f:
            line = line.strip()
            line_split = line.split("\t")
            abund_dic[line_split[0]] = line_split[1]
    bin_f.close()
    return  abund_dic

def tax(taxon, style):
    tax_dic = {}
    if str(style) == 'bin':
        bac_taxon = os.path.join(taxon, "gtdbtk.bac120.summary.tsv")
        if os.path.exists(bac_taxon):
            with open(bac_taxon) as bac_f:
                for line in bac_f:
                    line = line.strip()
                    line_split = line.split("\t",2)
                    tax_dic[line_split[0]] = line_split[1]
            bac_f.close()

        arc_taxon = os.path.join(taxon, "gtdbtk.ar53.summary.tsv")
        if os.path.exists(arc_taxon):
            with open(arc_taxon) as arc_f:
                for line in arc_f:
                    line = line.strip()
                    line_split = line.split("\t",2)
                    tax_dic[line_split[0]] = line_split[1]
            arc_f.close()

    elif str(style) == 'contig':
        with open(taxon) as f:
            for line in f:
                line = line.strip()
                line_split = line.split("\t")
                if len(line_split) > 1:
                    tax_dic[line_split[0]] = line_split[1]
                else:
                    tax_dic[line_split[0]] = 'Unclassfied'

        f.close()

    return tax_dic



def halocyc_parse(input,halocycgene,abundance,taxon,style,output):
    halocyc_dic = halocyc_gen(halocycgene)
    bins_abund_dic = abund(abundance)
    bins_tax = tax(taxon,style)
    if str(style) == 'bin':
        current_path = os.path.abspath(input)
        current_path_split = current_path.rsplit("/",2)[1]
        gene_output = open(output, 'w')
        with open(input) as input_f:
            for line in input_f:
                line = line.strip()
                line_split = line.split("\t")
                gene_output.writelines(line + "\t" + bins_abund_dic[current_path_split] + "\t" + halocyc_dic[line_split[2]] +"\t" + bins_tax[current_path_split] + "\n")
        input_f.close()

    if str(style) == 'contig':
        gene_output = open(output, 'w')
        with open(input) as input_f:
            for line in input_f:
                line = line.strip()
                line_split = line.split("\t")
                line_name = line_split[0].rsplit("_",1)[0]
                gene_output.writelines(line_split[0] + "\t" + line_split[1]+ "\t" + line_split[2] +"\t" + bins_abund_dic[line_name] + "\t" + halocyc_dic[line_split[2]] +"\t" + bins_tax[line_name] + "\n")
        input_f.close()



def main(args):
    parser = OptionParser()
    parser.add_option("-d", "--database",type="string",dest="database",action="store",help="database file halocyc.gene")
    parser.add_option("-i", "--input",type="string",dest="input",action="store",help="input file (halocyc_genename.txt)")
    parser.add_option("-o", "--output",type="string",dest="output",action="store",help="output file")
    parser.add_option("-t", "--taxon",type="string",dest="taxon",action="store",help="taxon annotation file (style for contig) or taxon output directory of gtdbtk (style for bin)")
    parser.add_option("-s", "--style",type="string",dest="style",action="store",help="style, 'bin' or 'contig',consistent with the type of protein sequences file used for annotation")
    parser.add_option("-a", "--abund",type="string",dest="abund",action="store",help="abundance file of bin or contig")
    (argvs, args)= parser.parse_args()

#    if not argvs.d and argvs.i and argvs.o and argvs.s and argvs.t and argvs.a:
#        parser.print_help()
#    else:
    halocyc_parse(argvs.input, argvs.database, argvs.abund, argvs.taxon, argvs.style, argvs.output)

if __name__ == "__main__":
    main(args)
