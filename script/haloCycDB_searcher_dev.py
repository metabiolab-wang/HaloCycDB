#!/usr/bin/env python
#-*-coding: utf8 -*-

from optparse import OptionParser
#from symbol import argument
from Bio import SeqIO
from multiprocessing import Pool, Manager
import sys
import os
import shutil
import re
import time


args = sys.argv


def batch_iterator(iterator, batch_size):
    """Returns lists of length batch_size.

    This can be used on any iterator, for example to batch up
    SeqRecord objects from Bio.SeqIO.parse(...), or to batch
    Alignment objects from Bio.Align.parse(...), or simply
    lines from a file handle.

    This is a generator function, and it returns lists of the
    entries from the supplied iterator.  Each list will have
    batch_size entries, although the final list may be shorter.
    """
    batch = []
    for entry in iterator:
        batch.append(entry)
        if len(batch) == batch_size:
            yield batch
            batch = []


def chunk_fasta(input,size):
    record_iter = SeqIO.parse(open(input), "fasta")
    current_path = os.getcwd()
    os.mkdir(os.path.join(current_path, "chunk_files"))
    dir_out = os.path.join(current_path, "chunk_files")
    base = os.path.basename(input)
    fname = os.path.splitext(base)[0]
    for i, batch in enumerate(batch_iterator(record_iter,size)):
        out_chunk_name = os.path.join(dir_out, "{0}_chunk{1}.fasta".format(fname, i+1))
        with open(out_chunk_name, "w") as handle:
            SeqIO.write(batch, handle, "fasta")


def pickseqs(input,output,seqs_name, lock):
    current_path = os.getcwd()
    dir_out = os.path.join(current_path, "chunk_files")
    seqs_input = os.path.join(dir_out, input)
    seqs_output = open(output, 'a')
    for seqs in SeqIO.parse(seqs_input, "fasta"):
        if str('>' + seqs.description) in seqs_name.keys():
            lock.acquire()
            seqs_output = open(output, 'a')
            seqs_output.writelines(">" + seqs.description + "\n")
            seqs_output.writelines(seqs.seq + "\n")
            seqs_output.close()
            lock.release()

def seq_output(input, filter_output, bresult_dict, method):
    seqs_output = open(filter_output, 'a')
    for seqs in SeqIO.parse(input, "fasta"):
        if str(seqs.id) in bresult_dict.keys():
            seqs_output.writelines(">" + seqs.id + "\n")
            seqs_output.writelines(seqs.seq + "\n")
    seqs_output.close()




def exist_tool(name):
    """check whether 'name' exist on PATH."""
    from shutil import which
    return which(name) is not None


def time_use(stime, etime):
    usetime = etime - stime
    Time_h = usetime // 3600
    Time_m = usetime // 60
    Time_s = usetime % 60
    print('......................................')
    print("Time used: " + str(int(Time_h)) + "h " + str(int(Time_m)) + "m " + str(round(Time_s,2)) + "s")


def blastseq(method, input, blast_output, database, evalue, similarity, process):
    database_dehalo = os.path.join(database, "dehalogenase/dehalogenase.fa")
    database_halo   = os.path.join(database, "halogenase/halogenase.fa")
    blast_output_dehalo = os.path.join(blast_output, 'dehalo.search')
    blast_output_halo   = os.path.join(blast_output, 'halo.search')
    os.system('mkdir ' + str(blast_output))

    if str(method) == 'usearch':
        if not exist_tool("usearch"):
            print("usearch can not find in system path, please installed before execute this script")
        else:
            print('Starting blast using usearch software')
            print('usearch -usearch_global ' + str(input) + ' -db ' + str(database_halo) + ' -id ' + str(similarity) + ' -blast6out ' + str(blast_output_halo) + ' -threads ' + str(process))
            os.system('usearch -usearch_global ' + str(input) + ' -db ' + str(database_halo) + ' -id ' + str(similarity) + ' -blast6out ' + str(blast_output_halo) + ' -threads ' + str(process))
            print('usearch -usearch_global ' + str(input) + ' -db ' + str(database_dehalo) + ' -id ' + str(similarity) + ' -blast6out ' + str(blast_output_dehalo) + ' -threads ' + str(process))
            os.system('usearch -usearch_global ' + str(input) + ' -db ' + str(database_dehalo) + ' -id ' + str(similarity) + ' -blast6out ' + str(blast_output_dehalo) + ' -threads ' + str(process))
            print("usearch for blast is finished!")

    elif str(method) == 'diamond':
        if not exist_tool("diamond"):
            print("diamond can not find in system path, please installed before execute this script")
        
        else:
            diamond_halo_db =  os.path.join(database, "halogenase/halogenase.dmnd")
            diamond_dehalo_db = os.path.join(database, "dehalogenase/dehalogenase.dmnd")
            
            if not os.path.exists(diamond_halo_db):
                print("Built the diamond format database for HaloCyc DB")
                print('diamond makedb --in ' + str(database_halo) + ' --db ' + str(diamond_halo_db))
                os.system('diamond makedb --in ' + str(database_halo) + ' --db ' + str(diamond_halo_db))
            
            print('Starting blast using diamond software')
            print('diamond blastp -k 1 -e ' + str(evalue) + ' -p ' + str(process) + ' -d ' + str(diamond_halo_db) + ' -q ' + str(input) + ' -a halogenase')
            os.system('diamond blastp -k 1 -e ' + str(evalue) + ' -p ' + str(process) + ' -d ' + str(diamond_halo_db) + ' -q ' + str(input) + ' -a halogenase')
            os.system('diamond view -a halogenase.daa -o tmp.halogenase.diamond.out')
            os.system("cat tmp.halogenase.diamond.out |awk -F '\\t' '{if ($3 > " + str(similarity * 100) + "){print $0}}' > " + str(blast_output_halo))
            os.system('rm tmp.halogenase.diamond.out halogenase.daa')
            print("diamond for blast is finished!")

            if not os.path.exists(diamond_dehalo_db):
                print("Built the diamond format database for HaloCyc DB")
                print('diamond makedb --in ' + str(database_dehalo) + ' --db ' + str(diamond_dehalo_db))
                os.system('diamond makedb --in ' + str(database_dehalo) + ' --db ' + str(diamond_dehalo_db))
            
            print('Starting blast using diamond software')
            print('diamond blastp -k 1 -e ' + str(evalue) + ' -p ' + str(process) + ' -d ' + str(diamond_dehalo_db) + ' -q ' + str(input) + ' -a dehalogenase')
            os.system('diamond blastp -k 1 -e ' + str(evalue) + ' -p ' + str(process) + ' -d ' + str(diamond_dehalo_db) + ' -q ' + str(input) + ' -a dehalogenase')
            os.system('diamond view -a dehalogenase.daa -o tmp.dehalogenase.diamond.out')
            os.system("cat tmp.dehalogenase.diamond.out |awk -F '\\t' '{if ($3 > " + str(similarity * 100) + "){print $0}}' > " + str(blast_output_dehalo))
            os.system('rm tmp.dehalogenase.diamond.out dehalogenase.daa')
            print("diamond for blast is finished!")
    
    elif str(method) == 'blastp':
        if not exist_tool("blastp"):
            print("blastp can not find in system path, please installed before execute this script")

        else:
            blastp_halo_db  = os.path.join(database, "halogenase/halogenase")
            blastp_dehalo_db = os.path.join(database, "dehalogenase/dehalogenase")

            if not os.path.exists((blastp_halo_db + '.psq')):
                print("Built the blast format database for halogenase DB")
                print('makeblastdb -in ' + str(database_halo) + ' -dbtype prot -parse_seqids -out ' + str(blastp_halo_db))
                os.system('makeblastdb -in ' + str(database_halo) + ' -dbtype prot -parse_seqids -out ' + str(blastp_halo_db))
            
            print('Starting blast using blastp software')
            print('blastp -query ' + str(input) + ' -out ' + str(blast_output_halo) + ' -db ' + str(blastp_halo_db) + ' -outfmt 6 -max_target_seqs 1 -evalue ' + str(evalue) + ' -num_threads ' + str(process))
            os.system('blastp -query ' + str(input) + ' -out tmp.blastp.out' + ' -db ' + str(blastp_halo_db) + ' -outfmt 6 -max_target_seqs 1 -evalue ' + str(evalue) + ' -num_threads ' + str(process))
            os.system("cat tmp.blastp.out |awk -F '\\t' '{if ($3 > " + str(similarity * 100) + "){print $0}}' > " + str(blast_output_halo))
            os.system('rm tmp.blastp.out')

            if not os.path.exists((blastp_dehalo_db + '.psq')):
                print("Built the blast format database for dehalogenase DB")
                print('makeblastdb -in ' + str(database_dehalo) + ' -dbtype prot -parse_seqids -out ' + str(blastp_dehalo_db))
                os.system('makeblastdb -in ' + str(database_dehalo) + ' -dbtype prot -parse_seqids -out ' + str(blastp_dehalo_db))
            
            print('Starting blast using blastp software')
            print('blastp -query ' + str(input) + ' -out ' + str(blast_output_dehalo) + ' -db ' + str(blastp_dehalo_db) + ' -outfmt 6 -max_target_seqs 1 -evalue ' + str(evalue) + ' -num_threads ' + str(process))
            os.system('blastp -query ' + str(input) + ' -out tmp.blastp.out' + ' -db ' + str(blastp_dehalo_db) + ' -outfmt 6 -max_target_seqs 1 -evalue ' + str(evalue) + ' -num_threads ' + str(process))
            os.system("cat tmp.blastp.out |awk -F '\\t' '{if ($3 > " + str(similarity * 100) + "){print $0}}' > " + str(blast_output_dehalo))
            os.system('rm tmp.blastp.out')
            print("blastp for blast is finished!")



def filterseq(input, filter_output, db_path, bresult, ratio, method):
    length = os.path.join(db_path, "halocycdb.len")

    length_dict={}
    with open(length) as lenf:
        for line in lenf:
            line = line.rstrip()
            seqname,gene_length = line.split('\t')
            length_dict[seqname]=gene_length
    lenf.close()
    
    bresult_dict={}
    with open(bresult) as f:
        for line in f:
            line = line.rstrip()
            blist = line.split('\t')
            if method == 'usearch':
                qname = re.split('\\s+', blist[0])[0]
                tname, id, blength = blist[1], blist[2], blist[3]
            else:
                qname, tname, id, blength = blist[0], blist[1], blist[2], blist[3]
            if float(id) > float(ratio * 100) and float(blength) > float(float(length_dict[tname]) * ratio):
                bresult_dict[qname] = line
    f.close()

    seq_output(input, filter_output, bresult_dict, method)

    filter_bresult = open('filter_bresult.txt', 'w')
    for key in bresult_dict.keys():
        filter_bresult.writelines(bresult_dict.get(key) + '\n')
    filter_bresult.close()
    os.system('mv filter_bresult.txt ' + bresult)


def hmmsearchseq(target_dir, database, evalue, domE, process):
    database_dehalo = os.path.join(database, 'dehalogenase/dehalogenase.hmm')
    database_halo   = os.path.join(database, 'halogenase/halogenase.hmm')
    hmm_output_dehalo = os.path.join(target_dir, 'dehalo.hresult')
    hmm_output_halo   = os.path.join(target_dir, 'halo.hresult')
    hmm_input_dehalo = os.path.join(target_dir, 'dehalogenase.fa')
    hmm_input_halo = os.path.join(target_dir, 'halogenase.fa')
    
    if not exist_tool("hmmsearch"):
            print("hmmserch can not find in system path, please installed before execute this script")
    else:
        try:
            if os.path.getsize(hmm_input_halo) > 0:
                print('Starting hmmsearch software')
                print('hmmsearch -E ' + str(evalue) + ' --domE ' + str(domE) + '  --domtblout ' + str(hmm_output_halo) + ' -o /dev/null ' + ' --cpu ' + str(process) + ' ' + str(database_halo) + ' ' + str(hmm_input_halo))
                os.system('hmmsearch -E ' + str(evalue) + ' --domE ' + str(domE) + '  --domtblout ' + str(hmm_output_halo) + ' -o /dev/null ' + ' --cpu ' + str(process) + ' ' + str(database_halo) + ' ' + str(hmm_input_halo))
            else:
                print("Can not find halogenase gene by hmmsearch!")
        except OSError as e:
            pass

        try:
            if os.path.getsize(hmm_input_dehalo) > 0:
                print('hmmsearch -E ' + str(evalue) + ' --domE ' + str(domE) + '  --domtblout ' + str(hmm_output_dehalo) + ' -o /dev/null ' + ' --cpu ' + str(process) + ' ' + str(database_dehalo) + ' ' + str(hmm_input_dehalo))
                os.system('hmmsearch -E ' + str(evalue) + ' --domE ' + str(domE) + '  --domtblout ' + str(hmm_output_dehalo) + ' -o /dev/null ' + ' --cpu ' + str(process) + ' ' + str(database_dehalo) + ' ' + str(hmm_input_dehalo))
            else:
                print("Can not find dehalogenase gene by hmmsearch!")
        except OSError as e:
            pass

        print("hmmsearch for filter is finished!")


def validseq(input, output_genename, output_fa, hresult, bresult, abund, method):
    hmmsearch_result={}
    with open(hresult) as hf:
        for line in hf:
            line = line.rstrip()
            if not line.startswith('#'):
                hqname = re.split('\\s+', line)[0]
                hmmsearch_result[hqname] = 1
    hf.close()

    blast_filter_result={}
    gene_name = {}
    with open(bresult) as bf:
        if method == 'usearch':
            for line in bf:
                line = line.rstrip()
                bqname = re.split('\\s', line)[0]
                if bqname in hmmsearch_result.keys():
                    blast_filter_result[bqname] = line
                    gene_name[bqname] = line.split('\t')[1]
        else:
            for line in bf:
                line = line.rstrip()
                bqname = line.split('\t')[0]
                if bqname in hmmsearch_result.keys():
                    blast_filter_result[bqname] = line
                    gene_name[bqname] = line.split('\t')[1]
    bf.close

#    if not os.path.exists(os.path.join(os.getcwd(), "chunk_files")):
#        print("Chunk the fasta file")
#        chunk_fasta(input, 10000)
#    files = [f for f in os.listdir(os.path.join(os.getcwd(), "chunk_files"))]
#    pool = Pool(process)
#    lock=Manager().Lock()
#    [pool.apply_async(pickseqs, args = (name, valid_output, hmmsearch_result, lock)) for name in files]
#    pool.close()
#    pool.join()
#    shutil.rmtree(os.path.join(os.getcwd(), "chunk_files"))

    seq_output(input, output_fa, blast_filter_result, method)
    
    if not abund:
        valid_ow = open(output_genename, 'w')
        for gname in gene_name.keys():
            valid_ow.writelines(gname + "\t" + gene_name.get(gname) + "\t" + gene_name.get(gname).split('_')[-1]+"\n")
        valid_ow.close()

    
    else:
        abund_dict={}
        with open(abund, 'r') as gf:
            for line in gf:
                if not line.startswith('Contig'):
                    line = line.rstrip()
                    seqname, abundt = line.split('\t')
                    abund_dict[seqname] = abundt
        gf.close()

        valid_ow = open(output_genename, 'w') 
        for gname in gene_name.keys():
            genename_new = gname.rsplit("_",1)[0]
            valid_ow.writelines(gname + "\t" + gene_name.get(gname) + "\t" + gene_name.get(gname).split('_')[-1]+ "\t" + abund_dict.get(genename_new) + "\n")
        valid_ow.close()
    
    bresult_new = open(bresult,'w')
    for line in blast_filter_result.keys():
        bresult_new.writelines(blast_filter_result.get(line) + "\n")
    bresult_new.close()


def parseseq(gene_table, info, parse_output):
    hresult_dict={}
    with open(gene_table) as gt:
        for line in gt:
            line = line.rstrip()
            hqname, abundt = line.split(' ')[0,1]
            hresult_dict[hqname] = abundt
    gt.close()

    info_dict={}
    with open(info) as inf:
        for line in inf:
            line = line.rstrip()
            gname = line.split(' ')[0]
            info_dict[gname] = line
    inf.close()
    
    poutput = open(parse_output, 'w')
    for line in hresult_dict.keys():
        poutput.write(info_dict[line] + "\t" + hresult_dict[line] + "\n")
    poutput.close()


def main(args):
    parser = OptionParser(prog='python HaloCycDB_searcher.py',
                            description="HaloCyc database searcher(v.1.0)",
                            usage="%prog [option] args...",
                            version="1.0")
    parser.add_option("-f", "--function",type="string",dest="function",action="store",
                        help="function for select to run: search|parse. search mode: to search and identify the halogenase/dehalogenase genes from protein sequences. parse mode: To maps the sequence IDs to gene information")
    parser.add_option("-i", "--input",type="string",dest="input",action="store",
                        help="protein sequences file derived from bins or contigs after open reading frame prediction")
    parser.add_option("-m", "--method",type="string",dest="method",action="store",
                        help="method selected for annotation: blastp|diamond|usearch; usearch used as default", default="usearch")
    parser.add_option("-o", "--output",type="string",dest="output",action="store",
                        help="output directory")
    parser.add_option("-d", "--database",type="string",dest="database",action="store",
                        help="path to HaloCycDB")
    parser.add_option("-e", "--evalue",type="int",dest="evalue",action="store",
                        help="evalue set for annotation, default=1e-5", default=1e-5)
    parser.add_option("-E", "--domE",type="int",dest="domE",action="store",
                        help="domain evalue set for hmmsearch, default=1e-5", default=1e-5)
    parser.add_option("-s", "--similarity",type="float",dest="similarity",action="store",
                        help="similarity(identity) between query and object sequence, range (0~1), default=0.4", default=0.4)
    parser.add_option("-r", "--ratio",type="float",dest="ratio",action="store",
                        help="The alignment length is at least a certain proportion of the subject gene length, default=0.4", default=0.4)
    parser.add_option("-a", "--abund",type="string",dest="abund",action="store",
                        help="sequence abundance apply for analysis")
    parser.add_option("-n", "--info",type="string",dest="info",action="store",
                        help="HaloCycDB information file (HaloCycDB.info)")
    parser.add_option("-p", "--process",type="int",dest="process",action="store",
                        help="number of processes for annotation, default=4", default=4)
    
    (argvs, args)= parser.parse_args()

    if not argvs.function:
        parser.print_help()
    
    elif str(argvs.function) == 'search':
        if not argvs.input:
            print("Please input query fasta file")
        elif not argvs.output:
            print("Please assign output file name")
        elif not argvs.database:
            print("Please assign the path of halocyc database")
        else:
            try:
                os.remove(argvs.output)
            except OSError:
                pass
            stime=time.time()
            blastseq(argvs.method, argvs.input, argvs.output, argvs.database, argvs.evalue, argvs.similarity, argvs.process)

            halo_fa = os.path.join(argvs.output,"halogenase.fa")
            halo_bresult = os.path.join(argvs.output, "halo.search")
            filterseq(argvs.input, halo_fa, argvs.database, halo_bresult, argvs.ratio, argvs.method)


            dehalo_fa = os.path.join(argvs.output,"dehalogenase.fa")
            dehalo_bresult = os.path.join(argvs.output, "dehalo.search")
            filterseq(argvs.input, dehalo_fa, argvs.database, dehalo_bresult, argvs.ratio, argvs.method)
            
            hmmsearchseq(argvs.output, argvs.database, argvs.evalue, argvs.domE, argvs.process)
            
            halo_genename = os.path.join(argvs.output, "halo_genename.txt")
            dehalo_genename = os.path.join(argvs.output, "dehalo_genename.txt")
            halo_filter_fa = os.path.join(argvs.output, "halo_filter.fa")
            dehalo_filter_fa = os.path.join(argvs.output, "dehalo_filter.fa")
            halo_hresult = os.path.join(argvs.output, "halo.hresult")
            dehalo_hresult = os.path.join(argvs.output, "dehalo.hresult")

            try:
                if os.path.getsize(halo_hresult) > 0 and os.path.getsize(halo_fa) > 0:         
                    validseq(halo_fa, halo_genename, halo_filter_fa, halo_hresult, halo_bresult, argvs.abund, argvs.method)
                else:
                    print("Genes of halogenase are not exist!")
            except OSError:
                pass
            
            try:
                if os.path.getsize(dehalo_hresult) > 0 and os.path.getsize(dehalo_fa) > 0:
                    validseq(dehalo_fa, dehalo_genename, dehalo_filter_fa, dehalo_hresult, dehalo_bresult, argvs.abund, argvs.method)
                else:
                    print("Genes of dehalogenase are not exist!")
            except OSError:
                pass
            
            try:
                if os.path.exists(halo_genename) and os.path.exists(dehalo_genename):
                    os.system('cat ' + halo_genename + ' '+ dehalo_genename + ' > ' + os.path.join(argvs.output, 'halocyc_genename.txt'))
                
                elif os.path.exists(halo_genename):
                    os.system('cat ' + halo_genename + ' > ' + os.path.join(argvs.output, 'halocyc_genename.txt'))
                
                elif os.path.exists(dehalo_genename):
                    os.system('cat ' + dehalo_genename + ' > ' + os.path.join(argvs.output, 'halocyc_genename.txt'))

                else:
                    print("Can not find genes of halocyc in " + argvs.input)
            except OSError:
                pass
            
            etime=time.time()
            time_use(stime, etime)

    elif str(argvs.function) == 'parse':
        if not argvs.gene_table:
            print("Please input gene table")
        elif not argvs.info:
            print("Please input HaloCyc database infomation (HaloCycDB.info)")
        elif not argvs.output:
            print("Please assign output file name")
        else:
            try:
                parse_output = argvs.output + '.txt'
                os.remove(parse_output)
            except OSError:
                pass
            stime=time.time()
            parseseq(argvs.gene_table, argvs.info, parse_output)
            etime=time.time()
            time_use(stime, etime)


if __name__ == "__main__":
    main(args)
