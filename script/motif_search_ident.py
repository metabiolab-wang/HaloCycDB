#!/usr/bin/env python
#-*-coding: utf8 -*-

#from Bio.Data import IUPACData
#from Bio import motifs
#from Bio.Seq import Seq
from Bio import SeqIO
import argparse
from optparse import OptionParser
#from symbol import argument
#import sys
import re
import os


def motif_search(genename, L2, L1, seqs, ident, threshold):
    if L1 == 'Halogenase':
        if L2 == 'Vanadium-Dependent Haloperoxidases':
            #### Vanadium dependant halogenase 
            vanadium_motif1 = re.compile(r'[RWC]\w{2,3}K\w{6,7}R', re.U)
            vanadium_motif2 = re.compile(r'P\w{2,4}[SAI]GH')
            mobj1 = re.search(vanadium_motif1, seqs)
            mobj2 = re.search(vanadium_motif2, seqs)
            if mobj1 or mobj2:
                return('TRUE')
            else:
                return('FALSE')

        elif L2 == 'Flavin-Dependent Halogenases':
            #### Flavin dependant halogenase
            flavin_motif1 = re.compile(r'W\w{1,2}[WG]\w{1,2}IP', re.U)
            flavin_motif2 = re.compile(r'Y\w{3}F\w{38}A', re.U)
            mobj1 = re.search(flavin_motif1, seqs)
            mobj2 = re.search(flavin_motif2, seqs)
            if mobj1 or mobj2:
                return('TRUE')
            else:
                return('FALSE')
        
        elif L2 == 'S-Adenosyl-L-Methionine-Dependent Halogenases':
            #### SAM-dependant halogenase
            sam_motif = re.compile(r'YP\w{1,2}TGT', re.U)
            mobj = re.search(sam_motif, seqs)
            if mobj:
                return('TRUE')
            else:
                return('FALSE')
        
        elif L2 == 'Nonheme Iron-Dependent Halogenases':
            #### Nonheme-dependant halogenase
            nonheme_motif1 = re.compile(r'WHQ[AVS]', re.U)
            nonheme_motif2 = re.compile(r'H\wA\w{3,6}H', re.U)
            nonheme_motif3 = re.compile(r'H\wG\w{9,11}H', re.U)
            nonheme_motif4 = re.compile(r'WH[WV]G', re.U)
            mobj1 = re.search(nonheme_motif1, seqs)
            mobj2 = re.search(nonheme_motif2, seqs)
            mobj3 = re.search(nonheme_motif3, seqs)
            mobj4 = re.search(nonheme_motif4, seqs)
            if mobj1 or mobj2 or mobj3 or mobj4:
                return('TRUE')
            else:
                return('FALSE')
    
    elif L1 == 'Dehalogenase':
        if L2 == 'Reductive dehalogenase':
            #### reducative dehalogenase
            #### pcpE cdrA indepent
            #rdh_other = ('CdrA', 'PcpE')
            rdh_other = ('CdrA')
            if genename in rdh_other:
                if float(ident) > threshold:
                    return('TRUE')
            else:
                rdh_c1_motif = re.compile(r'P\w{14}G\w{4}G', re.U)
                rdh_c2_motif = re.compile(r'Y\w{12}F\w{4}GY', re.U)
                rdh_c3_motif = re.compile(r'G\wGE\w{2}R', re.U)
                rdh_c4_motif = re.compile(r'T\w{7}[KR]P', re.U)

                rdh_FeS_motif1 = re.compile(r'C\w{2}C\w{2}C\w{3}CP', re.U)
                rdh_FeS_motif2 = re.compile(r'W\w{2}D\w{2}[KR]C', re.U)
                rdh_FeS_motif3 = re.compile(r'C\w{2}C\w{3}C[SP]', re.U)

                #c1_mobj = re.search(rdh_c1_motif, seqs)
                #c2_mobj = re.search(rdh_c2_motif, seqs)
                #c3_mobj = re.search(rdh_c3_motif, seqs)
                #c4_mobj = re.search(rdh_c4_motif, seqs)
                count=0
                patterns = [rdh_c1_motif, rdh_c2_motif, rdh_c3_motif, rdh_c4_motif]
                
                for patt in patterns:
                    if re.search(patt, seqs):
                        count +=1                

                fes1_mobj = re.search(rdh_FeS_motif1, seqs)
                fes2_mobj = re.search(rdh_FeS_motif2, seqs)
                fes3_mobj = re.search(rdh_FeS_motif3, seqs)

                if count >=2  and (fes1_mobj or fes2_mobj or fes3_mobj):
                    return('TRUE')
                else:
                    return('FALSE')
        
        elif L2 == 'Halohydrin dehalogenases':
            #### Halohydrin dehalogenase
            halohydrin_motif = re.compile(r'S\w{12,13}Y\w{3}R', re.U)
            mobj = re.search(halohydrin_motif, seqs)
            if mobj:
                return('TRUE')
            else:
                return('FALSE')
            
        elif L2 == 'Methyltransferase' or L2 == 'Dehydrochlorinase':
            #### Methyltransferase #MecE MecF CumA
            #### dehydrochlorinase #LinA1, LinA2, LinA1
            if float(ident) > threshold:
                return('TRUE')
        
        elif L2 == 'Glutathione S-transferase':
            #### Glutathione S-transferase
            glutathione_motif = re.compile(r'S\w{1,3}C', re.U)
            glutathione_zeta_exclude_motif = re.compile(r'I\w{3}F\w{2}LE\w{12}GD', re.U)

            mobj1 = re.search(glutathione_motif, seqs)
            mobj2 = re.search(glutathione_zeta_exclude_motif, seqs)

            if mobj1 and not mobj2:
                return('TRUE')
            else:
                return('FALSE')
        
        elif L2 == 'Hydrolytic dehalogenase' or L2 == 'Oxidative dehalogenase':
            #### hydrolytic halogenase and oxidative dehalogenase
            HAD_haloacid_L = ('haloalkanoic', 'L-DEX', 'DehH2', 'L-2-haloacid-II','L-2-haloacid-I', 'Bpro-0530',
                  'dehH109', 'HadL', 'Dhl-VII', 'DhlB', 'RPA4199', 'Bpro-4516','RSc1362', 'Adeh3811',
                 'Dhe', 'HdI-IVa', 'ST2570', 'Chd1', 'zgHAD', 'PH0459', 'RHA1-ro00230', 'RPA2507', 'PA0810',
                 'DehL')
            
            HAD_haloacid_D = ('DehI', 'DehD', 'Dehc2', 'DL-DEX', 'DehE', 'DhlC', 'DehHX', 'HadD')

            HAD_haloacid_chd = ('Chd')

            ABH_haloalkane_I_II = ('DhlA','DhmA', 'DhmA1', 'DhmA2', 'DrbA', 'Xf1965', 'CurN' 'JANN2620', 'DbeA', 'DadB', 'HanR', 'DhaA', 'LinB')

            ABH_fluoroacetate = ('BbdC', 'RPA1163', 'Alr0039', 'DehH1', 'FAc-DEX')

            oxidative_monooxygenases_I = ('TftD', 'HadA', 'TcpA', 'HcbB3', 'CnpA')

            oxidative_monooxygenases_II = ('DHPA', 'DHPB', 'PcpB', 'HcbA1', 'PcpA')

            hydrolytic_other1 = ('AtzA', 'TrzN')

            hydrolytic_other2 = ('Cis-caaD')

            hydrolytic_other3 = ('4-CBA-CoA-dehalogenase', 'FcbB', 'FcbB2', 'FcbB1')

            if genename in HAD_haloacid_L:
                # L sytle haloacid include "DehL"
                HAD_haloacid_L_motif = re.compile(r'FD.+(S\wN\w{2}|A\wH\w{2})D',re.U)
                mobj = re.search(HAD_haloacid_L_motif, seqs)
                if mobj:
                    return('TRUE')
                else:
                    return('FALSE')
                
            elif genename in HAD_haloacid_D:
                # D style haloacide
                #HAD_haloacid_D_motif1 = re.compile(r'FI\w{2}I', re.U)
                #HAD_haloacid_D_motif2 = re.compile(r'FL\w{2}[LS]', re.U)
                #HAD_haloacid_D_motif3 = re.compile(r'ML\w{2}L', re.U)
                #mobj1 = re.search(HAD_haloacid_D_motif1, line)
                #mobj2 = re.search(HAD_haloacid_D_motif2, line)
                #mobj3 = re.search(HAD_haloacid_D_motif3, line)
                #if mobj1 or mobj2 or mobj3:
                #    print("match2")
                #else:
                #    print("no match2")

                HAD_haloacid_D_motif = re.compile(r'VPW(V\wF|M\wV)', re.U)
                mobj = re.search(HAD_haloacid_D_motif, seqs)
                if mobj:
                    return('TRUE')
                else:
                    return('FALSE')

            elif genename in HAD_haloacid_chd:
                chd_motif =  re.compile(r'SYHGDH', re.U)
                mobj = re.search(chd_motif, seqs)
                if mobj:
                    return('TRUE')
                else:
                    return('FALSE')

            elif genename in ABH_haloalkane_I_II:
                ABH_haloalkane_I_II_motif = re.compile(r'DWG', re.U)
                mobj = re.search(ABH_haloalkane_I_II_motif, seqs)
                if mobj:
                    return('TRUE')
                else:
                    return('FALSE')
            
            elif genename in ABH_fluoroacetate:
                ABH_fluoroacetate_motif = re.compile(r'DRG', re.U)
                mobj = re.search(ABH_fluoroacetate_motif, seqs)
                if mobj:
                    return('TRUE')
                else:
                    return('FALSE')
                
            elif genename in oxidative_monooxygenases_I:
                oxidative_monooxygenases_I_motif = re.compile(r'DPQ\wDR.*HY[HE]', re.U)
                mobj = re.search(oxidative_monooxygenases_I_motif, seqs)
                if mobj:
                    return('TRUE')
                else:
                    return('FALSE')

            elif genename in oxidative_monooxygenases_II:
                if float(ident) > threshold:
                    return('TRUE')
            
            elif genename in hydrolytic_other1:
                hydrolytic_other1 = ('AtzA', 'TrzN')
                hydrolytic_other1_motif = re.compile(r'[NCD]\wH[TNQAVH]H', re.U)
                mobj = re.search(hydrolytic_other1_motif, seqs)
                if mobj:
                    return('TRUE')
                else:
                    return('FALSE')

            elif genename in hydrolytic_other2:
                hydrolytic_other2_motif = re.compile(r'R\w{2}R', re.U)
                mobj = re.search(hydrolytic_other2_motif, seqs)
                if mobj:
                    return('TRUE')
                else:
                    return('FALSE')
                
            elif genename in hydrolytic_other3:
                hydrolytic_other3_motif1 = re.compile(r'A\w{4}H', re.U)
                hydrolytic_other3_motif2 = re.compile(r'W\w{7}D', re.U)
                mobj1 = re.search(hydrolytic_other3_motif1, seqs)
                mobj2 = re.search(hydrolytic_other3_motif2, seqs)
                if mobj1 and mobj2:
                    return('TRUE')
                else:
                    return('FALSE')




def filter_with_motif(fasta_in, annotate_in, ident, threshold, output):

    fa = {}
    for seqs in SeqIO.parse(fasta_in, "fasta"):
        fa[seqs.id] = seqs.seq

    idt = {}
    with open(ident) as ident_f:
        for line in ident_f:
            line = line.strip()
            idt[re.split('\\s+', line.split('\t')[0])[0]] = line.split('\t')[2]

    fname= os.path.basename(annotate_in).split(".")[0]

    seqs_output = open(output + 'halocyc_filter.fa', 'w')
    parse_output = open(output + fname + '_filter.txt', 'w')
    with open(annotate_in) as annotate_f:
        for line in annotate_f:
            line = line.strip()
            lines = line.split('\t')
            if motif_search(lines[2], lines[4], lines[5], str(fa[lines[0]]), idt[lines[0]], (threshold * 100)) == 'TRUE':
                seqs_output.writelines(">" + lines[0] + "\n")
                seqs_output.writelines(fa[lines[0]] + '\n')
                parse_output.writelines(line + '\n')
    
    seqs_output.close()
    parse_output.close()
    annotate_f.close()
                

if __name__ == "__main__":
    parser = OptionParser()
    parser.add_option("-i", "--input",type="string",dest="input",action="store",help="input fasta file")
    parser.add_option("-a", "--annotate",type="string",dest="annotate",action="store",help="halocyc parse file")
    parser.add_option('-o', "--output",type="string",dest="output",action="store",help="output directory, the output file 'halocyc_genename_parse_filter.txt' and 'halocyc_filter.fa' will be obtained in the output directory")
    parser.add_option('-b', "--bresult",type="string",dest="bresult",action="store",help="annotation result file")
    parser.add_option('-t', "--threshold",type="float",dest="threshold",action="store",help="threshold of sequence similarity, range (0~1), default=0.4", default=0.4)

    (argvs,args)= parser.parse_args()

    if not argvs.input and argvs.annotate and argvs.output and argvs.bresult:
        parser.print_help()
    else:
        filter_with_motif(argvs.input, argvs.annotate, argvs.bresult, argvs.threshold, argvs.output)
