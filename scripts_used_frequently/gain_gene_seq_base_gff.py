# #!/usr/bin/env python3
# -*- encoding:utf-8 -*-


def Get_antisense_strand(line):
    ''' convert the nuclear strand to antisense strand '''
    line = line.upper()
    anti_list = []
    for i in line :
        if i == 'A':
            anti_list.append('T')
        elif i == 'T':
            anti_list.append('A')
        elif i == 'C':
            anti_list.append('G')
        elif i == 'G':
            anti_list.append('C')
        elif i == 'N':
            anti_list.append('N')
    anti_line = ''.join(anti_list)
    return anti_line


def fasta_to_dict(genome_fasta):
    ''' Save fasta file to dictionary '''
    fasta_dict = {}
    with open(genome_fasta) as file:
        for line in file:
            line = line.rstrip()
            if line.startswith('>'):
                dict_key = line
                fasta_dict[dict_key] = ''
            else:
                fasta_dict[dict_key] += line
    return fasta_dict


def Get_seq_at_specific_position(seq,start,end):
    ''' Extract a sequence based on its start and stop positions '''
    seq_fa=seq[start-1:end]
    return seq_fa


def Get_infomation_from_gff(gff3_file,genome_fasta,item,output_file):
    ''' Extract sequences from the genome based on gff files '''
    gff = open(gff3_file,'r')
    result = open(output_file,'w')
    items = item + '	'                         #TODO
    fa_dict = fasta_to_dict(genome_fasta)
    for gff_line in gff:
        gff_line = gff_line.rstrip()
        if items in gff_line:
            gff_line_list = gff_line.split('\t')
            seq_id = '>' + gff_line_list[0]
            source = gff_line_list[1]        #Optionally
            start = int(gff_line_list[3])
            end = int(gff_line_list[4])
            if gff_line_list[5] == '.':      #Optionally
                score = '.'
            if gff_line_list[5] is float:
                score = float(gff_line_list[5])
            if gff_line_list[5] is int:
                score = int(gff_line_list[5])
            strand = gff_line_list[6]
            phase = gff_line_list[7]         #Optionally
            attributes_list = gff_line_list[8].split(';')
            gene_id = '>' + attributes_list[0].rstrip()[3:]
            seq_fa = fa_dict[seq_id]
            gene_seq = Get_seq_at_specific_position(seq_fa,start,end)
            if strand == '+':
                result.write(gene_id+'\n')
                result.write(gene_seq+'\n')
            elif strand == '-':
                coding_strand_seq = Get_antisense_strand(gene_seq)
                result.write(gene_id + '\n')
                result.write(coding_strand_seq[::-1] + '\n')
    result.close()

def main():
    ''' Setting command line parameters in UNIX / Linux '''
    from optparse import OptionParser
    script_usage = "%prog [-g] gff.file [-f] fasta.file [-o] out.fa [-i] seq_item[gene,mRNA,CDS]"
    script_description = "use this script to according to the gene-id to find the corresponding sequences from fasta.file base on the position and antisense/positive-strand descripted in gtf.file."
    optpar = OptionParser(usage=script_usage, description=script_description)
    optpar.add_option('-g', '--gff', dest='gff_file',
                      help='input the anotition-file(filename.gff)')
    optpar.add_option('-f', '--genome', dest='fasta_seq',
                      help='input the genome-fasta that comtained the whole sequences')
    optpar.add_option('-o', '--out', dest='out_fasta',
                      help='please define a filename to place your result seq')
    optpar.add_option('-i', '--item', dest='item_name',
                      help='input the target seq(gene,mRNA,CDS)you want to gain')
    options,args = optpar.parse_args()
    gff = options.gff_file
    genome_fa = options.fasta_seq
    item = options.item_name
    out_file = options.out_fasta
    Get_infomation_from_gff(gff, genome_fa, item, out_file)

if __name__ == '__main__':
    main()
