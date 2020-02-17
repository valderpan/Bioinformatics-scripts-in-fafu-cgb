# #!/usr/bin/env python3
# -*- encoding:utf-8 -*-

def read_fasta(file_path=''):
    '''
    Loading FASTA file and return a iterative object
    '''
    line = ''
    try:
        fasta_handle = open(file_path,'r')
    except:
        raise IOError('Your input FASTA is no right!')

    #make sure the file is not empty
    while True:
        line = fasta_handle.readline()
        if line[0] == '':
            return
        if line[0] == '>':
            break

    # when the file is not empty , try to load the file
    while True:
        if line[0] != '>':
            raise ValueError("Records in FASTA file must start with '>'")
        title = line[1:].rstrip()
        lines = []
        line = fasta_handle.readline()
        while True:
            if not line:
                break
            elif line[0] == '>':
                break
            lines.append(line.rstrip())
            line = fasta_handle.readline()

        yield title,''.join(lines).replace(" ","").replace("\r","")
        if not line:
            return

    fasta_handle.close()
    assert False,"Your input FASTA file have format problem."


## 计算fasta中每条chr的序列长度
SP1_genemo={}
for chr_name,chr_seq in read_fasta(file_path='D:\Result\SP1.fa'): # 可以尝试读取hg19.fa
    SP1_genemo[chr_name] =chr_seq
    # print(chr_name)
    # print(chr_seq)
    print(SP1_genemo.keys())
    # print(SP1_genemo.values())
    for chr_name in SP1_genemo.keys():
        print(chr_name,len(SP1_genemo.get(chr_name)))
