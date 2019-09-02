# #!/usr/bin/env python3
# -*- encoding:utf-8 -*-

'''
将fastq文件转换为fasta文件
'''
with open(r'D:\Result\cigar-test-1.fq') as input_fastq:
    for index,line in enumerate(input_fastq):
        if index % 4 == 0 :   #index从0开始,此行命令表示fastq文件的第一行
            # print(index)
            print('>'+line.strip()[1:])  #可以尝试一下删掉.strip()[1:] 其中[1：]表示不要@
        elif index % 4 == 1:
            for i in range(0,len(line.strip()),40):
                # print(i)  结果显示0 40 80 120
                # print(line.strip())
                print(line.strip()[i:i+40]) #0-40\40-80\80-120
        else:
            continue
