# #!/usr/bin/env python3
# -*- encoding:utf-8 -*-

## 将转换生成的fa文件保存为文件
output_file = open(r'D:\Result\cigar-test-1.fa','w')  #创建一个.fa只写模式的文件作为输出文件
with open(r'D:\Result\cigar-test-1.fq') as input_fastq:
    for index,line in enumerate(input_fastq):
        if index % 4 == 0 :
            output_file.write('>'+line.strip()[1:]+'\n')
            # 将print换为.write，但是print是自动换行的，但这个不能，因此要加换行符！！
        elif index % 4 == 1:
            for i in range(0,len(line.strip()),40):
                output_file.write(line.strip()[i:i+40]+'\n')
                # 同样加换行符！
        else:
            continue
