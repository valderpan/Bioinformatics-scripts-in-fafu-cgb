脚本使用：（这里以chr9为例子进行演示）

1.将保留下来的10个group与剩余所有group进行划分，保留下来的10个group称为conserverd_group,剩余group成为nonconserverd_group,按照20一组保存到excel文件中
$ python Roc22_chrCorrect.py -g g9.xlsx -q 6,7,8,13,31,49,62,70,73,103
参数说明：
-g 输入chr文件(例如g9.xlsx)
-q 保留下来的group 的id

2.生成shell脚本来移动每个group的blocks
$ python Roc22_chrCorrect.py -x1 nonconserved_group_1-20.xlsx -x2 nonconserved_group_21-40.xlsx -x3 nonconserved_group_41-60.xlsx -x4 nonconserved_group_61-80.xlsx -x5 nonconserved_group_81-end.xlsx -O run_script.sh
参数说明：
-x1,-x2,-x3,-x4,-x5分别为生成的nonconserverd_group文件
-O 生成shell脚本

3.运行脚本
$ nohup bash run_script.sh &

4.结果解读
生成的最后的xlsx文件即为最终结果文件，chr9最终文件为90.xlsx
remain*.xlsx 文件为每个nonconserverd_group中的未分配到的片段


注意：此脚本目前不适用于chr8染色体(断裂)！！！