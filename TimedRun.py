#!/usr/bin/env python3
# -*- encoding:utf-8 -*-
# @Author : Haoran Pan
# date: 2021/1/5


import time
import Fontcolor as Font
import runshell
import sys


def TimedRunScript(command_file,sleeptime=0):
    if sleeptime:
        time.sleep(int(sleeptime))
    Font.Log_output('Total sleep time : {} seconds'.format(sleeptime))
    Font.Log_output('The shell program starts running !')
    with open(command_file) as f:
        lines = (line.strip() for line in f)
        for line in lines:
            runshell.run_command(line)
            time.sleep(1)
    Font.Log_output('The shell program has finished running !')


if __name__ == '__main__':
    if len(sys.argv) < 2 or len(sys.argv) > 3:
        Font.Scripts_tip('Run shell program after sleep time')
        print('Usage:')
        print('\tpython {0} <CommandFile> <sleeptime[default:0]>'.format(sys.argv[0]))
    elif len(sys.argv) == 2:
        TimedRunScript(sys.argv[1])
    elif len(sys.argv) == 3:
        TimedRunScript(sys.argv[1],sys.argv[2])