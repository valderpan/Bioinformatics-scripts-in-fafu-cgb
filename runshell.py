#!/usr/bin/env python3
# -*- encoding:utf-8 -*-
# @Author : Haoran Pan
# date: 2021/1/5

import Fontcolor as Font
import subprocess

def run_command(cmd):
    Font.Log_output(cmd)
    return_code = subprocess.call(cmd, shell=True)
    if return_code != 0:
        Font.Error_commond(cmd)
