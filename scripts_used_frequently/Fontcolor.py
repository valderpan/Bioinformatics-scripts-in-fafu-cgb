#!/usr/bin/env python3
# -*- encoding:utf-8 -*-
# @Author : Haoran Pan
# date: 2020/11/22

from datetime import datetime
from colorama import init, Fore, Back, Style
init(autoreset=False)


class Colored(object):
    # 前景色:红色 背景色:默认
    def red(self, s):
        return Fore.RED + s + Fore.RESET
    # 前景色:绿色 背景色:默认
    def green(self, s):
        return Fore.GREEN + s + Fore.RESET
    # 前景色:黄色 背景色:默认
    def yellow(self, s):
        return Fore.YELLOW + s + Fore.RESET
    # 前景色:蓝色 背景色:默认
    def blue(self, s):
        return Fore.BLUE + s + Fore.RESET
    # 前景色:洋红色 背景色:默认
    def magenta(self, s):
        return Fore.MAGENTA + s + Fore.RESET
    # 前景色:青色 背景色:默认
    def cyan(self, s):
        return Fore.CYAN + s + Fore.RESET
    # 前景色:白色 背景色:默认
    def white(self, s):
        return Fore.WHITE + s + Fore.RESET
    # 前景色:黑色 背景色:默认
    def black(self, s):
        return Fore.BLACK
    # 前景色:白色 背景色:绿色
    def magenta_green(self, s):
        return Fore.MAGENTA + Back.GREEN + s + Fore.RESET + Back.RESET

# now = datetime.now().replace(microsecond=0)
# color = Colored()
# print(color.red('{} I am red!'.format(now)))
# print(color.green('{} | I am gree!'.format(now)))
# print(color.yellow('I am yellow!'))
# print(color.blue('I am blue!'))
# print(color.magenta('I am magenta!'))
# print(color.cyan('I am cyan!'))
# print(color.white('I am white!'))
# print(color.magenta_green('I am white green!'))

def Log_output(input_str):
    '''
    这里采用cyan色进行log的输出
    :param input_str:
    :return:
    '''
    color = Colored()
    now = datetime.now().replace(microsecond=0)
    print(color.cyan(str(now) + '\t' + str(input_str)))

def Error_output(input_str):
    '''
    这里采用magenta色进行Error的输出
    :return:
    '''
    color = Colored()
    now = datetime.now().replace(microsecond=0)
    # print(color.magenta(str(now) + '\t' + str(input_str)))
    print(color.magenta("Error: [{0}] Error occurred when running the following command: {1} ".format(now,input_str)))