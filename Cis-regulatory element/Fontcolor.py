#!/usr/bin/env python3
# -*- encoding:utf-8 -*-
# @Author : Haoran Pan
# date: 2020/12/19

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


def Log_output(input_str):
    '''
    这里采用cyan色进行log的输出
    :param input_str:
    :return:
    '''
    color = Colored()
    now = datetime.now().replace(microsecond=0)
    print(color.cyan(str("[{}]".format(now)) + '\t' + str(input_str)))


def Error_commond(input_str):
    '''
    这里采用magenta色进行Error command的输出
    :return:
    '''
    color = Colored()
    now = datetime.now().replace(microsecond=0)
    # print(color.magenta(str(now) + '\t' + str(input_str)))
    print(color.magenta("[{0}] Error: Error occurred when running the following command: {1} ".format(now,input_str)))


def Error_output(input_str):
    "采用magenta色用于对捕获错误进行输出"
    color = Colored()
    now = datetime.now().replace(microsecond=0)
    print(color.magenta("[{0}] Error: {1}".format(now,input_str)))


def Tips_output(input_str):
    '''
    这里采用yellow色对提示信息进行输出
    :param input_str:
    :return:
    '''
    color = Colored()
    now = datetime.now().replace(microsecond=0)
    print(color.yellow(str("[{}]".format(now)) + '\t' + str(input_str)))


def Scripts_tip(input_str):
    '''
    这里采用green色对脚本的描述进行输出
    :param input_str:
    :return:
    '''
    color = Colored()
    print(color.green(str("[Description:]")+' '+str(input_str)))