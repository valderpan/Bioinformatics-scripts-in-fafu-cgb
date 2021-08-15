#!/usr/bin/env python3
# -*- encoding:utf-8 -*-
# @Author : Haoran Pan
# date: 2021/8/15

import logging
from rich.logging import RichHandler

FORMAT = "%(message)s"
logging.basicConfig(
    level="NOTSET", format=FORMAT, datefmt="[%Y-%m-%d %H:%M:%S]", handlers=[RichHandler(rich_tracebacks=True)]
)
log = logging.getLogger("rich")


def debug_out(input_str):
    log.debug(input_str)

def info_out(input_str):
    log.info(input_str)

def warning_out(input_str):
    log.warning(input_str)

def error_out(input_str):
    log.error(input_str)

def critical_out(input_str):
    log.critical(input_str)
