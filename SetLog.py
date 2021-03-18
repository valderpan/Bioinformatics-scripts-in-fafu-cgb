#!/usr/bin/env python3
# -*- encoding:utf-8 -*-
# @Author : Haoran Pan

import sys
from loguru import logger

logger.remove()
logger.add(sys.stdout,
           format="<green>{time:[YYYY-MM-DD HH:mm:ss]}</green>  | <level>{level}</level>  | <level>{message}</level>",
           colorize=True)

def debug_out(input_str):
    logger.debug(input_str)

def info_out(input_str):
    logger.info(input_str)

def warning_out(input_str):
    logger.warning(input_str)

def error_out(input_str):
    logger.error(input_str)

def critical_out(input_str):
    logger.critical(input_str)
