#!/usr/bin/env python3
"""
print information during searching the descent direction
"""

def print_info(fnm, str1):
    fp = open(fnm,'a')
    fp.write('%s\n'% str1)
    fp.close()