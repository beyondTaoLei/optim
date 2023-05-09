#!/usr/bin/env python3
"""
useful functions
Usage:
    [from util import fun] in another *.py file
https://stackoverflow.com/questions/20309456/call-a-function-from-another-file
"""
import os, shutil
import numpy as np
from numpy import linalg as LA

def BLOCK_LOW(rank,size,n):
    """
    rank: processor index
    size: number of processors
    n: the number of tasks
    return the index, from 0, not 1!
    """
    return rank*n//size

def BLOCK_HIGH(rank,size,n):
    return BLOCK_LOW(rank+1,size,n)-1

def BLOCK_SIZE(rank,size,n):
    return BLOCK_HIGH(rank,size,n)-BLOCK_LOW(rank,size,n)+1

