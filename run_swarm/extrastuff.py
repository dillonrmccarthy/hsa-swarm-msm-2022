#!/usr/bin/python3
from time import sleep
import random
import numpy as np

class Bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKCYAN = '\033[96m'
    OKCYANB = '\u001b[44;1m'
    b = '\033[1;4m'
    u = '\033[4m'
    ur = '\033[24m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'

def printdots(text, delay=.2):
    print(end=text)
    n_dots = 0
    #iii = 0

    while True:
        #if iii == 6:
        #    break
        #if n_dots == 3:
        if n_dots == 3:
            print(end='\b\b\b', flush=True)
            print(end='   ',    flush=True)
            print(end='\b\b\b', flush=True)
            break
            #n_dots = 0
        else:
            print(end='.', flush=True)
            n_dots += 1
        sleep(delay)
        #iii += 1


def gen_seeds(gen):
    gen = int(gen)
    if gen % 2 == 0:
        seeds = [random.randint(1111,4999) for _ in range(5)]
    else:
        seeds = [random.randint(5000,9999) for _ in range(5)]
    return seeds


if __name__ == '__main__':
    np.random.seed(1324134)
    x = np.random.randint(255, size=(5, 3))
    #print(x)
    ind = [0,2,4]
    x = x[ind]
    #print(x)



    x = np.random.randint(5, size=(6,))
    print(x)
    x = x-1
    print(x)
