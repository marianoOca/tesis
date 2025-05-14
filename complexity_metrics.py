import rpy2.robjects as ro
import re

#from rpy2.robjects import pandas2ri
#from miscellaneous import *


class ComplexitySelector:
    def __init__(self, c_id:str):
        self.extension = ".txt" #funciones que devuelven un sólo parámetro
        if   c_id == "i":
            self.prefix = "icalc_"
            self.name = "icalc"
            self.function = icalc
        elif c_id == "d":
            self.prefix = "discr_"
            self.name = "discrepancia"
            self.function = discrepancy
        elif c_id == "d2":
            self.prefix = "disc2_"
            self.name = "discrepancia en bloque de 2"
            self.function = lambda seq : discrepancy(seq, 2)
        elif c_id == "d3":
            self.prefix = "disc3_"
            self.name = "discrepancia en bloque de 3"
            self.function = lambda seq : discrepancy(seq, 3)
        elif c_id == "d4":
            self.prefix = "disc4_"
            self.name = "discrepancia en bloque de 4"
            self.function = lambda seq : discrepancy(seq, 4)
        elif c_id == "b":
            self.prefix = "decom_"
            self.name = "Block Decomposition Method"
            self.function = bdm
            self.extension = ".csv" #funciones que devuelven más de un parámetro (como list)
        else:
            raise ValueError("Unknown complexity id.")


### ICALC ###

def icalc(seq:str) -> float:
    b = []
    for i in range(len(seq)):
        mr = 0
        for j in range(i, mr, -1):
            if seq[j-mr:j] != seq[i-mr+1:i+1]:
                continue
            while i - mr >= 0 and seq[i-mr] == seq[j-mr-1]:
                mr += 1
        b.append(mr)
    
    res = 0
    for i in range(len(b)):
        res += 1.0 / (1.0 + b[i])
    return res / len(b)

#from Bio import SeqIO
#def show_icalc(seq_record:SeqIO.SeqRecord):
#    print(seq_record.description)
#    print("Icalc: ", icalc(str(seq_record.seq)))


### DISCREPANCY ###

def Kadane_for_2blocks(seq:str, pos:str, neg:str) -> int:
    res = 0
    maxEnding = 0
    pos_index = 0
    neg_index = 0

    for i in range(len(seq)-len(pos)+1):
        if seq[i:i+len(pos)] == pos and i >= pos_index:
            to_add = 1
            pos_index = i + len(pos)
        elif seq[i:i+len(pos)] == neg and i >= neg_index:
            to_add = -1
            neg_index = i + len(pos)
        else:
            to_add = 0
        maxEnding = max(maxEnding + to_add, to_add)

        res = max(res, maxEnding)

    return res

def discrepancy(seq:str, block_size:int = 1) -> int:
    #alphabet = set("TGAC")                  #para ácidos nucléicos
    alphabet = set("ACDEFGHIKLMNPQRSTVWY")   #para aminoácidos
    working_alph = alphabet.copy()
    res = 0

    for _ in range(block_size-1):
        working_alph = {i + j for i in working_alph for j in alphabet}

    for i in working_alph:
        remaining_alphabet = working_alph - {i}
        for j in remaining_alphabet:
            res = max(res, Kadane_for_2blocks(seq, i, j))

    return res


### BENNETT & Kolmogorov / Block Decomposition Method ###

def bdm(seq:str) -> list:
    ro.r['source']('OACC-master/bennett.R')     # Carga el script en R
    bennett = ro.r['bennett']                   # Accede a la funcuón de R
    result = bennett(seq)                       # la función devuelve un data.frame de R

    # Extrae la columna 'value' como lista Python del data.frame de R
    value_column = list(result.rx2('values'))
    # Extrae los valores numéricos usando expresiones regulares
    numbers = [float(re.search(r'[0-9.]+', str(item)).group()) for item in value_column]
    
    return numbers
