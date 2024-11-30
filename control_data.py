import random
from miscellaneous import *


### SHOW ENTRY ###

def show_entry(seq_record:SeqIO.SeqRecord):
    print(seq_record.description)
    print(repr(seq_record.seq))
    print("Size: ", len(seq_record))


### RESIZE ###

def copy_with_size(dataset, out_file, size:int):
    f = open(out_file, "w")
    for seq_record in SeqIO.parse(dataset, "fasta"):
        f.write(">" + seq_record.description + "\n")
        write_with_size(str(seq_record.seq), f, size)
    f.close()


### SORT ###

def save_to_list(seq_record:SeqIO.SeqRecord, seq_list:list):
    seq_list.append(seq_record)

def sort_to_file(dataset, out_file):
    seq_list = []
    map_bio(dataset, save_to_list, seq_list)
    seq_list.sort(key = len) #ordena según su largo
    
    f = open(out_file, "w")
    for seq_record in seq_list:
        f.write(">" + seq_record.description + "\n")
        write_with_size(str(seq_record.seq), f)
    f.close()


### SHUFFLE ###

def shuffle_seq(seq_record:SeqIO.SeqRecord, out_file):
    out_file.write(">" + seq_record.description + "\n")
    seq_list = list(str(seq_record.seq))
    random.shuffle(seq_list)
    seq_ran = ''.join(seq_list)
    write_with_size(seq_ran, out_file)

def shuffle_to_file(dataset, out_file, seed = 1):
    random.seed(seed)
    f = open(out_file, "w")
    map_bio(dataset, shuffle_seq, f)
    f.close()


### RANDOM ###

# Fuente: https://es.wikipedia.org/wiki/Formato_FASTA#Representación_de_la_secuencia
def random_seq(seq_record:SeqIO.SeqRecord, out_file):
    out_file.write(">" + seq_record.description + "\n")
    size = len(seq_record)
    #sample = list("TGAC")                  #para ácidos nucléicos
    sample = list("ACDEFGHIKLMNPQRSTVWY")   #para aminoácidos
    random_list = [random.choice(sample) for _ in range(size)]
    random_seq = ''.join(random_list)
    write_with_size(random_seq, out_file)

def random_to_file(dataset, out_file, seed = 1):
    random.seed(seed)
    f = open(out_file, "w")
    map_bio(dataset, random_seq, f)
    f.close()
