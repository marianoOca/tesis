import random
from miscellaneous import *


### COUNT ENTRIES ###

def count_entries(in_file):
    res = []
    map_bio(in_file, lambda _, res : res.append(1), res)
    return len(res)


### SHOW ENTRY ###

def show_entry(seq_record:SeqIO.SeqRecord):
    print(seq_record.description)
    print(repr(seq_record.seq))
    print("Size: ", len(seq_record))


### RESIZE ###

def copy_with_size(in_file, out_file, size:int):
    file = open(out_file, "w")
    for seq_record in SeqIO.parse(in_file, "fasta"):
        file.write(">" + seq_record.description + "\n")
        write_with_size(str(seq_record.seq), file, size)
    file.close()


### FILTER ###

def filter_to_file(in_file, out_file, criteria): #ej: citeria = lambda seq_record : len(seq_record) >= 50
    file = open(out_file, "w")
    for seq_record in SeqIO.parse(in_file, "fasta"):
        if (criteria(seq_record)):
            file.write(">" + seq_record.description + "\n")
            write_with_size(str(seq_record.seq), file)
    file.close()


### SORT ###

def save_to_list(seq_record:SeqIO.SeqRecord, seq_list:list):
    seq_list.append(seq_record)

def sort_to_file(in_file, out_file):
    seq_list = []
    map_bio(in_file, save_to_list, seq_list)
    seq_list.sort(key = len) #ordena según su largo
    
    file = open(out_file, "w")
    for seq_record in seq_list:
        file.write(">" + seq_record.description + "\n")
        write_with_size(str(seq_record.seq), file)
    file.close()


### SAVE SIZES ###

def size_to_list(in_file) -> list:
    print("\nCalculating sizes")
    size_list = []
    map_bio(in_file, lambda seq_record, res : res.append(len(seq_record)), size_list)
    return size_list


### SHUFFLE ###

def shuffle_seq(seq_record:SeqIO.SeqRecord, out_file):
    out_file.write(">" + seq_record.description + "\n")
    seq_list = list(str(seq_record.seq))
    random.shuffle(seq_list)
    seq_ran = ''.join(seq_list)
    write_with_size(seq_ran, out_file)

def shuffle_to_file(in_file, out_file, seed = 1):
    random.seed(seed)
    file = open(out_file, "w")
    map_bio(in_file, shuffle_seq, file)
    file.close()


### RANDOM ###

# Fuente: https://es.wikipedia.org/wiki/Formato_FASTA#Representación_de_la_secuencia
def random_seq(seq_record:SeqIO.SeqRecord, out_file):
    out_file.write(">" + seq_record.description + "\n")
    size = len(seq_record)
    #alphabet = list("TGAC")                  #para ácidos nucléicos
    alphabet = list("ACDEFGHIKLMNPQRSTVWY")   #para aminoácidos
    random_list = [random.choice(alphabet) for _ in range(size)]
    random_seq = ''.join(random_list)
    write_with_size(random_seq, out_file)

def random_to_file(in_file, out_file, seed = 1):
    random.seed(seed)
    file = open(out_file, "w")
    map_bio(in_file, random_seq, file)
    file.close()


### CONTROL DATASET ###
def single_char_seq(seq_record:SeqIO.SeqRecord, out_file):
    out_file.write(">" + seq_record.description + "\n")
    A_list = "A" * len(seq_record)
    write_with_size(A_list, out_file)

def single_char_to_file(in_file, out_file):
    file = open(out_file, "w")
    map_bio(in_file, single_char_seq, file)
    file.close()

def sorted_seq(seq_record:SeqIO.SeqRecord, out_file):
    out_file.write(">" + seq_record.description + "\n")
    seq_list = list(str(seq_record.seq))
    seq_list.sort()
    seq_ran = ''.join(seq_list)
    write_with_size(seq_ran, out_file)

def sorted_to_file(in_file, out_file):
    file = open(out_file, "w")
    map_bio(in_file, sorted_seq, file)
    file.close()