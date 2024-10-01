# Requeriments: pip install Bio

from Bio import SeqIO
import random


### AUX ###
def save_list_to_file(l:list, out_file):
    file = open(out_file, 'w')
    for n in l:
        file.write(f"{n}\n")
    file.close()

def read_list_from_file(in_file):
    res = []
    file =  open(in_file, 'r')
    for line in file:
        num = float(line.strip())
        if num.is_integer():
            num = int(num)
        res.append(num)
    file.close()
    return res

def write_with_size(seq:str, out_file, size:int = 100):
    i = size
    while i < len(seq):
        out_file.write(seq[i-size:i] + "\n")
        i = i + size
    out_file.write(seq[i-size:] + "\n\n")

def size_to_list_aux(seq_record:SeqIO.SeqRecord, res:list):
    res.append(len(seq_record))

def size_to_list(dataset):
    size_list = []
    map_bio(dataset, size_to_list_aux, size_list)
    return size_list

def filter_to_file(dataset, out_file):
    f = open(out_file, "w")
    for seq_record in SeqIO.parse(dataset, "fasta"):
        if (len(seq_record) >= 50):
            f.write(">" + seq_record.description + "\n")
            write_with_size(str(seq_record.seq), f)
    f.close()

### MAP DATASET ###

def map_bio(dataset, function, auxVar = False):
    if auxVar == False:
        for seq_record in SeqIO.parse(dataset, "fasta"):
            function(seq_record)
    else:
        for seq_record in SeqIO.parse(dataset, "fasta"):
            function(seq_record, auxVar)


### sHOW ENTRY ###

def show_entry(seq_record:SeqIO.SeqRecord):
    print(seq_record.description)
    print(repr(seq_record.seq))
    print("Size: ", len(seq_record))


### RESIZE ###

def copy_with_size(dataset, out_file, size):
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
    seq_ran = ''.join(seq_list) #no se xq no funciona con str(seq_list)
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
    #sample = list("TGAC") #este es para ls_orchid.fasta
    sample = list("ACDEFGHIKLMNPQRSTVWY") #basado en proteinnet7.fasta, 20 aminoácidos
    random_list = [random.choice(sample) for _ in range(size)]
    random_seq = ''.join(random_list)
    write_with_size(random_seq, out_file)

def random_to_file(dataset, out_file, seed = 1):
    random.seed(seed)
    f = open(out_file, "w")
    map_bio(dataset, random_seq, f)
    f.close()


### ICALC ###

def icalc(seq:str):
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

def show_icalc(seq_record:SeqIO.SeqRecord):
    print(seq_record.description)
    print("Icalc: ", icalc(str(seq_record.seq)))

def icalc_to_list_aux(seq_record:SeqIO.SeqRecord, res:list):
    res.append(icalc(str(seq_record.seq)))

def icalc_to_list(dataset):
    icalc_list = []
    map_bio(dataset, icalc_to_list_aux, icalc_list)
    return icalc_list


## EXPERIMENT ##

def generate_working_files(origin_dataset, exp : str, cuantity : int):
    for i in range(1, cuantity + 1):
        if i < 10:
            sufix = "_" + exp + "0" + str(i) + ".fasta"
        else:
            sufix = "_" + exp + str(i) + ".fasta"

        if exp == "s":
            shuffle_to_file(origin_dataset + ".fasta", origin_dataset + sufix, i)
        if exp == "r":
            random_to_file(origin_dataset + ".fasta", origin_dataset + sufix, i)


def calculate_icalc_from_files(origin_dataset, cuantity : int):
    for i in range(1, cuantity + 1):
        if i < 10:
            working_dataset = origin_dataset + "0" + str(i)
        else:
            working_dataset = origin_dataset + str(i)
        
        icalc_list =  icalc_to_list(working_dataset + ".fasta")
        save_list_to_file(icalc_list, "icalc_" + working_dataset + ".txt")

#exp = s: shuffle y exp = r:random
def experiment(dataset, exp : str = "s_and_r", cuantity : int = 10):
    if exp != "s_and_r":
        generate_working_files(dataset, exp, cuantity)
        calculate_icalc_from_files(dataset + "_" + exp, cuantity)
    else:
        icalc_list =  icalc_to_list(dataset + ".fasta")
        save_list_to_file(icalc_list, "icalc_" + dataset + ".txt")

        generate_working_files(dataset, "s", cuantity)
        calculate_icalc_from_files(dataset + "_s", cuantity)
        generate_working_files(dataset, "r", cuantity)
        calculate_icalc_from_files(dataset + "_r", cuantity)
