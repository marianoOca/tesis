# Requeriments: pip install Bio

from Bio import SeqIO
import multiprocessing as mp
import random


### AUX ###
def save_list_to_file(l:list, out_file):
    file = open(out_file, 'w')
    for n in l:
        file.write(f"{n}\n")
    file.close()

def read_list_from_file(in_file) -> list:
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

def size_to_list(dataset) -> list:
    print("\nCalculating sizes")
    size_list = []
    map_bio(dataset, lambda seq_record, res : res.append(len(seq_record)), size_list)
    return size_list

def filter_to_file(dataset, out_file):
    f = open(out_file, "w")
    for seq_record in SeqIO.parse(dataset, "fasta"):
        if (len(seq_record) >= 50):
            f.write(">" + seq_record.description + "\n")
            write_with_size(str(seq_record.seq), f)
    f.close()

def make_name(prefix:str, num:int, sufix:str = "") -> str:
    if num < 10:
        middle = "0" + str(num)
    else:
        middle = str(num)
    return prefix + middle + sufix

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

def show_icalc(seq_record:SeqIO.SeqRecord):
    print(seq_record.description)
    print("Icalc: ", icalc(str(seq_record.seq)))


### DISCREPANCY ###

def Kadane_for_2(seq:str, pos:str, neg:str) -> int:
    res = 0
    maxEnding = 0

    for i in range(len(seq)):
        to_add = 1 if seq[i] == pos else (-1 if seq[i] == neg else 0)
        maxEnding = max(maxEnding + to_add, to_add)

        res = max(res, maxEnding)

    return res

def discrepancy(seq:str) -> int:
    alphabet = {"A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y"}
    res = 0

    for i in alphabet:
        remaining_alphabet = alphabet - {i}
        for j in remaining_alphabet:
            res = max(res, Kadane_for_2(seq, i, j))

    return res


## EXPERIMENT ## SORT OF CURRY ##

def mapeable_to_file(l:list): #[exp:srt, ori:str, dest:str, seed:int]
    exp = l[0]
    origin_file = l[1]
    destination_file = l[2]
    seed = l[3]
    if   exp == "s":
        shuffle_to_file(origin_file, destination_file, seed)
    elif exp == "r":
        random_to_file(origin_file, destination_file, seed)

class Info:
    def __init__(self, complexity:str):
        if   complexity == "i":
            self.prefix = "icalc_"
            self.name = "icalc"
            self.function = icalc
        elif complexity == "d":
            self.prefix = "discr_"
            self.name = "discrepancia"
            self.function = discrepancy

def complexity_to_list(dataset, complexity:str) -> list:
    res_list = []
    f = Info(complexity).function
    map_bio(dataset, lambda seq_record, res : res.append(f(str(seq_record.seq))), res_list)
    return res_list

def complexity_to_file(l:list): #[ori:str, dest:str, complexity:str]
    working_dataset = l[0]
    destination_file = l[1]
    complexity = l[2]
    res = complexity_to_list(working_dataset, complexity)
    save_list_to_file(res, destination_file)


## EXPERIMENT ## MULTIPROCESSING ##

def multiprocess(function, data:list, message:str) -> list:
    print(message)
    
    cant_processes = mp.cpu_count()
    pool = mp.Pool(processes = cant_processes)      #generamos workers como núcleos del procesador tengamos

    results = []
    for i in range(0, len(data), cant_processes):
        print("Generating " + str(i+1) + " to " + str(i+cant_processes if i+cant_processes < len(data) else len(data)) + " of "+ str(len(data)))
        chunk = data[i:i+cant_processes]            #se divide el workload en bloques de tamaño cant_processes, llamados chunks
        chunk_results = pool.map(function, chunk)   #se procesa cada chunk sincrónicamente
        results.extend(chunk_results)
    pool.close()
    pool.join()

    return results


## EXPERIMENT ##

def generate_working_files(origin_dataset, exp:str, cuantity:int):
    files_to_generate = []
    for i in range(cuantity):
        destination_file = make_name(origin_dataset + "_" + exp, i+1, ".fasta")
        files_to_generate.append([exp, origin_dataset + ".fasta", destination_file, i+1])

    multiprocess(mapeable_to_file, files_to_generate, "\nGenerating " + ("shuffled" if exp == "s" else "random") + " files:")

#cuantity = 0: se está trabajando sobre el archivo original
def calculate_complexity_from_files(origin_dataset, complexity:str, cuantity:int = 0):
    files_to_process = []
    info = Info(complexity)

    if cuantity == 0:
        print("\nCalculating " + info.name + " from orginial dataset")
        complexity_to_file([origin_dataset + ".fasta", "results/" + info.prefix + origin_dataset + ".txt", complexity])
    else:
        for i in range(cuantity):
            working_dataset = make_name(origin_dataset, i+1, ".fasta")
            destination_file = make_name(info.prefix + origin_dataset, i+1, ".txt")
            files_to_process.append([working_dataset, "results/" + destination_file, complexity])

        multiprocess(complexity_to_file, files_to_process, "\nCalculating " + info.name + " for " + str(cuantity) + " files:")


#exp = "s": shuffle y exp = "r":random
#complexity = "i":icalc y complexity = "d":discrepancia
#gen indica si se debe generar los datos de prueba random y/o shuffled
def experiment(dataset, complexity:str, exp:str = "s_and_r", gen:bool = False, cuantity:int = 10):
    if exp != "s_and_r":
        if gen:
            generate_working_files(dataset, exp, cuantity)
        calculate_complexity_from_files(dataset + "_" + exp, complexity, cuantity)
    else:
        if gen:
            sizes = size_to_list(dataset + ".fasta")
            save_list_to_file(sizes, "sizes_" + dataset + ".txt")
            generate_working_files(dataset, "s", cuantity)
            generate_working_files(dataset, "r", cuantity)

        calculate_complexity_from_files(dataset, complexity)
        calculate_complexity_from_files(dataset + "_s", complexity, cuantity)
        calculate_complexity_from_files(dataset + "_r", complexity, cuantity)
