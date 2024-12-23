from Bio import SeqIO


### AUX ###

def write_with_size(seq:str, out_file, size:int = 100):
    i = size
    while i < len(seq):
        out_file.write(seq[i-size:i] + "\n")
        i = i + size
    out_file.write(seq[i-size:] + "\n\n")

def make_name(prefix:str, num:int, sufix:str = "") -> str:
    if num < 10:
        middle = "0" + str(num)
    else:
        middle = str(num)
    return prefix + middle + sufix


## LIST <--> FILE ##

def save_list_to_file(l:list, out_file):
    file = open(out_file, 'w')
    
    if type(l[0]) != list:
        for n in l:
            file.write(f"{n}\n")
    else:
        for row in l:
            file.write(",".join(map(str, row)) + "\n")
    file.close()

def read_list_from_file(in_file) -> list:
    res = []
    file =  open(in_file, 'r')
    if in_file[len(in_file)-4:] == ".csv":
        for line in file:
            row_srt = line.strip().split(",")
            row_floats = [float(n) for n in row_srt]
            row_num = [int(n) if n.is_integer() else n for n in row_floats]
            res.append(row_num)
    else:
        for line in file:
            num = float(line.strip())
            if num.is_integer():
                num = int(num)
            res.append(num)
    file.close()
    return res

### MAP DATASET ###

def map_bio(in_file, function, auxVar = False):
    if auxVar == False:
        for seq_record in SeqIO.parse(in_file, "fasta"):
            function(seq_record)
    else:
        for seq_record in SeqIO.parse(in_file, "fasta"):
            function(seq_record, auxVar)