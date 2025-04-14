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
        for line in file:         #".txt" case
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


### PLOTING ###
import numpy as np

def get_boxplot_lines(data:list) -> list:
    median = np.median(data)
    upper_quartile = np.percentile(data, 75)
    lower_quartile = np.percentile(data, 25)

    iqr = upper_quartile - lower_quartile
    upper_whisker = upper_quartile+(1.5*iqr)
    lower_whisker = lower_quartile-(1.5*iqr)

    upper_whisker = max([e for e in data if e < upper_whisker])
    lower_whisker = min([e for e in data if e > lower_whisker])

    return [upper_whisker, upper_quartile, median, lower_quartile, lower_whisker]

def get_parameters_for(selector:int, dataset_name:str):
    if   selector < 4:
        if   selector == 1:
            complexity = "Icalc"
            prefix = "icalc_"
            x_range = [0.05, 0.75]
        elif selector == 2:
            complexity = "Discrepancia"
            prefix = "discr_"
            x_range = [0, 450]
        elif selector == 3:
            complexity = "Discrepancia en Bloque 2"
            prefix = "disc2_"
            x_range = [0, 450]
        source = "results/" + prefix + dataset_name

        original = read_list_from_file(source + ".txt")
        shuffled_results = [original]
        random_results = [original]

        for i in range(1,11):
            shuffled_result = make_name(source + "_s", i, ".txt")
            random_result   = make_name(source + "_r", i, ".txt")

            shuffled_results.append(read_list_from_file(shuffled_result))
            random_results.append(read_list_from_file(random_result))
    else:
        decom_data = read_list_from_file("results/decom_usp_f.csv")

        original = [row[selector - 4] for row in decom_data]
        shuffled_results = [original]
        random_results = [original]
    
        for i in range(1,11):
            shuffled_result = make_name("results/decom_usp_f_s", i, ".csv")
            random_result   = make_name("results/decom_usp_f_r", i, ".csv")

            decom_data = read_list_from_file(shuffled_result)
            shuffled_results.append([row[selector - 4] for row in decom_data])
            #decom_data = read_list_from_file(random_result)
            random_results.append([row[selector - 4] for row in decom_data])

        if selector == 4:
            complexity = "Kolmogorov"
            x_range = [0, 7000]
        elif selector == 5:
            complexity = "Bennett"
            x_range = [0, 17000]
        elif selector == 6:
            complexity = "Shannon Entropy"
            x_range = [0, 4.5]
        elif selector == 7:
            complexity = "2nd Order Entropy"
            x_range = [0, 9]
        elif selector == 8:
            complexity = "Compression length (using gzip)"
            x_range = [0, 4500]
    
    return [complexity, shuffled_results, random_results, x_range]