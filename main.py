# Requeriments: pip install Bio
# Requeriments: pip install rpy2
# Requeriments R: install(acss)

import multiprocessing as mp

from miscellaneous import *
from control_data import *
from complexity_metrics import *


## SOME KIND OF CURRY ## NOT FOR USER ##

def handle_data_generation_from_list(l:list): #[exp:srt, ori:str, dest:str, seed:int]
    exp = l[0]
    origin_file = l[1]
    destination_file = l[2]
    seed = l[3]
    if   exp == "s":
        shuffle_to_file(origin_file, destination_file, seed)
    elif exp == "r":
        random_to_file(origin_file, destination_file, seed)
    else:
        raise ValueError("exp must be \"r\" for random or \"s\" for shuffled.")

def handle_complexity_from_list(l:list): #[ori:str, dest:str, complexity_id:str, mode:str]
    working_dataset = l[0]
    destination_file = l[1]
    complexity_id = l[2]
    mode = l[3]
    if  mode == "performance":
        res = complexity_to_list(working_dataset, complexity_id)
        save_list_to_file(res, destination_file)
    elif mode == "feedback":
        complexity_to_file_with_feedback(working_dataset, destination_file, complexity_id)
    else:
        raise ValueError("mode must be \"performance\" or \"feedback\".")


## COMPLEXITY TO LIST AND TO FILE ##

def complexity_to_list(dataset, complexity_id:str) -> list:
    res_list = []
    f = ComplexitySelector(complexity_id).function
    map_bio(dataset, lambda seq_record, res : res.append(f(str(seq_record.seq))), res_list)
    return res_list

def complexity_to_file_with_feedback(working_dataset, out_file, complexity_id):
    f = ComplexitySelector(complexity_id).function
    for seq_record in SeqIO.parse(working_dataset, "fasta"):
        res = f(str(seq_record.seq))
        file = open(out_file, 'a')
        if type(res) != list:
            file.write(f"{res}\n")
        else:
            file.write(",".join(map(str, res)) + "\n")
        file.close()


## MULTIPROCESSING ##

def multiprocess(function, data:list) -> list:
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

def generate_working_files(dataset_name, exp:str, cuantity:int):
    files_to_generate = []
    for i in range(cuantity):
        destination_file = make_name("data/" + dataset_name + "_" + exp, i+1, ".fasta")
        files_to_generate.append([exp, dataset_name + ".fasta", destination_file, i+1])

    print("\nGenerating " + ("shuffled" if exp == "s" else "random") + " files:")
    multiprocess(handle_data_generation_from_list, files_to_generate)

#cuantity = 0: se está trabajando sobre el archivo original
def complexity_from_files(dataset_name, complexity_id:str, cuantity:int = 0, mode:str = "performance"):
    files_to_process = []
    info = ComplexitySelector(complexity_id)
    if cuantity == 0:
        print("\nCalculating " + info.name + " from orginial dataset")
        handle_complexity_from_list(["data/" + dataset_name + ".fasta", "results/" + info.prefix + dataset_name + info.extension, complexity_id, mode])
    else:
        for i in range(cuantity):
            working_dataset = make_name("data/" + dataset_name, i+1, ".fasta")
            destination_file = make_name("results/" + info.prefix + dataset_name, i+1, info.extension)
            files_to_process.append([working_dataset, destination_file, complexity_id, mode])

        print("\nCalculating " + info.name + " for " + str(cuantity) + " files:")
        multiprocess(handle_complexity_from_list, files_to_process)


#exp = "s": shuffle y exp = "r":random
#complexity_id = "i":icalc y complexity_id = "d":discrepancia
#gen indica si se debe generar los datos de prueba random y/o shuffled
#mode: "performance" para mayor performance pero el resultado sólo se va a ver al final y
#      "feedback" para ir viendo los resultados a medida que se computan, pero va a llevar más tiempo en computar
def experiment(dataset_name, complexity_id:str, exp:str = "s_and_r", gen:bool = False, cuantity:int = 10, mode:str = "performance"):
    if exp != "s_and_r":
        if gen:
            generate_working_files(dataset_name, exp, cuantity)
        complexity_from_files(dataset_name + "_" + exp, complexity_id, cuantity)
    else:
        if gen:
            sizes = size_to_list(dataset_name + ".fasta")
            save_list_to_file(sizes, "sizes_" + dataset_name + ".txt")
            generate_working_files(dataset_name, "s", cuantity)
            generate_working_files(dataset_name, "r", cuantity)

        complexity_from_files(dataset_name, complexity_id, 0, mode = mode)
        complexity_from_files(dataset_name + "_s", complexity_id, cuantity, mode)
        complexity_from_files(dataset_name + "_r", complexity_id, cuantity, mode)
