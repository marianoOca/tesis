
from misc_utils import read_list_from_file, make_name
import matplotlib as mpl
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

def get_norm(scale:int, sizes:list, daset1:list, daset2:list, bins_x:int, bins_y:int, range:list):
    # necesario para que ambos gráficos sean consistentes con los colores y el colorbar, de forma din'amica
    counts1, _, _ = np.histogram2d(sizes, daset1, bins=(bins_x, bins_y), range=range)
    counts2, _, _ = np.histogram2d(sizes, daset2, bins=(bins_x, bins_y), range=range)
    vmax = max(counts1.max(), counts2.max())
    vmax = int(-(-vmax // scale)) * scale
    norm = mpl.colors.Normalize(vmin=0, vmax=vmax)
    return norm

def get_parameters_for(selector:int, dataset_name:str):
    if   selector < 4:
        if   selector == 1:
            complexity = "Icalc"
            prefix = "icalc_"
            y_range = [0.05, 0.75]
        elif selector == 2:
            complexity = "Discrepancia"
            prefix = "discr_"
            y_range = [0, 450]
        elif selector == 3:
            complexity = "Discrepancia en Bloque 2"
            prefix = "disc2_"
            y_range = [0, 450]
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
            decom_data = read_list_from_file(random_result)
            random_results.append([row[selector - 4] for row in decom_data])

        if selector == 4:
            complexity = "Kolmogorov"
            y_range = [0, 7000]
        elif selector == 5:
            complexity = "Bennett"
            y_range = [0, 17000]
        elif selector == 6:
            complexity = "Shannon Entropy"
            y_range = [0, 4.5]
        elif selector == 7:
            complexity = "2nd Ord. Entropy"
            y_range = [0, 9]
        elif selector == 8:
            complexity = "Compression Length"
            y_range = [0, 4500]
    
    return [complexity, shuffled_results, random_results, y_range]