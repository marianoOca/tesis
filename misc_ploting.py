
from misc_utils import read_list_from_file, make_name
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from matplotlib.ticker import FuncFormatter
import matplotlib as mpl
import numpy as np

def thousands(x, pos):
    return f"{int(x):,}".replace(",", "\u202F") #también puede ser " " ó "\u200A" ó "\u202F" ó "\u2009"

def set_thousands_formatter(ax):
    formatter = FuncFormatter(thousands)

    _, xmax = ax.get_xlim()
    _, ymax = ax.get_ylim()

    if xmax >= 1000:
        ax.xaxis.set_major_formatter(formatter)
    if ymax >= 1000:
        ax.yaxis.set_major_formatter(formatter)

def set_grid(axes):
    for ax in axes:
        if isinstance(ax, mpl.axes.Axes): #esto es necesario porque axs y plt manejan la grilla de forma distinta
            current_ax = ax
        else:
            current_ax = ax.gca()
        current_ax.set_axisbelow(True)
        current_ax.grid(which='major', linestyle='-', alpha=0.4)
        current_ax.minorticks_on()
        current_ax.grid(which='minor', linestyle=':', alpha=0.3)
        set_thousands_formatter(current_ax)

def draw_boxplot(dataset:list) -> str:
    lines = get_boxplot_lines(dataset)
    return "|" +str(lines[0])[:7]+ "--[" +str(lines[1])[:7]+ "|" +str(lines[2])[:7]+ "]" +str(lines[3])[:7]+ "--|" +str(lines[4])[:7]

def print_boxplot_lines(shuffled_dataset, random_dataset = None, original_dataset = None) -> None:
    if random_dataset is None:
        print("Boxplot lines: " + draw_boxplot(shuffled_dataset))
    else:
        print("Boxplot lines Shuffled: " + draw_boxplot(shuffled_dataset))
        print("Boxplot lines Random  : " + draw_boxplot(random_dataset))
    if original_dataset is not None:
        print("Boxplot lines Original: " + draw_boxplot(original_dataset))

def get_boxplot_lines(data:list) -> list:
    median = np.median(data)
    upper_quartile = np.percentile(data, 75)
    lower_quartile = np.percentile(data, 25)

    iqr = upper_quartile - lower_quartile
    upper_whisker = upper_quartile+(1.5*iqr)
    lower_whisker = lower_quartile-(1.5*iqr)

    upper_whisker = max([e for e in data if e <= upper_whisker])
    lower_whisker = min([e for e in data if e >= lower_whisker])

    return [upper_whisker, upper_quartile, median, lower_quartile, lower_whisker]

def get_upper_limit(data:list, scale:int = 100) -> int:
    return int(-(-get_boxplot_lines(data)[0] // scale)) * scale #usa múltiplos de scale para el techo

def get_norm(scale:int, sizes:list, daset1:list, daset2:list, bins_x:int, bins_y:int, range:list):
    # necesario para que ambos gráficos sean consistentes con los colores y el colorbar, de forma din'amica
    counts1, _, _ = np.histogram2d(sizes, daset1, bins=(bins_x, bins_y), range=range)
    counts2, _, _ = np.histogram2d(sizes, daset2, bins=(bins_x, bins_y), range=range)
    vmax = max(counts1.max(), counts2.max())
    vmax = int(-(-vmax // scale)) * scale
    norm = mpl.colors.Normalize(vmin=0, vmax=vmax)
    return norm

def get_colormap() -> mpl.colors.ListedColormap:
    #Set colormap (un quilombo)
    magma = mpl.colormaps['magma'].resampled(256)
    newcolors = magma(np.linspace(0, 1, 100000))
    transparent = np.array([256/256, 256/256, 256/256, 0]) # color en formato rgb: [R/256, G/256, B/256, alpha], lo pongo en 0 para que sea transparente
    newcolors[:1, :] = transparent
    return ListedColormap(newcolors)

def get_parameters_for(selector:int, dataset_name:str) -> list:
    if   selector < 4:
        if   selector == 1:
            complexity = "Icalc"
            prefix = "icalc_"
            y_range = [0, 0.72]
        elif selector == 2:
            complexity = "Discrepancia"
            prefix = "discr_"
            y_range = [0, 450]
        elif selector == 3:
            complexity = "Discrepancia en Bloque 2"
            prefix = "disc2_"
            y_range = [0, 450]
        source = "results/" + prefix + dataset_name

        original_resuts = read_list_from_file(source + ".txt")
        shuffled_results = []
        random_results = []

        for i in range(1,11):
            shuffled_result = make_name(source + "_s", i, ".txt")
            random_result   = make_name(source + "_r", i, ".txt")

            shuffled_results.append(read_list_from_file(shuffled_result))
            random_results.append(read_list_from_file(random_result))

        single_char_results = read_list_from_file("results/" + prefix + "singl_" + dataset_name + ".txt")
        sorted_results = read_list_from_file("results/" + prefix + "sorte_" + dataset_name + ".txt")
    else:
        decom_data = read_list_from_file("results/decom_usp_f.csv")

        original_resuts = [row[selector - 4] for row in decom_data]
        shuffled_results = []
        random_results = []

        for i in range(1,11):
            shuffled_result = make_name("results/decom_usp_f_s", i, ".csv")
            random_result   = make_name("results/decom_usp_f_r", i, ".csv")

            decom_data = read_list_from_file(shuffled_result)
            shuffled_results.append([row[selector - 4] for row in decom_data])
            decom_data = read_list_from_file(random_result)
            random_results.append([row[selector - 4] for row in decom_data])

        decom_data = read_list_from_file("results/decom_singl_usp_f.csv")
        single_char_results = [row[selector - 4] for row in decom_data]
        decom_data = read_list_from_file("results/decom_sorte_usp_f.csv")
        sorted_results = [row[selector - 4] for row in decom_data]

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

    return [complexity, original_resuts, shuffled_results, random_results, single_char_results, sorted_results, y_range]