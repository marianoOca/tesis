from main import *

if __name__ == "__main__":

    experiment("usp_987", "d2", gen = True, quantity= 5, mode="feedback")
    #this will
    #   - generate 5 shuffled and 5 random control datasets
    #   - generate a single character and a sorted control dataset
    #   - calculate Discrepancy in block 2 for usp_987.fasta and the generated control datasets
    #     in feedback mode, this means every individual result will be file saved after computed
    #   - save the results in txt files in the results folder adding disc2_ at the beginning of 
    #     name of every processed file
    #Note: this requires data/usp_987.fasta to exists
    
## MORE TOOLKIT USE EXAMPLES ##

    #make_name("usp_f_s", 2, ".fasta")
    #Will return "usp_f_s02.fasta"

    #filter_to_file("data/uniprot_sprot.fasta", "data/usp_f.fasta" ,lambda seq_record : len(seq_record) >= 50)
    #From data/uniprot_sprot.fasta copies to data/usp_f.fasta only the entries that have length equal or longer than 50

    #sort_to_file("data/usp_f.fasta", "data/usp_f.fasta")
    #Will sort all the entries by length from data/usp_f.fasta and save them in data/usp_f.fasta

    #count_entries("data/usp_f.fasta")
    #Will return the amount of fasta entries in the file data/usp_f.fasta

    #generate_control_files("usp_987")
    #Generates the files data\singl_usp_987.fasta and data\sorte_usp_987.fasta with single character and sorted versions of usp_987
    
    #generate_working_files("usp_f", "r")
    #Will generate 10 random versions of data\usp_f.fasta, named from data\usp_f_r01.fasta to data\usp_f_r10.fasta

    #complexity_from_files("usp_f", "i")
    #Calculates Icalc complexity from data\usp_f.fasta and saves the results in results\icalc_usp_f.txt

#Deeper explanations and more examples are described in the appendix of the thesis available at github.com/marianoOca/tesis