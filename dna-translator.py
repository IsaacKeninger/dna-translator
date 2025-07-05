# Isaac Keninger
# CS 65
# Final Project

aminoAcids = {
    "UUU": "Phenylalanine", "UUC": "Phenylalanine", "UUA": "Leucine", "UUG": "Leucine",
    "CUU": "Leucine", "CUC": "Leucine", "CUA": "Leucine", "CUG": "Leucine",
    "AUU": "Isoleucine", "AUC": "Isoleucine", "AUA": "Isoleucine", "AUG": "Methionine",
    "GUU": "Valine", "GUC": "Valine", "GUA": "Valine", "GUG": "Valine",
    "UCU": "Serine", "UCC": "Serine", "UCA": "Serine", "UCG": "Serine",
    "CCU": "Proline", "CCC": "Proline", "CCA": "Proline", "CCG": "Proline",
    "ACU": "Threonine", "ACC": "Threonine", "ACA": "Threonine", "ACG": "Threonine",
    "GCU": "Alanine", "GCC": "Alanine", "GCA": "Alanine", "GCG": "Alanine",
    "UAU": "Tyrosine", "UAC": "Tyrosine", "UAA": "STOP", "UAG": "STOP",
    "CAU": "Histidine", "CAC": "Histidine", "CAA": "Glutamine", "CAG": "Glutamine",
    "AAU": "Asparagine", "AAC": "Asparagine", "AAA": "Lysine", "AAG": "Lysine",
    "GAU": "Aspartic Acid", "GAC": "Aspartic Acid", "GAA": "Glutamic Acid", "GAG": "Glutamic Acid",
    "UGU": "Cysteine", "UGC": "Cysteine", "UGA": "STOP", "UGG": "Tryptophan",
    "CGU": "Arginine", "CGC": "Arginine", "CGA": "Arginine", "CGG": "Arginine",
    "AGU": "Serine", "AGC": "Serine", "AGA": "Arginine", "AGG": "Arginine",
    "GGU": "Glycine", "GGC": "Glycine", "GGA": "Glycine", "GGG": "Glycine"
    }

aminoAcidsAbv = {
    "UUU": "Phe", "UUC": "Phe", "UUA": "Leu", "UUG": "Leu",
    "CUU": "Leu", "CUC": "Leu", "CUA": "Leu", "CUG": "Leu",
    "AUU": "Ile", "AUC": "Ile", "AUA": "Ile", "AUG": "Met",
    "GUU": "Val", "GUC": "Val", "GUA": "Val", "GUG": "Val",
    "UCU": "Ser", "UCC": "Ser", "UCA": "Ser", "UCG": "Ser",
    "CCU": "Pro", "CCC": "Pro", "CCA": "Pro", "CCG": "Pro",
    "ACU": "Thr", "ACC": "Thr", "ACA": "Thr", "ACG": "Thr",
    "GCU": "Ala", "GCC": "Ala", "GCA": "Ala", "GCG": "Ala",
    "UAU": "Tyr", "UAC": "Tyr", "UAA": "STOP", "UAG": "STOP",
    "CAU": "His", "CAC": "His", "CAA": "Glu", "CAG": "Glu",
    "AAU": "Asn", "AAC": "Asn", "AAA": "Lys", "AAG": "Lys",
    "GAU": "Asp", "GAC": "Asp", "GAA": "Glu", "GAG": "Glu",
    "UGU": "Cys", "UGC": "Cys", "UGA": "STOP", "UGG": "Trp",
    "CGU": "Arg", "CGC": "Arg", "CGA": "Arg", "CGG": "Arg",
    "AGU": "Ser", "AGC": "Ser", "AGA": "Arg", "AGG": "Arg",
    "GGU": "Gly", "GGC": "Gly", "GGA": "Gly", "GGG": "Gly"}

import random

def generateRandom_DNA(length = 10):
    """
    This generates a randomly generated strand of DNA
    Parameters: Length, (An integer), the length of the randomly generated DNA strand
    return: A list of DNA
    """
    nucleotides = ["A", "C", "T", "G"]
    DNA_list = []
    for i in range(length):
        num = (random.randint(0,3))
        DNA_list.append(nucleotides[num])
    return DNA_list

# : In transcription, A on DNA pairs with U on mRNA, C pairs with G, G pairs with C, and T on DNA pairs with A on mRNA. ''

def convert_to_RNA(sequence):
    """
    This translates a given DNA strand into RNA
    Parameters: Sequence, (A list), A dna strand
    return: A list of RNA 
    """
    rna = []
    for character in sequence:
        if character == "A":
            rna.append("U")
        elif character == "C":
            rna.append("G")
        elif character == "G":
            rna.append("C")
        elif character == "T":
            rna.append("A")
    
    return rna

def aminoSort(seq):
    """
    Translate an rna strand into amino acids
    Paramaters: Seq, (A list), an rna strand
    Returns: A list of Amino acids
    """
    aminos = []
    count1 = 0
    count2 = 3
    for i in range(len(seq)):
        
        # Join tool learned from https://www.geeksforgeeks.org/python-remove-empty-strings-from-list-of-strings/
        aminos.append(''.join(seq[count1:count2]))
        count1 += 3
        count2 += 3
        clean_aminos = list(filter(None,aminos)) # CITE THIS
        
        # This code removes any remainiding amino acids that cannot be matched
        if len(clean_aminos[-1]) != 3:
            clean_aminos = clean_aminos[:-1]

    return clean_aminos

def aminoNameSort(seq):
    """
    This translates the rna amino acid groups into their amino acid names
    Parameters: seq, (a list), a list of rna grouped into amino acids
    returns: A list of amino acid names
    """
    names = []
    for value in seq:
        if value in aminoAcids:
            names.append(aminoAcids[value])
    return names

# I needed to abbreviate the amino names so they would appear cleanly on the bar graph
def aminoNameSortAbv(seq):
    """
    This translates the rna amino acid groups into their amino acid names in abbreviated form
    Parameters: seq, (a list), a list of rna grouped into amino acids
    returns: A list of amino acid names abbreviated
    """
    names = []
    for value in seq:
        if value in aminoAcidsAbv:
            names.append(aminoAcidsAbv[value])
    return names

# Finds how many of each amino acid is present
def findOccurences(seq):
    """
    Finds how many of each amino acid is found in the sequence
    Parameters: seq, a list, sequence of amino acids
    returns: a dictionary of amino acid keys with the # of occurences in the sequence as values
    """
    occurences = {}
    for value in seq:
        if value not in occurences:
            occurences[value] = 1
        else:
            occurences[value] += 1
    return occurences

# I used 50 as an argument for the amount of nucleotides but this is customizable
# the default parameter is set to 10
DNA_sequence = generateRandom_DNA(50)
print("DNA STRAND", ''.join(DNA_sequence))

RNA_sequence = convert_to_RNA(DNA_sequence)
print("RNA STRAND:", ''.join(RNA_sequence))

aminos = aminoSort(RNA_sequence)
print("AMINOS LIST: ", aminos)    
    
aminoNames = aminoNameSort(aminos)
print("AMINO NAMES LIST: ", aminoNames)

aminoOccurences = findOccurences(aminoNames)
print("AMINO OCCURENCES: ", aminoOccurences)


aminoNamesAbv = aminoNameSortAbv(aminos)
aminoOccurencesAbv = findOccurences(aminoNamesAbv)



# Here is a visualization of the DNA sequence and Amino acid occurences
# I used Chat GPT to create these graphics.
# I replaced the data Chat-gpt created with the data used prior.

import matplotlib.pyplot as plt
import tkinter as tk
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import matplotlib.pyplot as plt

def create_visualizations(dna_sequence, amino_occurences):
    # Create the main Tkinter window
    root = tk.Tk()
    root.title("DNA and Amino Acid Visualizations")

    # Create the figure for the visualizations
    fig, axs = plt.subplots(2, 1, figsize=(15, 12))

    # DNA Visualization (Top Plot)
    colors = ['blue', 'red', 'green', 'yellow']  # Colors for A, T, C, G
    color_map = {"A": "blue", "T": "red", "C": "green", "G": "yellow"}
    dna_colors = [color_map[nucleotide] for nucleotide in dna_sequence]
    axs[0].bar(range(len(dna_sequence)), [1] * len(dna_sequence), color=dna_colors)
    axs[0].set_title("DNA Visualization")
    axs[0].set_ylabel("Base Pair")
    axs[0].set_xticks([])
    axs[0].set_yticks([])

    # Bar Graph for Amino Acid Occurrences (Bottom Plot)
    amino_acids = list(amino_occurrences.keys())
    occurrences = list(amino_occurrences.values())
    axs[1].bar(amino_acids, occurrences, color='green', width=0.5)
    axs[1].set_title("Amino Acid Occurrences")
    axs[1].set_xlabel("Amino Acids")
    axs[1].set_ylabel("Occurrences")

    # Embed the Matplotlib figure in the Tkinter window
    canvas = FigureCanvasTkAgg(fig, master=root)
    canvas_widget = canvas.get_tk_widget()
    canvas_widget.pack()

    # Start the Tkinter event loop
    root.mainloop()

# Example DNA sequence and amino acid occurrences
# Call the function to display visualizations
create_visualizations(DNA_sequence, aminoOccurencesAbv)
