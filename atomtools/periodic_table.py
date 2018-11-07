import pkgutil
import numpy as np

class periodic_table:
    size = 86
    vatomic_name = []
    vatomic_number = []
    vatomic_mass = []
    vatomic_outermost_nelectron = []

    # File of atom_tools/Atom.txt
    periodic_file = pkgutil.get_data('atomtools', 'Atoms.txt').decode().split("\n")
    for line in periodic_file:
        splited_line = line.split()
        vatomic_number.append(int(splited_line[0]))
        vatomic_name.append(splited_line[1])
        vatomic_mass.append(float(splited_line[3]))
        vatomic_outermost_nelectron.append(int(splited_line[2]))

        if vatomic_number[-1]==size:
            break


def name2num(name):
    """
    INPUT: atom name\n
    OUTPUT: atom number\n
    """
    # print("name ", name)
    num = periodic_table.vatomic_number[periodic_table.vatomic_name.index(name)]
    return num


def num2name(Z):
    """
    INPUT: atom number\n
    OUTPUT: atom name\n
    """
    name = periodic_table.vatomic_name[periodic_table.vatomic_number.index(Z)]
    return name


def num2group(Z):
    """
    INPUT: atom number\n
    OUTPUT: atom group in periodic table\n
    """
    group = periodic_table.vatomic_outermost_nelectron[periodic_table.vatomic_number.index(Z)]
    return group


def num2mass(Z):
    """
    INPUT: atom number\n
    OUTPUT: atomic weight\n
    """
    mass = periodic_table.vatomic_mass[periodic_table.vatomic_number.index(Z)]

    return mass

def name2mass(name):
    return periodic_table.vatomic_mass[periodic_table.vatomic_number.index(name2num(name))]
