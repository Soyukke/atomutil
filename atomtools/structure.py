import numpy as np
from numpy.linalg import norm
import os
from atomtools import periodic_table
import copy
from copy import deepcopy
import subprocess
import itertools

class structure:
    def __init__(self, natom, ntypat, typat, znucl, xcart, xred, rprim, selective=[]):
        self.natom = natom
        self.ntypat = ntypat
        self.typat = typat
        self.znucl = znucl
        self.xcart = xcart
        self.xred = xred
        self.rprim = rprim
        self.selective = selective

    def copy(self):
        """
        INPUT: None\n
        OUTPUT: Copy of object\n 
        """
        natom = deepcopy(self.natom)
        ntypat = deepcopy(self.ntypat)
        typat = deepcopy(self.typat)
        znucl = deepcopy(self.znucl)
        xcart = deepcopy(self.xcart)
        xred = deepcopy(self.xred)
        rprim = deepcopy(self.rprim)
        selective = deepcopy(self.selective)
        return structure(natom=natom, ntypat=ntypat, typat=typat, znucl=znucl, xcart=xcart, xred=xred, rprim=rprim, selective=selective)


    def remove(self, index):
        """
        INPUT: index of atom\n
        OUTPUT: index atom removed structure\n  
        """
        natom = self.natom - 1
        removed_typat = self.typat[index]
        typat = np.delete(self.typat, index)
        znucl = deepcopy(self.znucl)
        ntypat = deepcopy(self.ntypat)
        xcart = np.delete(self.xcart, index, axis=0)
        xred = np.delete(self.xred, index, axis=0)
        rprim = deepcopy(self.rprim)
        # 削除する原子と同じ種類の原子が他にないときはntypat と　znuclから消す
        if not (removed_typat in typat):
            znucl = np.delete(znucl, removed_typat-1)
            ntypat -= 1
            # 削除したらtypatも正しくする
            unique_typat = np.unique(typat)
            for i in range(len(unique_typat)):
                typat[typat==unique_typat[i]] = i + 1

        return structure(natom=natom, ntypat=ntypat, typat=typat, znucl=znucl, xcart=xcart, xred=xred, rprim=rprim)

    def rotate(self, vn, radian):
        """
        INPUT : vn: rotate axis, radian
        OUTPUT: rotated structure, rotate center is mass center
        """
        # Center to 0, 0 ,0
        p_center = gravity_center(structure_=self)
        newstruct = shift(structure_=self, vshift=-p_center)

        newxcart = rodrigues_rotation(vp=self.xcart, vn=vn, a=radian)
        newxred = cartesian2direct(rprim=self.rprim, xcart=newxcart)
        newstruct.xcart = newxcart
        newstruct.xred = newxred

        # Shift center to origin
        newstruct = shift(structure_=newstruct, vshift=p_center)
        return newstruct
        
    def volume(self):
        """
        INPUT: None\n
        OUTPUT: Volume of lattice\n
        """
        return np.abs(np.dot(self.rprim[0], np.cross(self.rprim[1], self.rprim[2])))

    def molar_mass(self):
        """
        INPUT : None\n
        OUTPUT: Total mol mass\n
        """
        mass = 0
        for i in range(self.natom):
            mass += self.atomic_mass(i)
        return mass
        

    def atomic_name(self, index):
        """
        INPUT: index\n
        OUPTUT: atomic_name of the index\n
        """
        return periodic_table.num2name(self.znucl[self.typat[index]-1])

    def atomic_number(self, index):
        """
        INPUT: index\n
        OUPUT: atomic number of the index
        """
        return self.znucl[self.typat[index]-1]

    def atomic_mass(self, index):
        """
        INPUT : index\n
        OUTPUT: atomic mass of the index 
        """
        return periodic_table.num2mass(self.znucl[self.typat[index]-1])

    def atomic_group(self,index):
        """
        INPUT : index\n
        OUTPUT: atomic group of the index\n
        """
        return periodic_table.num2group(self.atomic_number(index=index))

    def selective_bool(self):
        """
        INPUT: None
        OUTPUT: selective as bool
        """
        selective = np.empty(shape=np.shape(self.selective), dtype=np.bool)
        selective[np.array(self.selective) == 'T'] = True
        selective[np.array(self.selective) == 'F'] = False
        return selective

    def near_atomic_index(self, index, radius):
        """
        INPUT : index: center atom, radius: search space
        OUTPUT: vindex
        """
        vindex = np.array(list(range(self.natom)))
        norms = norm(self.xcart[index] - self.xcart, axis=1)
        # 自分のindexを取り除く
        vindex_near_bool = np.logical_and(0<norms, norms <= radius)
        vindex_near = vindex[vindex_near_bool]
        return vindex_near

def poscar2structure_selective(filename="POSCAR"):
    """
    INPUT: POSCAR file name
    OUTPUT: structure
    FUNCTION: 
    """
    n_line = 0         # line number
    r = 1.0        # mulutipricity
    rprim = np.zeros([3, 3], dtype=np.float)    # lattice vector
    natom = 0   # the number of atoms
    ntypat = 0        # the number of atom types
    znucl = []        # atom number array
    typat = []        # atom number index array
    coordinate = ""        # coordinate type Direct or Cartesian
    flag_sd = 0             # flag_selective_dynamics
    xcart = []          # atom positions
    vselective_flag = []      # selective dynamics
    with open(filename, "r") as f:
        for line in f:
            sl = line.split()  # splitted line
            if n_line == 1:
                r = np.float(sl[0])
            # lattice vectors line
            elif n_line >= 2 and n_line <= 4:
                rprim[n_line-2] = np.array(sl, dtype=np.float)
            # atom name line
            elif n_line == 5:
                znucl = np.array([periodic_table.name2num(sl[i])
                                  for i in range(len(sl))])
                ntypat = len(znucl)
            # the number of each atom's line
            elif n_line == 6:
                # the number of atoms one znucle
                n_each_atom = np.array(sl, dtype=np.int)
                natom = np.sum(n_each_atom)
                xcart = np.zeros([natom, 3], dtype=np.float)
                # make typat
                for i in range(len(n_each_atom)):
                    typat = np.append(typat, np.array(
                        [i+1 for j in range(n_each_atom[i])]))
                typat = np.array(typat, dtype=np.int)
            # Selective dynamics or not
            elif n_line == 7 and (line[0]=='s' or line[0]=='S'): 
                flag_sd = 1
            # Coordinate line
            elif n_line == 7+flag_sd:
                coordinate = line.replace("\n", "")
            # x, y, z line
            elif n_line >= 8+flag_sd:
                # save positions
                if n_line-8-flag_sd == natom:
                    break
                idx_xcart = n_line-(8+flag_sd)
                xcart[idx_xcart] = np.array(sl[0:3], dtype=np.float)
                if flag_sd==1:
                    vselective_flag.append(sl[3:6])

            n_line += 1

    # Magnificate rprim
    rprim = r * rprim
    # change coordinate Direct to Cartesian
    # [dx, dy, dz][va] = [cx, cy, cz]
    #             [vb]
    #             [vc]
    if coordinate[0] == "D" or coordinate[0] == "d":
        xred = deepcopy(xcart)
        xcart = direct2cartesian(rprim=rprim, xred=xred)
    elif coordinate[0] == "C" or coordinate[0] == "c":
        xred = cartesian2direct(rprim=rprim, xcart=xcart)

    struct0 = structure(natom=natom, ntypat=ntypat, typat=typat, znucl=znucl, xcart=xcart, xred=xred, rprim=rprim)
    return struct0, vselective_flag


def poscar2structure(filename="POSCAR"):
    """
    INPUT: POSCAR file name
    OUTPUT: structure
    FUNCTION: 
    """
    struct0, vselective_flag = poscar2structure_selective(filename=filename)
    return struct0

def xyz2structure(filename="test.xyz", cell=None, margin=None, is_xyz2=False):
    """
    INPUT: fn1: XYZ file, cell: unit cell vectors 3x3 (a1, a2, a3) in angst, margin: is active when not use "cell", is_xyz: When line 1 and 2 is removed, True\n
    OUTPUT: structure\n
    """
    # atom_position = ("ATOM", X, Y, Z)
    atom_position = []
    with open(filename, "r") as f:
        cnt = 0
        for line in f:
            if (is_xyz2==False and cnt >= 2) or is_xyz2:
                atom_position.append(line.split()[0:4])
            cnt+=1
    natom = len(atom_position)
    atom_position = np.reshape(atom_position, (natom, 4))
    # Vector of atomic number 
    vatom_number = np.array([periodic_table.name2num(atom) for atom in atom_position[:, 0]])
    xcart = np.array(atom_position[:, 1:4], dtype=np.float)
    # idx sort by number of atom
    vidx = my_argsort(vatom_number)
    vatom_number = vatom_number[vidx]
    xcart = xcart[vidx]
    print(vidx)
    znucl = np.unique(vatom_number)
    ntypat = len(znucl)

    typat = []
    [typat.append(znucl.tolist().index(atom_number)+1) for atom_number in vatom_number]
    typat = np.array(typat)

    # Set cell and margin of unitcell
    if cell is None:
        x_min=np.min(xcart[:, 0])
        y_min=np.min(xcart[:, 1])
        z_min=np.min(xcart[:, 2])

        x_max=np.max(xcart[:, 0])
        y_max=np.max(xcart[:, 1])
        z_max=np.max(xcart[:, 2])

        x_len = x_max - x_min
        y_len = y_max - y_min
        z_len = z_max - z_min
        if margin is None:
            margin = np.array([5, 5, 5])
        elif type(margin) is list:
            margin = np.array(margin)
        xcart += margin/2 - np.array([x_min, y_min, z_min])
        rprim = np.array([[x_len+margin[0], 0, 0], [0, y_len+margin[1] , 0], [0, 0, z_len+margin[2]]])
    else:
        rprim = cell

    xred = cartesian2direct(rprim=rprim, xcart=xcart)
    return structure(natom=natom, ntypat=ntypat, typat=typat, znucl=znucl, xcart=xcart, xred=xred, rprim=rprim)


def mol2structure(filename="test.mol", cell=None, margin=None, mode="POSITION"):
    """
    INPUT: filename: mol file, cell: unit cell vectors 3x3 (a1, a2, a3) in angst, margin: is active when not use "cell", mode: "POSITION" or "BOND"\n
    OUTPUT: structure\n
    """
    if mode == "POSITION":
        atoms = []
        positions = []
        with open(filename, "r") as fh:
            cnt = 0
            for line in fh:
                spl = line.split()
                if cnt > 3 and len(spl) == 16:
                    atoms.append(spl[3])
                    positions.append(spl[0:3])
                cnt += 1
        natom = len(atoms)
        typat2 = [periodic_table.name2num(atom) for atom in atoms]
        # sort
        idxs = my_argsort(varray=typat2)
        typat2 = np.array(typat2)[idxs]
        positions = np.array(positions, dtype=np.float)
        positions = positions[idxs]
        znucl = np.unique(typat2)
        ntypat = len(znucl)
        typat = np.copy(typat2)
        for i in range(ntypat):
            typat[typat2 == znucl[i]] = i+1
        xcart = positions

        # Set cell and margin of unitcell
        if cell is None:
            x_min = np.min(xcart[:, 0])
            y_min = np.min(xcart[:, 1])
            z_min = np.min(xcart[:, 2])

            x_max = np.max(xcart[:, 0])
            y_max = np.max(xcart[:, 1])
            z_max = np.max(xcart[:, 2])

            x_len = x_max - x_min
            y_len = y_max - y_min
            z_len = z_max - z_min
            if margin is None:
                margin = np.array([5, 5, 5])
            elif type(margin) is list:
                margin = np.array(margin)
            xcart += margin/2 - np.array([x_min, y_min, z_min])
            rprim = np.array(
                [[x_len+margin[0], 0, 0], [0, y_len+margin[1], 0], [0, 0, z_len+margin[2]]])
        else:
            rprim = cell


    xred = cartesian2direct(rprim=rprim, xcart=xcart)
    return structure(natom=natom, ntypat=ntypat, typat=typat, znucl=znucl, xcart=xcart, xred=xred, rprim=rprim)




def generate_poscar(structure_, filename, isSelective=False, coordinate="Cartesian", vselective_flag=[]):
    """
    INPUT: atomic infomation, filename natom, znucl, ntypat, typat, xcart, rprim, fn\n
    OUTPUT: none\n
    FUNCTION: generate POSCAR\n
    """
    # Argument sort
    vidx = my_argsort(structure_.typat)
    poscar = str(structure_.natom) + " atoms system\n"
    poscar += str(1.0)+"\n"
    for i in range(3):
        poscar += "{0: .16E}    {1: .16E}    {2: .16E}\n".format(*structure_.rprim[i])
    poscar += " ".join([periodic_table.num2name(znucl) for znucl in  structure_.znucl]) + "\n"
    for i in range(1, structure_.ntypat+1):
        poscar += str(np.sum(structure_.typat == i)) + " "
    poscar += "\n"

    if isSelective:
        poscar += "Selective dynamics\n"

    poscar += coordinate + "\n"
    temp_x = []  # matrix of atomic position
    if coordinate == "Direct" or coordinate == "direct":
        temp_x = structure_.xred
    else:
        temp_x = structure_.xcart
    for i in range(structure_.natom):
        poscar += "{0: .16E}    {1: .16E}    {2: .16E}".format(*temp_x[vidx][i])
        if isSelective:
            poscar += " {0} {1} {2}".format(*vselective_flag[i])
        poscar += "\n"
    generate_file(filename, poscar)


def generate_xyz(structure_, filename):
    """
    INPUT: structure, filename: .xyz
    OUTPUT: None
    FUNCTION: Generate xyz format file
    """
    with open(filename, "w") as f:
        f.write(str(structure_.natom)+"\n")
        f.write("comment line\n")
        for i in range(structure_.natom):
            f.write("{0:2} {1: 1.16E} {2: 1.16E} {3: 1.16E}\n".format(structure_.atomic_name(index=i), *structure_.xcart[i]))

def generate_mol(structure_, filename):
    """
    INPUT: structure, filename: .mol
    OUTPUT: None
    FUNCTION: Generate mol format file
    """

    natom = structure_.natom
    xcart = structure_.xcart
    # Bond
    bond_max = 1.6
    bond_table = []
    for i in range(natom):
        for j in range(i):
            if norm(xcart[i]-xcart[j]) < bond_max:
                bond_table.append([j+1, i+1])
    nbond = len(bond_table)

    countline = "{0:3d}{1:3d}{2:3d}{3:3d}{4:3d}{5:3d}{6:3d}{7:3d}{8:3d}{9:3d}{10:3d}".format(*[natom, nbond, 0, 0, 0, 0, 0, 0, 0, 0, 1]) + " V2000\n"
    with open(filename, "w") as f:
        f.write("id\n")
        f.write("generate_mol()\n")
        f.write("comment\n")
        f.write(countline)
        for i in range(natom):
            f.write("{0:10.4f}{1:10.4f}{2:10.4f}{3:^3}".format(*xcart[i], structure_.atomic_name(index=i)))
            for j in range(12):
                f.write("{0:3d}".format(0))
            f.write("\n")
        for i in range(nbond):
            f.write("{0:3d}{1:3d}{2:3d}".format(*bond_table[i], 1))
                        
            for j in range(4):
                f.write("{0:3d}".format(0))
            f.write("\n")


        f.write("M  END\n")

def generate_lammpsdata(structure_, filename):
    """
    INPUT : structure, filename
    OUTPUT: lammps 's atomic data format
    """
    natom = structure_.natom
    znucl = structure_.znucl
    ntypat = structure_.ntypat
    typat = structure_.typat
    rprim = structure_.rprim
    xcart = structure_.xcart
    # 格子ベクトルのノルム
    snorm = norm(rprim, axis=1)
    # 格子ベクトルのなす角度
    xy = 90 - np.mod(180/np.pi*np.arccos(np.dot(rprim[0], rprim[1])/(snorm[0]*snorm[1])), 180)
    yz = 90 - np.mod(180/np.pi*np.arccos(np.dot(rprim[1], rprim[2])/(snorm[1]*snorm[2])), 180)
    zx = 90 - np.mod(180/np.pi*np.arccos(np.dot(rprim[2], rprim[0])/(snorm[2]*snorm[0])), 180)
    # 質量
    masses = []
    for i in range(ntypat):
        masses.append(periodic_table.num2mass(znucl[i]))
    # 原子数とタイプ数とユニットセル
    data  = "{:d} atoms\n".format(natom)
    data += "{:d} atom types\n".format(structure_.ntypat)
    data += "{: .9f} {: .9f} xlo xhi\n".format(0, snorm[0])
    data += "{: .9f} {: .9f} ylo yhi\n".format(0, snorm[1])
    data += "{: .9f} {: .9f} zlo zhi\n".format(0, snorm[2])
    data += "{0: .9f} {1: .9f} {2: .9f} xy yz zx\n".format(xy, yz, zx)
    data += "\n"
    data += "Masses\n\n"
    for i in range(ntypat):
        data += "{:d} {:f}\n".format(i+1, masses[i])
    data += "\n"
    # 原子座標とか
    data += "Atoms\n\n"
    for i in range(natom):
        data += "{0:d} {1:d} {2:} {3:.9f} {4: .9f} {5: .9f}\n".format(i+1, typat[i], 0.0, *xcart[i])
    generate_file(fn=filename, strs=data)


def generate_file(fn, strs):
    """
    INPUT: output filename, writing string\n
    FUNCTION: generate file, and back up if already exit fn\n
    """
    with open(fn, "w") as f:
        f.write(strs)
        print("generated " + os.path.abspath(fn))

def my_argsort(varray):
    """
    INPUT: varray: array\n
    FUNCTION: return index array of sorted varray\n
    """
    vidx = list(range(len(varray)))
    vidx_varray = []
    vunique = np.unique(varray).tolist()
    for unique in vunique:
        [vidx_varray.append(idx) for idx in vidx if varray[idx] == unique]
    return vidx_varray

def direct2cartesian(rprim, xred):
    """
    INPUT: rprim: lattice vector, xcart: set of position\n
    OUTPUT: xcart\n
    """
    # change coordinate Direct to Cartesian
    # [dx, dy, dz][va] = [cx, cy, cz]
    #             [vb]
    #             [vc]
    return np.matmul(xred, rprim)

def cartesian2direct(rprim, xcart):
    """
    INPUT: rprim: lattice vector, xred: set of position\n
    OUTPUT: xred\n
    """
    # change coordinate Direct to Cartesian
    # [dx, dy, dz][va] = [cx, cy, cz]
    #             [vb]
    #             [vc]
    return np.matmul(xcart, np.linalg.inv(rprim))


def rodrigues_rotation(vp, vn, a):
    """
    rodrigues rotation
    rotate vp, a [rad] around vn vector\n
    INPUT: position vector: vp, rotate axis: vn, rotate radian: a\n
    OUTPUT: rotated position (3, )\n
    ([rod1][x1, ..., xN])^T = [x, y, z]([rod1])^T  = [x1', y1', z1']
    ([rod2][y1, ..., xN])              ([rod2])      [.,     .,    ]
    ([rod3][z1, ..., zN])              ([rod3])      [xN', yN', zN']
    """
    rod1 = [np.cos(a)+vn[0]**2 * (1-np.cos(a)), vn[0]*vn[1]*(1-np.cos(a)) -
            vn[2]*np.sin(a), vn[0]*vn[2]*(1-np.cos(a))+vn[1]*np.sin(a)]
    rod2 = [vn[1]*vn[0]*(1-np.cos(a))+vn[2]*np.sin(a), np.cos(a)+vn[1]
            ** 2 * (1-np.cos(a)), vn[1]*vn[2]*(1-np.cos(a))-vn[0]*np.sin(a)]
    rod3 = [vn[2]*vn[0]*(1-np.cos(a))-vn[1]*np.sin(a), vn[2]*vn[1] *
            (1-np.cos(a))+vn[0]*np.sin(a), np.cos(a)+vn[2]**2 * (1-np.cos(a))]
    rod = np.array([rod1, rod2, rod3])
    return np.matmul(vp, np.transpose(rod))

def gravity_center(structure_):
    """ 
    INPUT: structure_: class structure\n
    OUTPUT:gravity_center\n
    """
    vatomic_number = [structure_.znucl[typat-1] for typat in structure_.typat]
    vatomic_mass = [periodic_table.num2mass(atomic_number) for atomic_number in vatomic_number]
    # print(vatomic_number)
    return np.dot(vatomic_mass, structure_.xcart)/np.sum(vatomic_mass)
    
    
def add(structure1, structure2):
    """
    INPUT : structure1: structure, structure2: structure
    OUTPUT: structure1 + structure2
    """
    assert structure1.volume() > structure2.volume(), "structure1.volume() must be larger than structure2.volume()"

    natom = structure1.natom + structure2.natom
    # typatから原子番号の列に変換して結合
    types3 = [structure1.znucl[typat-1] for typat in structure1.typat]
    [types3.append(structure2.znucl[typat-1]) for typat in structure2.typat]
    znucl = np.unique(np.append(structure1.znucl, structure2.znucl))
    ntypat = len(znucl)
    typat = []
    [typat.append(znucl.tolist().index(type3)+1) for type3 in types3]
    typat = np.array(typat)
    vidx = my_argsort(typat)
    typat = typat[vidx]
    rprim = structure1.rprim
    # Add xcart
    xcart = np.r_[structure1.xcart, structure2.xcart]
    # Add xred
    xred = np.r_[structure1.xred, cartesian2direct(rprim=rprim, xcart=structure2.xcart)]
    xcart = xcart[vidx]
    xred = xred[vidx]

    return structure(natom=natom, ntypat=ntypat, typat=typat, znucl=znucl, xcart=xcart, xred=xred, rprim=rprim)

def shift(structure_, vshift):
    """
    INPUT: vshift: shift vector
    OUTPUT: Shifted structure
    FUNCTION: add vshift all xcart and recalclate xred
    """
    structure2 = structure_.copy()
    structure2.xcart = structure2.xcart + vshift
    structure2.xred = cartesian2direct(rprim=structure2.rprim, xcart=structure2.xcart)
    return structure2

def find_H2O(structure_):
    """
    INPUT:    structure_: structure object
    OUTPUT:   index of H2O atom
    """
    # O - H bond length is almost 1 angst
    bond_OH = 1.2
    sO_index = []   # Set of O
    sH_index = []   # Set of H
    # For calculate unit cell image
    sR = []
    for l in range(-1, 2):
        for m in range(-1, 2):
            for n in range(-1, 2):
                if l!=0 or m!=0 or n!=0:
                    sR.append(l*structure_.rprim[0] + m*structure_.rprim[1] + n*structure_.rprim[2])

    # Set of H2O indexs
    sH2O_index = []
    # Find O and H index
    for i in range(structure_.natom):
        atomic_name = structure_.atomic_name(index=i)
        if atomic_name == 'O':
            sO_index.append(i)
        elif atomic_name == 'H':
            sH_index.append(i)
    # Find H2O
    for O_index in sO_index:
        x_O = structure_.xcart[O_index]
        nH = 0
        nearH_index = []

        for H_index in sH_index:
            x_H = structure_.xcart[H_index]
            # Condider unitcell image
            min_norm = norm(x_O - x_H)
            for R in sR:
                x_O_temp = x_O + R
                if min_norm > norm(x_O_temp - x_H):
                    min_norm = norm(x_O_temp - x_H)

            if min_norm < bond_OH:
                nH += 1
                nearH_index.append(H_index)
        # H2O
        if nH == 2:
            sH2O_index.append([O_index, nearH_index[0], nearH_index[1]])
    
    return sH2O_index

def addmol2h2o(structure1, structure2, length=1.2):
    """
    INPUT : structure1: molecule, structure2: H2O
    OUTPUT: structure1 + structure2
    """
    sH2O_index = find_H2O(structure_=structure2)
    assert len(sH2O_index)*3 == structure2.natom, "structure2.natom must be len(H2O)*3"

    # Find to remove H2O indexs
    vidx = []
    for i in range(structure2.natom):
        xcart = structure2.xcart[i]
        norms = norm(structure1.xcart - xcart, axis=1)
        if True in (norms < length):
            vidx.append(i)

    vidx_remove = []
    for i in range(len(sH2O_index)):
        for index in sH2O_index[i]:
            if index in vidx:
                vidx_remove.extend(sH2O_index[i])
                break
    vidx_remove = np.sort(vidx_remove)[::-1]
    
    structure3 = structure2.copy()
    for idx in vidx_remove:
        structure3 = structure3.remove(index=idx)
    return add(structure1=structure3, structure2=structure1)



def supercell(structure_, nx, ny, nz):
    """
    INPUT:    structure_: structure object, nx: multiplicity of x, ny: ... y, nz: ...z\n
    OUTPUT:   supercell structure
    """

    natom = structure_.natom * nx * ny * nz
    ntypat = structure_.ntypat
    rprim0 = structure_.rprim
    rprim = np.array([nx * rprim0[0], ny * rprim0[1], nz * rprim0[2]])
    xcart = []
    xcart0 = structure_.xcart
    typat0 = structure_.typat
    typat = []
    znucl = structure_.znucl
    
    for z in range(nz):
        for y in range(ny):
            for x in range(nx):
                for i in range(structure_.natom):
                    subxcart = xcart0[i]
                    subtypat = typat0[i]
                    xcart.append(subxcart + np.matmul(np.array([x, y, z]), rprim0))
                    typat.append(subtypat)

    typat = np.array(typat)
    xcart = np.array(xcart)
    vidx = my_argsort(typat)
    # Sort by atomic number
    typat = typat[vidx]
    xcart = xcart[vidx]
    xred = cartesian2direct(rprim=rprim, xcart=xcart)
    

    return structure(natom=natom, ntypat=ntypat, typat=typat, znucl=znucl, xcart=xcart, xred=xred, rprim=rprim, selective=[])

