import numpy as np
from numpy.linalg import norm
from atomtools.structure import rodrigues_rotation
from atomtools import material_colors
from atomtools import periodic_table
import os

def write_sphere(radius, stacks, slices, position, filename, mode="w", color_name=None):
    vnormal = []
    vvertex = []
    vvface = []
    # Normal
    for i in range(stacks+1):
        for j in range(slices):
            if i==0 or i==stacks:
                vnormal.append([0, 0, np.cos(i*np.pi/stacks)])
                break
            x = np.cos(2*np.pi*j/slices)*np.sin(i*np.pi/stacks)
            y = np.sin(2*np.pi*j/slices)*np.sin(i*np.pi/stacks)
            z = np.cos(i*np.pi/stacks)
            vnormal.append([x, y, z])
    # Vertex
    vvertex = radius * np.array(vnormal) + np.array(position)
    nvertex = len(vvertex)
    # Face
    for i in range(stacks):
        for j in range(slices):
            vvface.append([])
            nvvface = len(vvface)-1
            if i==0:
                f1 = -1
                f2 = -(np.mod(j, slices)+2)
                f3 = -(np.mod(j+1, slices)+2)
                vvface[nvvface].append([f1, f1])
                vvface[nvvface].append([f2, f2])
                vvface[nvvface].append([f3, f3])
            elif i==stacks-1:
                f1 = -nvertex
                f2 = -nvertex + np.mod(j, slices) + 1
                f3 = -nvertex + np.mod(j+1, slices) + 1
                vvface[nvvface].append([f1, f1])
                vvface[nvvface].append([f2, f2])
                vvface[nvvface].append([f3, f3])
            else:
                cnt = (i-1)*slices + j + 2
                if j==slices-1:
                    f1 = -cnt
                    f2 = -(cnt + slices)
                    f3 = -(cnt + 1)
                    f4 = -(cnt + 1 - slices)
                else:
                    f1 = -cnt
                    f2 = -(cnt + slices)
                    f3 = -(cnt + slices + 1)
                    f4 = -(cnt + 1)
                vvface[nvvface].append([f1, f1])
                vvface[nvvface].append([f2, f2])
                vvface[nvvface].append([f3, f3])
                vvface[nvvface].append([f4, f4])

    with open(filename, mode) as f:
        if not (color_name is None):
            f.write("usemtl " + color_name + "\n") 
        # Write vertex
        for vertex in vvertex:
            f.write("v {0: 5.5f} {1: 5.5f} {2: 5.5f}\n".format(*vertex))
        # Write normal
        for normal in vnormal:
            f.write("vn {0: 5.5f} {1: 5.5f} {2: 5.5f}\n".format(*normal))
        # Write face
        for vface in vvface:
            f.write("f ")
            for face in vface:
                f.write("{0:d}//{1:d} ".format(*face))     
            f.write("\n")


def write_cylinder(radius, nresolution, vp1, vp2, filename, mode="w", color_name=None):
    # Normal vector of Top circle
    vp1 = np.array(vp1)
    vp2 = np.array(vp2)
    n21 = vp2 - vp1
    height = norm(n21)
    n21 = n21 / height
    # Normal vector of side
    m = np.zeros(3)
    if n21[1] != 0:
        alpha = n21[0] / n21[1]
        m[0] = 1 / np.sqrt(1 + alpha**2)
        m[1] = -alpha * m[0]
        m[2] = 0
    elif n21[0] != 0:
        alpha = n21[1] / n21[0]
        m[1] = 1 / np.sqrt(1 + alpha**2)
        m[0] = -alpha * m[1]
        m[2] = 0
    elif n21[2] != 0:
        alpha = n21[0] / n21[2]
        m[0] = 1 / np.sqrt(1 + alpha**2)
        m[2] = -alpha * m[0]
        m[1] = 0

    # Vertex 1 of side
    vbase = radius * m

    vvertex_top = np.zeros([nresolution, 3])
    vvertex_bot = np.zeros([nresolution, 3])
    vnormal_side = np.zeros([nresolution, 3])
    # Calc vertex and normal
    for i in range(nresolution):
        vnormal_side[i] = rodrigues_rotation(
            vp=vbase, vn=n21, a=2*np.pi*i/nresolution)
        vvertex_bot[i] = radius * vnormal_side[i] + vp1
        vvertex_top[i] = radius * vnormal_side[i] + vp2
    vface = []
    for i in range(len(vnormal_side)):
        # []
        vface.append([i+1, np.mod(i, nresolution)+1+nresolution, np.mod(i+1, nresolution)+1+nresolution, np.mod(i+1, nresolution)+1, np.mod(i, nresolution
        )+1, np.mod(i+1, nresolution)+1])
    vface = np.array(-1 * np.array(vface)).tolist()
    with open(filename, mode) as f:
        if not (color_name is None):
            f.write("usemtl " + color_name + "\n") 
        for vertex_top in vvertex_top:
            f.write("v {0: 5.5f} {1: 5.5f} {2: 5.5f}\n".format(*vertex_top))
        for vertex_bot in vvertex_bot:
            f.write("v {0: 5.5f} {1: 5.5f} {2: 5.5f}\n".format(*vertex_bot))
        for normal_side in vnormal_side:
            f.write("vn {0: 5.5f} {1: 5.5f} {2: 5.5f}\n".format(*normal_side))
        for face in vface:  # Side vertex
                f.write("f {0:d}//{4:d} {1:d}//{4:d} {2:d}//{5:d} {3:d}//{5:d}\n".format(*face))
        f.write("vn {0: 5.5f} {1: 5.5f} {2: 5.5f}\n".format(*-n21))
        f.write("f")
        for i in range(nresolution):
            f.write(" {:d}//-1".format(-(i+1)))
        f.write("\n")
        f.write("vn {0: 5.5f} {1: 5.5f} {2: 5.5f}\n".format(*n21))
        f.write("f")
        for i in range(nresolution):
            f.write(" {:d}//-1".format(-(i+1+nresolution)))
        f.write("\n")

def write_obj(fn_obj, fn_mtl):
    with open(fn_obj, "w") as f:
        f.write("mtllib " + fn_mtl + "\n")

def write_mtl(filename, color_name, color, mode):
    with open(filename, mode) as f:
        f.write("newmtl " + color_name + "\n")
        f.write("Ka {0: 1.5f} {1: 1.5f} {2: 1.5f}\n".format(*color.ambient))
        f.write("Kd {0: 1.5f} {1: 1.5f} {2: 1.5f}\n".format(*color.diffuse))
        f.write("Ks {0: 1.5f} {1: 1.5f} {2: 1.5f}\n".format(*color.specular))
        f.write("Ns {: 4.5f}\n".format(1000*color.shininess/128)) # Convert value Ns: 0~1000, GL_SHININESS: 0~128

def structure2obj(filename, structure, nresolution=50):
    sbond = []

    for i in range(structure.natom):
        for j in range(i):
            length = norm(structure.xcart[i] - structure.xcart[j])
            if length < 1.6:
                sbond.append([i, j])

    dirname = os.path.dirname(filename)
    if dirname == '':
        dirname = '.'
    fn_obj = os.path.basename(filename)
    fn_mtl = fn_obj.split('.')[0] + ".mtl"
    fn_full_mtl = dirname + "/" + fn_mtl

    # Make mtl
    write_mtl(filename=fn_full_mtl, color_name="silver", color=material_colors.silver, mode="w")
    for atomic_number in structure.znucl:
        atomic_name = periodic_table.num2name(atomic_number)
        write_mtl(filename=fn_full_mtl, color_name=atomic_name, color=material_colors.name2color(atomic_name), mode="a")
    # Make obj 
    write_obj(fn_obj=filename, fn_mtl=fn_mtl)
    # Write atom by sphere
    for i in range(structure.natom):
        name = periodic_table.num2name(structure.znucl[structure.typat[i]-1])
        write_sphere(radius=0.3, stacks=nresolution, slices=nresolution,position=structure.xcart[i], filename=filename, mode="a", color_name=name)
    # Write Bond by cylinder 
    for bond in sbond:
        write_cylinder(radius=0.3, nresolution=nresolution,vp1=structure.xcart[bond[0]], vp2=structure.xcart[bond[1]], filename=filename, mode="a", color_name="silver")
