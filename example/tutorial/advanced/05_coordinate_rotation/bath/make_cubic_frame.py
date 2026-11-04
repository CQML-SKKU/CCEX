#!/usr/bin/env python3
"""Write basic/01_first_run's natural-abundance bath in the crystal (cubic) frame.

basic/01_first_run/bath/bath_Diamond_NV_natab_1 was generated in a cell whose z axis is
the NV axis [111] (genbath/UnitCell_Dia_orth). This script rotates every position about
the NV site so that the NV axis points along [111] instead, which is how a bath written
in the conventional cubic cell would look. CCEX's frame.rotation with qubit_axis [1,1,1]
undoes exactly this (R below is its matrix), so the rotated bath with the rotation on must
reproduce the original bath.

    python3 make_cubic_frame.py           # writes bath_Diamond_NV_natab_cubic_1 next to this script
"""
import math, os

HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.join(HERE, '..', '..', 'basic/01_first_run', 'bath', 'bath_Diamond_NV_natab_1')
DST = os.path.join(HERE, 'bath_Diamond_NV_natab_cubic_1')
R0 = [73.73895, 75.72808, 73.97033]          # the NV site (bath_Diamond_NV_natab_defect)

def unit(v):
    n = math.sqrt(sum(c * c for c in v)); return [c / n for c in v]
def cross(u, v):
    return [u[1]*v[2]-u[2]*v[1], u[2]*v[0]-u[0]*v[2], u[0]*v[1]-u[1]*v[0]]

# CCEX's rotation for qubit_axis q: rows ez = q, ex = q x z, ey = ez x ex, with z = [0,0,1] of the bath file.
b, q = unit([0, 0, 1]), unit([1, 1, 1])
ez = q; ex = unit(cross(q, b)); ey = cross(ez, ex)
R = [ex, ey, ez]
Rt = [[R[j][i] for j in range(3)] for i in range(3)]    # its inverse

lines = open(SRC).read().splitlines()
out = [lines[0]]                                          # the count header stays
for line in lines[1:]:
    p = line.split()
    if len(p) != 4: out.append(line); continue
    v = [float(p[i]) - R0[i] for i in range(3)]
    w = [sum(Rt[i][j] * v[j] for j in range(3)) + R0[i] for i in range(3)]
    out.append('%.17g\t%.17g\t%.17g\t%s' % (w[0], w[1], w[2], p[3]))
open(DST, 'w').write('\n'.join(out) + '\n')
print('wrote', DST, len(out) - 1, 'spins')
