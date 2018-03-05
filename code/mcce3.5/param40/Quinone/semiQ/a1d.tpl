#1,4-dimethyl,9,10-anthraquinone        
#Agnes 10/09  
CONFLIST A1D        A1DBK A1D-1 A1DDM

NATOM    A1DDM      0
NATOM    A1DBK      0
NATOM    A1D-1      30

IATOM    A1D-1  C1  0
IATOM    A1D-1  C1M 1
IATOM    A1D-1  C2  2
IATOM    A1D-1  H2  3
IATOM    A1D-1  C3  4
IATOM    A1D-1  H3  5
IATOM    A1D-1  C4  6
IATOM    A1D-1  C4M 7
IATOM    A1D-1  C5  8
IATOM    A1D-1  H5  9
IATOM    A1D-1  C6  10
IATOM    A1D-1  H6  11
IATOM    A1D-1  C7  12
IATOM    A1D-1  H7  13
IATOM    A1D-1  C8  14
IATOM    A1D-1  H8  15
IATOM    A1D-1  C9  16
IATOM    A1D-1  C10 17
IATOM    A1D-1  O9  18
IATOM    A1D-1  O10 19
IATOM    A1D-1  C11 20
IATOM    A1D-1  C12 21
IATOM    A1D-1  C13 22
IATOM    A1D-1  C14 23
IATOM    A1D-1 1H1M 24
IATOM    A1D-1 2H1M 25
IATOM    A1D-1 3H1M 26
IATOM    A1D-1 1H4M 27
IATOM    A1D-1 2H4M 28
IATOM    A1D-1 3H4M 29

ATOMNAME A1D-1    0  C1
ATOMNAME A1D-1    1  C1M
ATOMNAME A1D-1    2  C2 
ATOMNAME A1D-1    3  H2 
ATOMNAME A1D-1    4  C3 
ATOMNAME A1D-1    5  H3 
ATOMNAME A1D-1    6  C4 
ATOMNAME A1D-1    7  C4M
ATOMNAME A1D-1    8  C5 
ATOMNAME A1D-1    9  H5 
ATOMNAME A1D-1   10  C6 
ATOMNAME A1D-1   11  H6 
ATOMNAME A1D-1   12  C7 
ATOMNAME A1D-1   13  H7 
ATOMNAME A1D-1   14  C8 
ATOMNAME A1D-1   15  H8 
ATOMNAME A1D-1   16  C9 
ATOMNAME A1D-1   17  C10 
ATOMNAME A1D-1   18  O9
ATOMNAME A1D-1   19  O10
ATOMNAME A1D-1   20  C11
ATOMNAME A1D-1   21  C12
ATOMNAME A1D-1   22  C13 
ATOMNAME A1D-1   23  C14
ATOMNAME A1D-1   24 1H1M
ATOMNAME A1D-1   25 2H1M
ATOMNAME A1D-1   26 3H1M
ATOMNAME A1D-1   27 1H4M
ATOMNAME A1D-1   28 2H4M
ATOMNAME A1D-1   29 3H4M


#1.Basic conformer Information: name, pka, em, rxn.
#23456789A123456789B123456789C
PROTON   A1D-1      0

PKA      A1D-1      0.0

ELECTRON A1D-1      1

EM       A1D-1      0.0

RXN      A1D-1      -14.351


#2.Structure connectivity
#NEUTRAL-----------
#ONNECT   conf atom  orbital  ires conn ires conn ires conn ires conn ires conn
#ONNECT |-----|----|---------|----|----|----|----|----|----|----|----|----|----|
CONNECT  A1D-1  C1  sp2       0     C1M 0     C2  0     C13
CONNECT  A1D-1  C1M sp3       0     C1  0    1H1M 0    2H1M 0    3H1M
CONNECT  A1D-1  C2  sp2       0     H2  0     C1  0     C3
CONNECT  A1D-1  H2  s         0     C2
CONNECT  A1D-1  C3  sp2       0     C2  0     H3  0     C4
CONNECT  A1D-1  H3  s         0     C3
CONNECT  A1D-1  C4  sp2       0     C3  0     C4M 0     C14
CONNECT  A1D-1  C4M sp3       0     C4  0    1H4M 0    2H4M 0    3H4M
CONNECT  A1D-1  C5  sp2       0     C12 0     H5  0     C6
CONNECT  A1D-1  H5  s         0     C5
CONNECT  A1D-1  C6  sp2       0     C5  0     C7  0     H6
CONNECT  A1D-1  H6  s         0     C6
CONNECT  A1D-1  C7  sp2       0     C8  0     C6  0     H7
CONNECT  A1D-1  H7  s         0     C7
CONNECT  A1D-1  C8  sp2       0     C7  0     C11 0     H8
CONNECT  A1D-1  H8  s         0     C8
CONNECT  A1D-1  C9  sp2       0     C11 0     O9  0     C13
CONNECT  A1D-1  C10 sp2       0     C14 0     O10 0     C12
CONNECT  A1D-1  O9  s         0     C9
CONNECT  A1D-1  O10 s         0     C10
CONNECT  A1D-1  C11 sp2       0     C8  0     C9  0     C12
CONNECT  A1D-1  C12 sp2       0     C11 0     C10 0     C5 
CONNECT  A1D-1  C13 sp2       0     C1  0     C9  0     C14
CONNECT  A1D-1  C14 sp2       0     C4  0     C10 0     C13
CONNECT  A1D-1 1H1M s         0     C1M
CONNECT  A1D-1 2H1M s         0     C1M
CONNECT  A1D-1 3H1M s         0     C1M
CONNECT  A1D-1 1H4M s         0     C4M
CONNECT  A1D-1 2H4M s         0     C4M
CONNECT  A1D-1 3H4M s         0     C4M


#3.Atom Parameters: Partial charges and Radii
#23456789A123456789B123456789C
RADIUS   A1D    C1  1.70
RADIUS   A1D    C1M 1.70
RADIUS   A1D    C2  1.70
RADIUS   A1D    H2  1.00
RADIUS   A1D    C3  1.70
RADIUS   A1D    H3  1.00
RADIUS   A1D    C4  1.70
RADIUS   A1D    C4M 1.70
RADIUS   A1D    C5  1.70
RADIUS   A1D    H5  1.00
RADIUS   A1D    C6  1.70
RADIUS   A1D    H6  1.00
RADIUS   A1D    C7  1.70
RADIUS   A1D    H7  1.00
RADIUS   A1D    C8  1.70
RADIUS   A1D    H8  1.00
RADIUS   A1D    C9  1.70
RADIUS   A1D    C10 1.70
RADIUS   A1D    O9  1.40
RADIUS   A1D    O10 1.40
RADIUS   A1D    C9  1.70
RADIUS   A1D    C10 1.70
RADIUS   A1D    C9  1.70
RADIUS   A1D    C10 1.70
RADIUS   A1D   1H1M 1.00
RADIUS   A1D   2H1M 1.00
RADIUS   A1D   3H1M 1.00
RADIUS   A1D   1H4M 1.00
RADIUS   A1D   2H4M 1.00
RADIUS   A1D   3H4M 1.00




#NEUTRAL------
#23456789A123456789B123456789C
# opt ub3lyp/lanl2dz nosymm geom=connectivity pop=chelpg scf(maxcycle=600)   Agnes 10/09
CHARGE   A1D-1  C1   0.24
CHARGE   A1D-1  C1M -0.38
CHARGE   A1D-1  C2  -0.27
CHARGE   A1D-1  H2   0.12
CHARGE   A1D-1  C3  -0.27
CHARGE   A1D-1  H3   0.12
CHARGE   A1D-1  C4   0.23
CHARGE   A1D-1  C4M -0.35
CHARGE   A1D-1  C5  -0.07
CHARGE   A1D-1  H5   0.07
CHARGE   A1D-1  C6  -0.17
CHARGE   A1D-1  H6   0.09
CHARGE   A1D-1  C7  -0.17
CHARGE   A1D-1  H7   0.09
CHARGE   A1D-1  C8  -0.08
CHARGE   A1D-1  H8   0.07
CHARGE   A1D-1  C9   0.30
CHARGE   A1D-1  C10  0.32
CHARGE   A1D-1  O9  -0.57
CHARGE   A1D-1  O10 -0.57
CHARGE   A1D-1  C11 -0.02
CHARGE   A1D-1  C12 -0.03
CHARGE   A1D-1  C13 -0.13
CHARGE   A1D-1  C14 -0.13
CHARGE   A1D-1 1H1M  0.07
CHARGE   A1D-1 2H1M  0.11
CHARGE   A1D-1 3H1M  0.11
CHARGE   A1D-1 1H4M  0.10
CHARGE   A1D-1 2H4M  0.10
CHARGE   A1D-1 3H4M  0.07





#ParaNam|Res  |Atom|Param/toggle
TRANS    A1D          t

#====================================
#        Res    #
#23456789012345678901234567890123
#-------|-----|----|----|----|----|
#SPIN     A1D   0     C9 - C10- C1
#SPIN     A1D   1     C8 - C5 - C10
#SPIN     A1D   2     C1 - C4 - C9


#=========================================================================
#        Res    #      Axis     Rotated_Atoms
#23456789012345678901234567890123
#-------|-----|----|---------|----|----|----|----|----|----|----|
ROTAMER  A1D   0     C9 - C10  WHOLE_CONF
ROTAMER  A1D   1     C11- C13  WHOLE_CONF