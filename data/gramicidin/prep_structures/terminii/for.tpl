CONFLIST FOR        FORBK 

NATOM    FORBK      3

IATOM    FORBK  CA  0
IATOM    FORBK 1HA  1
IATOM    FORBK  O   2

ATOMNAME FORBK    0  CA 
ATOMNAME FORBK    1 1HA 
ATOMNAME FORBK    2  O  




#1.Basic Conformer Information: name, pka, em, rxn.
#23456789A123456789B123456789C

#2.Structure Connectivity
#23456789A123456789B123456789C123456789D123456789E123456789F123456789G123456789H123456789I
#23456789A123456789B123456789C123456789D123456789E123456789F123456789G123456789H123456789I
#ONNECT   conf atom  orbital  ires conn ires conn ires conn ires conn
#ONNECT |-----|----|---------|----|----|----|----|----|----|----|----|----|----|----|----|
CONNECT  FORBK  CA  sp3       -1    N   0    1HA  0    O            
CONNECT  FORBK 1HA  s         0     CA
CONNECT  FORBK  O   sp2       0     CA 

#3.Atom Parameters: Partial Charges and Radii
# Radii from "Bondi, J.Phys.Chem., 68, 441, 1964."
RADIUS   FOR    CA  2.00
RADIUS   FOR   1HA  0.00
RADIUS   FOR    O   1.40

CHARGE   FORBK  CA    0.510
CHARGE   FORBK  O    -0.510
