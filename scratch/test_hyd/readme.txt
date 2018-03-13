Amber parameters and input files
================================

Marc Baaden

1. Prep files
=============

1.1. Ions
---------

Filename               PDBS                    Comment
...................... ....................... ...................... 
H3O.bin                H3O+.pdb                Hydronium (+1)
Eu3+.bin               EU3+.pdb                Europium (3+)
Yb3+.bin               YB3+.pdb                Ytterbium (3+)
La3+.bin               La3+.pdb                Lanthanum (3+)
NO3.bin                NO3-.pdb                Nitrate (-)
PCL.bin                PCL-.pdb                Perchlorate (-)
TFL.bin                TFL-.pdb                Perchlorate (-)

1.2. Organic Molecules
----------------------

Filename               PDBS                    Comment
...................... ....................... ...................... 
ACN.bin                HNO3.pdb                Nitric acid
CA1.bin                TETA.pdb                Scaffold for
                                               calixarenes
CON.bin                DPDA.pdb,DPPM.pdb       Calix-cone part
DPP.bin                DPDA.pdb,DPPM.pdb       diphenylphosphineoxide
                                               arm
TBU.bin                DPDA.pdb,DPPM.pdb       tert-butyl-unit
AMI.bin                TETA.pdb,DPDA.pdb       Amide moiety for
                                               calixarenes
PHO.bin                1tbp.pdb                Neutral PO from TBP
POP.bin                tbph+.pdb               Protonated POH from
                                               TBPH+ (+1)
OBU.bin                1tbp.pdb,tbph+.pdb      O-butyl moiety of TBP

2. PDB files
============

2.1. Ions
---------

Filename      residues      reference     charge        description
............. ............. ............. ............. ............. 
H3O+.pdb      H3O           [2001.3]      +1            Hydronium ion
EU3+.pdb      M             [2000.6]      +3            Europium
YB3+.pdb      M             [2000.6]      +3            Ytterbium
NO3-.pdb      NO3           [2001.3]      -1            Nitrate
PCL-.pdb      PCL           [2000.6]      -1            Perchlorate
TFL-.pdb      TFL           [2000.6]      -1            Triflate

2.2. Organic Molecules
----------------------

Filename      residues      reference     charge        description
............. ............. ............. ............. ............. 
HNO3.pdb      ACN           [2001.3]      0             Nitric acid
1tbp.pdb      PHO, OBU      [2001.3]      0             tri-n-butyl-phosphate
tbph+.pdb     POP, OBU      [2001.3]      +1            tri-n-butyl-phosphate
DPPM.pdb      CON, DPP, TBU [2000.4][2001.2] 0          tetraphosphineoxide
                                                        calix[4]
DPDA.pdb      CON, DPP, TBU [2000.4]      0             di-phosphineoxide-diamide
                                                        calix[4]

3. Errata
=========

In [2001.3] there is an error concerning the TBP and TBPH+ charges.
The Second CT in OS-CT-CT-.. has to bear a +0.06 charge and not a
-0.06 one.

4. References
=============

[2000.4]

     M. Baaden, G. Wipff, M. R. Yaftian, M. Burgard and D. Matt;
     "Cation coordination by calix{4}arenes bearing amide and/or
     phosphine oxide pendant groups: how many arms are needed to bind
     Li+ vs Na+ ? A combined NMR and molecular dynamics study.", J.
     chem. Soc., Perkin Trans. 2., 2000, 1315-1321.

[2000.6]

     M. Baaden, F. Berny, G. Wipff and C. Madic; "A molecular dynamics
     and quantum mechanics study of M3+ lanthanide cation solvation by
     acetonitrile: the role of cation size, counterions and
     polarization effects investigated.", J. Phys. Chem. A, 104, 2000,
     7659-7671.

[2001.2]

     M. Baaden, M. Burgard, C. Boehme and G. Wipff; "Lanthanide cation
     binding to a phosphoryl-calix{4}arene: the importance of solvent
     and counterions investigated by molecular dynamics and quantum
     mechanical simulations", Phys. Chem. Chem. Phys., 3, 2001,
     1317-1325

[2001.3]

     M. Baaden, M. Burgard, and G. Wipff; "TBP at the water - oil
     interface: the effect of TBP concentration and water acidity
     investigated by molecular dynamics simulations", J. Phys. Chem.
     B, 105, 2001, 11131-11141

