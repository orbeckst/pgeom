pgeom
=====

Copyright (c) 2000-2010 Oliver Beckstein <oliver.beckstein@bioch.ox.ac.uk>  
Published under the GNU Public Licence, version 3 or any higher one.  

The code in the directory ``nr`` is from [Numerical Recipes in C, 2nd
edition][NR C] and is *not free*. You _should_ have a look at the
online copy of the book and then type the required code from the
screen in your favourite text editor. The files in the ``nr``
directory are only provided for you to verify that your manual copy of
the code is correct (e.g. by using the UNIX ``diff`` utility).

[NR C]: http://www.nrbook.com/c/


Overview
--------

``pgeom`` is simple C code to generate model pores. Pores are built
from concentric rings of atomic spheres. A pore is defined by its
external radius (the maximum outer radius of the whole assembly), and
the lengths and radii of the mouth and pore regions.

The top mouth region is taken to be the mirror image of the bottom
mouth region, which sandwich the pore region. Each region is defined
by one length and one radius. The mouth radius is linearly decreased to
match up with the pore radius. With _R_<sub>M</sub> >= _R_<sub>P</sub>
this yields a funnel or hour glass shape.

``pgeom`` also contains experimental code that is meant to calculate
analytical pore volumes (see Appendix B of my thesis 
([pdf][thesis pdf])). It takes the generated structure and treats each 
atom as a CH<sub>4</sub> van der Waals united atom (in the Gromacs force field)
and then tries to calculate (essentially) the canonical partition
function for a single water molecule in this external potential.

[thesis pdf]: http://sbcb.bioch.ox.ac.uk/oliver/download/Thesis/OB_thesis_2sided_v3.pdf


Citing
------

If you use this code in published work, please cite

> O. Beckstein, Philip C. Biggin and Mark S. P. Sansom, _A hydrophobic
> gating mechanism for nanopores_, J. Phys. Chem. B **105** (2001),
> 12902-12905. DOI:[10.1021/jp012233y][]

Thank you!

[10.1021/jp012233y]: http://pubs.acs.org/doi/abs/10.1021/jp012233y



Installation
------------

Set BIN_DIR in the Makefile. 

    cd pgeom
    make
    make install

You can also compile a double precision version with -DDOUBLE (or
uncomment the 'DOUBLE := yes' in the Makefile) but that is only useful
for the pore volume stuff.

It will install:

- ``pgeom`` (some stuff has been moved in a small library that is also
  installed)
- numerical recipes routines (in directory ``./nr``) are used for pore
  volume stuff but required for compilation nevertheless; they should be
  compiled automatically
   

Usage
-----

### Hydrophobic pores ###

To generate my standard model pores (the ones used in my papers), the
commandline looks as follows (snipped from my master Makefile):
    RADIUS  := varies from 1.5 to R_INNER
    R_OUTER := 18
    R_INNER := 10
    L_PORE  := 8
    L_MOUTH := 4
    
    pgeom -v -R $(R_OUTER) -P $(RADIUS)  $(L_PORE) 
             -M $(R_INNER) $(L_MOUTH) \
          -x -o pore.pdb -s pore.itp

The benchmark pore of radius 5.5 Angstrom (right at the switching
transition):

    ./pgeom -v -R 18 -P 5.5 8 -M 10 4 -x -o pore.pdb -s pore.itp

Radii are taken to be 'solvent accessible radii', ie the centre of the
wall atom lies on a circle of radius R+rA where rA is the van der
Waals radius of the wall atom (~1.95 Å for CH4).

It is also possible to keep all the arguments in a file (say,
pore.txt) and run

    pgeom -f pore.txt -o pore.pdb -s pore.itp

In this case the input file would look like

    # example for a 5.5 A hydrophobic pore
    # (values are white-space separated)
    RADIUS 18
    
    # R1, R2: radii at top/bottom; L: length
    #      R1     R2   L
    MOUTH  10     5.5  4
    PORE    5.5   5.5  8
    MOUTH   5.5  10    4

The pore looks like this:

    *****    |     upper MOUTH
    ******** |     central PORE region
    ******** |
    *****    |     lower MOUTH


After that the pore has to be embedded into a slab of dummy atoms by a
series of gromacs commands and solvated with water (and ions if
desired). I do this in a Makefile; email me and I will send it to
you. And older version of all the necessary files is available at
<http://sbcb.bioch.ox.ac.uk/oliver/download/HyGate/>.


### Hydrophilic pores ###

I just added the appropriate charges to the itp file; the dipole
moments were taken to mimick the peptide backbone dipole moment. At an
atomic distance of 4 Å this requires a charge of +/- 0.38e. A
semi-automatic approach is taken by Scripts/porecharger.pl which takes
a list of atom numbers (alternating between positively and negatively
charged ones) and the charge and changes the itp file accordingly. You
still need to write down the list of atoms by using eg rasmol or vmd
to select the atoms.


Bugs
----

* The program can seg-fault if the pore is too big (in my defence:
  this was my first C-program ever and I only later learned about
  dynamic memory allocation...). Should be fixed eventually...
* The analytical pore volume code is highly **experimental** and
  should not be relied on to produce sensible results.
