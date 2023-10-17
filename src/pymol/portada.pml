# Fetch protein and remove non-necessary atoms
fetch 2yet
remove solvent
remove chain B
remove resn NAG+Cl+Na+MAN+XYP+BCN+ACE
origin chain A

# Color per secondary structure
color 0xf8f8f8, ss l+''
color 0xc69b48, ss s
color 0x193441, ss h
set cartoon_discrete_colors, 1

# Color active site residues
show sticks, (resi 1+86+175 and not name C)

# Color side chain
color magenta, (not name CA+O+N)

# Color all atoms of side chains normally but carbons
color atomic, (not elem C)

# Change sphere radius of Cu ion
set sphere_scale, 0.4

# Bonds
select uno, (resi 1 and (name ND1 or name N))+(resi 86 and name NE2)+(resi 175 and name OH)
select dos, name cu
distance bond, uno, dos
hide labels, bond
color black, bond

# Animation
movie.add_roll(8.0,axis='y',start=1)