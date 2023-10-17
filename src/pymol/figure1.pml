# Fetch protein and remove non-necessary atoms (HOH and Cl)
fetch 5n05
remove solvent
select chlorine_nag, resn Cl+NAG
remove chlorine_nag

# Color per secondary structure
color white, ss l+''
color cyan, ss s
color deepblue, ss h

set cartoon_discrete_colors, 1

# Color active site residues
show sticks, (resi 1+78+164 and not name C)

# Color side chain
color magenta, (not name CA+O+N)

# Color all atoms of side chains normally but carbons
color atomic, (not elem C)

# Color substrate
color atomic, resn bgc
color limegreen, (resn bgc and elem C)

# Change sphere radius of Cu ion
set sphere_scale, 0.4

# Bonds
select uno, (resi 1 and (name ND1 or name N))+(resi 78 and name NE2)+(resi 164 and name OH)
select dos, name cu
distance bond, uno, dos
hide labels, bond
color black, bond


# Transparency for picture of active site
remove (not resi 1+78+164 and not resn bgc and not name cu)
color red, name CA
remove (resi 78+164 and name N)
