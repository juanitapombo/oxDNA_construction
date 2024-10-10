# oxDNA_construction
Scripts to construct starting configurations for oxDNA simulations


# brush_planar_origami.py: 
This script will take in an oxDNA configuration and topology file for a one layered DNA origami and output the configuration and topology file for the same origami decorated with DNA brushes in a 3D Hilbert curve (to speed up relaxation). The length and sequence of the brushes can be specified. This script will append a brush onto any staple 3'end and is not written to modify the grafting density of the same design.
