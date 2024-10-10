## This script works for any planar (onelayered) origami centered at the origin and aligned perpendicularly to the x axis
## This script can be easily adjusted to handle multilayered but still planar origami 

from datetime import datetime
current_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")

import numpy as np
import random

path='/Users/juanitapombog/Desktop/origami_placement/haozhi_designs_9924/haozhi_designs_9924/brushed_corrected/50nm_brushed/50nuc_brushes/norep_planes/'
topFile = f'{path}50nmDISK.top'
trjFile = f'{path}50nmDISK.dat'
brush_length = 50
#sequence = ['A', 'T']  

############# FUNCTIONS #############

def read_trjFile(trjFile):
    with open(trjFile) as f:
        content = f.readlines()
    return content

def extract_frame_data(content, start_idx):
    frame_data = []
    idx = start_idx
    MDstep, box, energy = None, None, None
    while idx < len(content):
        if content[idx].startswith("t ="):
            MDstep = content[idx]
        elif content[idx].startswith("b ="):
            box = content[idx]
        elif content[idx].startswith("E ="):
            energy = content[idx]
        else:
            while idx < len(content) and not content[idx].startswith("t =") and not content[idx].startswith("b =") and not content[idx].startswith("E ="):
                frame_data.append(content[idx].strip())
                idx += 1
            break
        idx += 1
    return MDstep, box, energy, frame_data, idx

def parse_frame(frame_data, ids):
    coords = [frame_data[i].split() for i in ids]
    return np.array(coords)

def read_topFile(topFile):
    """Reads the topology file and returns the number of bases, strands, and the topology data."""
    with open(topFile) as f:
        content = f.readlines()
    split_content = [line.split() for line in content]
    nbases = split_content[0][0]
    nstrands = split_content[0][1]
    topology = split_content[1:]
    return nbases, nstrands, topology

def count_bases_per_strand(top):
    """Counts the number of bases for each strand."""
    strand_counts = {}
    for line in top:
        strand_id = line[0]  
        if strand_id in strand_counts:
            strand_counts[strand_id] += 1
        else:
            strand_counts[strand_id] = 1
    return strand_counts

def classify_strands(strand_counts, cutoff=70):
    """Classifies strands as scaffold or staple based on base count and cutoff."""
    ## cutoff=70: scaffolds should be > 70 and staples are almost always < 70 
    scaffold_strand = []
    staple_strand = []
    
    for strand_id, count in strand_counts.items():
        if count > cutoff:
            scaffold_strand.append(strand_id)
        else:
            staple_strand.append(strand_id)
    
    return scaffold_strand, staple_strand

def get_3primes_staples(top, scaffold_strand):
    """Finds the locations for 3-prime ends and excludes scaffold strands."""
    three_prime_top = []
    three_prime_index = []

    for index, line in enumerate(top):
        if line[2] == '-1':
            if line[0] not in scaffold_strand:  # Strand is not scaffold
                three_prime_top.append(line)
                three_prime_index.append(index)

    return three_prime_top, three_prime_index

def normalize_vector(v):
    norm = np.linalg.norm(v)
    if norm == 0:
        return v
    return v / norm

def write_topFile(nbases, nstrands, topology, path):
    first_row = f"{nbases} {nstrands}"
    topFile = [first_row] + [" ".join(row) for row in topology]
    with open(path, 'w') as f:
        for line in topFile:
            f.write(line + '\n')
    return topFile

def write_datFile(MDstep, box, energy, frame_data, path):
    datFile = [MDstep.strip(), box.strip(), energy.strip()] + frame_data#[" ".join(row) for row in frame_data]
    with open(path, 'w') as f:
        for line in datFile:
            f.write(line + '\n')
    return datFile

def fix_topology(top_list):
    counter = 0         
    for i in range(1, len(top_list)):
        if top_list[i-1][2] != '-1':
            top_list[i-1][2] = str(-1 + counter)
        
        top_list[i-1][3] = str(counter + 1)
        
        if top_list[i][0] != top_list[i-1][0]:
            top_list[i-1][3] = '-1'
            top_list[i][2] = '-1' 
     
        counter += 1

    top_list[-1][2] = str(-1 + counter)
    top_list[-1][3] = '-1'
    return top_list

def writeTopConf(topList, confList, systemName, path, MDstep, box, energy):

    nbases = len(topList)
    nstrands = len(set(row[0] for row in topList))
    if len(confList) == len(topList): 
        write_topFile(nbases, nstrands, topList, f'{path}{systemName}.top')
        write_datFile(MDstep, box, energy, confList, f'{path}{systemName}.dat')
        print(f'created {systemName} top & conf files')       
    else: 
        print(f'ERROR: {systemName} top & conf files are not the same length')

def firstorder_3DHilbert(length):
    
    curve = []

    current_point = np.array([0, 0, 0], dtype=float)
    curve.append(current_point.copy())
    
    segment_moves = [
    np.array([1, 0, 0], dtype=float),
    np.array([0, 0, 1], dtype=float),  
    np.array([0, 1, 0], dtype=float),  
    np.array([0, 0, -1], dtype=float), 
    np.array([1, 0, 0], dtype=float),
    np.array([0, 0, 1], dtype=float),
    np.array([0, -1, 0], dtype=float),
    np.array([0, 0, -1], dtype=float)
    ]

    full_segments = length // len(segment_moves)

    for i in range(full_segments):
        for move in segment_moves:
            current_point += move
            curve.append(current_point.copy()) 
        
    remainder = length % len(segment_moves)
    
    for j in range(remainder):
        current_point += segment_moves[j % len(segment_moves)]  # Cycle through moves
        curve.append(current_point.copy())

    return curve
############# MAIN CODE TO FIND 3PRIME ENDS ON STAPLE STRANDS FROM OXDNA STRUCTURE  & #############
          ############# PLACE DNA BRUSHES FOLLOWING A HILLBERT CURVE #############
          
## for now I will assume that the design is always centered at the origin 

curve = firstorder_3DHilbert(brush_length)

# check for desired brush sequence, if it doesn't exist default to polyT brushes
try:
    sequence
except NameError:
    brush_seq = ['T'] * brush_length
    print(f'no sequence providing, defaulting to polyT brushes of length = {brush_length} bps')
else:
    brush_seq = sequence

nbases, nstrands, top = read_topFile(topFile)
strand_counts = count_bases_per_strand(top)
scaffold_strand, staple_strand = classify_strands(strand_counts)
three_prime_top, three_prime_index = get_3primes_staples(top, scaffold_strand)
three_prime_index_safe = three_prime_index
trj = read_trjFile(trjFile)
idx = 0 
MDstep, box, energy, frame_data, idx = extract_frame_data(trj,idx)

pos_xdir_indices = []
for index in three_prime_index:
    line = frame_data[index]
    parts = line.split() 
    x_coordinate = float(parts[0])
    if x_coordinate > 0:
        pos_xdir_indices.append(index)
        
neg_xdir_indices = []
for index in three_prime_index:
    line = frame_data[index]
    parts = line.split() 
    x_coordinate = float(parts[0])
    if x_coordinate < 0:
        neg_xdir_indices.append(index)
        
#three_prime_index = [2592] # any 3' index for debugging 

three_prime_coords = parse_frame(frame_data, three_prime_index)
coords = three_prime_coords[:, :3].astype(float)
indices_coords_3prime = np.column_stack((three_prime_index, coords))

### THE BELOW SCRIPT DECORATES THE ORIGAMI WITH BRUSHES FOR ALL STAPLE THREE PRIME ENDS ### 
        ####### CREATE TOPOLOGY FILE WITH BRUSHES ######## 

flag = False
new_topology = []
previous_strand = None
last_row_index = None
counter = 0

brush_indices = []

for i, line in enumerate(top):

        current_strand = line[0]
    
        # append original line if belongs in strand without a 3' end or belongs in scaffold 
        if i not in three_prime_index and not flag:
            new_topology.append([line[0], line[1], str(int(line[2]) + counter), str(int(line[3]) + counter)])  
            last_row_index = len(new_topology) - 1
        
        if current_strand != previous_strand: 
            # if there is a 3' end and it is not the first staple with 3'
            if flag and last_row_index is not None: 
                # last row of staple make sure it ends with 5' 
                new_topology[last_row_index][3] = '-1'
                if i not in three_prime_index: 
                    new_topology.append([line[0], line[1], '-1', str(int(line[3]) + counter)]) 
                # next row after staple make sure starts with 3' and base indices are corrected for after addition 
            if not flag and last_row_index is not None:
                #print(i)
                new_topology[last_row_index][3] = '-1'
            flag = False 
        
        if i in three_prime_index: 
            # add brushes of any desired length at 3' staple ends 
            for l in range(1,brush_length+1):
                if l == 1 and counter == 0: 
                    new_entry = [line[0], brush_seq[l - 1], '-1', str(i + 1)]
                    new_topology.append(new_entry)
                    counter += 1 
                elif l == 1 and counter != 0: 
                    new_entry = [line[0], brush_seq[l - 1], '-1', str(i + 1 + counter)]
                    new_topology.append(new_entry)
                    counter += 1 
                elif l > 1:
                    new_entry = [line[0], brush_seq[l - 1], str(i - 1 + counter), str(i + 1 + counter)]
                    new_topology.append(new_entry)
                    counter += 1 
                    
                # Record the index where the new brush base is added to read when creating the dat file 
                brush_indices.append(len(new_topology) - 1)
                
            # include the line that contains the previous 3' end that is now bonded to the brush and belongs in the core 
            # of the origami 
            prev_3prime = [line[0], line[1], str(i - 1 + counter), str(i + 1 + counter)]
            new_topology.append(prev_3prime)
            flag = True  
            
        if i not in three_prime_index and flag:
            new_entry = [line[0], line[1], str(i - 1 + counter), str(i + 1 + counter)]
            new_topology.append(new_entry)
            last_row_index = len(new_topology) - 1
    
        previous_strand = current_strand

if flag is True:
    new_topology[-1][3] = '-1'


        ####### CREATE CONFIGURATION FILE WITH BRUSHES ######## 

dist = 1 # mostly arbitrary for now 
strand = 0 
counter = 0 
previous_index = None
x_sum = False 
y_sum = False
new_strand = 0 


new_frame_data = frame_data[:]
for index in sorted(brush_indices, reverse=True):
    # for the first base, take the original index and extract the line 
    # based on wether x or y is larger 
        while True:
        
            if previous_index is None or index == (previous_index - 1): 
                if strand == 0: 
                    start = three_prime_index[-1 + (-new_strand)] #+ 1 - brush_length
                    
                    line = new_frame_data[three_prime_index[-1 + (-new_strand)]] # correct bomk @ 909PM EST
                    parts = line.split() 
                    x_new, y_new = float(parts[0]), float(parts[1])  
                    
                    # this assumes 100% that the origami will be centered at the origin and perpendicular to the x axis
                    #I should modify the incoming top file to center at the origin before running this section of code             

                    x_new += dist        

                    new_line = [f"{x_new}", f"{y_new}"] + parts[2:]    
                    new_frame_data.insert(start, " ".join(new_line))

                    x_sum = True
            
                # for the remaining bases in the strand, add to the larger component 
                # progressively for each base 
                else:
                    
                    if x_sum == True and y_sum == False: 
                        x_new += dist               ## this is correct!! do not change this 
                    
                    
                    new_line = [f"{x_new}", f"{y_new}"] + parts[2:]    
                    new_frame_data.insert(start, " ".join(new_line))
                                        
                strand += 1
                counter += 1
                previous_index = index           
                break # exit while loop if true and go into the next index 
           
           
            else: # When the index belongs to a new brush 
                new_strand += 1
                x_sum = False 
                y_sum = False
                strand = 0
                previous_index = None
                continue # go through the main if loop for the same index once we start a new strand 

## write new topology and configuration files 
box_multiplier = 1
box = f"b = {' '.join([str(float(x) * box_multiplier) for x in box.split()[2:]])}\n"

new_topology = fix_topology(new_topology)

#writeTopConf(new_topology, new_frame_data, 'linear_brushes', path, MDstep, box, energy)
#raise SystemExit("Stopping script here")

origami_indices = []
for index in range(len(new_topology)):
    if index not in brush_indices:
        origami_indices.append(index)

origami_indices_str = ",".join(map(str, origami_indices))
with open(f"{path}{brush_length}_origami_indices.txt", "w") as file:
    file.write(origami_indices_str)        

### add brushes in a hilbert curve to relax faster ### 

strand = 0 
counter = 0 
previous_index = None
z_offset = 1 # fully arbitrary 

hilb_frame_data = new_frame_data[:]
alt_sign = 1 # to offset by z_offset in both negative and positive z directions wrt the original 3'

for index in sorted(brush_indices, reverse=True):

    while True:
        
        if previous_index is None or index == (previous_index - 1): 
            z_offset = random.random()
            if strand == 0: 
                                
                line = new_frame_data[index]
                parts = line.split() 
                x_init, y_init, z_const = float(parts[0]), float(parts[1]), float(parts[2])  
                
                x_new = x_init + curve[strand+1][0]
                y_new = y_init + np.sign(y_new) * curve[strand+1][1]
                z_new = z_const + np.sign(y_new) * curve[strand+1][2]
                
                new_line = [f"{x_new}", f"{y_new}", f"{z_new}"] + parts[3:]
                hilb_frame_data[index] = " ".join(new_line)
                
            else: 
                
                x_new = x_init + curve[strand+1][0]                
                y_new = y_init + np.sign(y_init) * curve[strand+1][1]
                z_new = z_const + np.sign(y_new) * curve[strand+1][2] 
                
                new_line = [f"{x_new}", f"{y_new}", f"{z_new}"] + parts[3:]
                hilb_frame_data[index] = " ".join(new_line)
                        
            strand += 1
            counter += 1 
    
            previous_index = index 
            alt_sign *= -1
            break
        
        else: # When the index belongs to a new brush 
            strand = 0
            previous_index = None
            continue

box_multiplier = 3
box = f"b = {' '.join([str(float(x) * box_multiplier) for x in box.split()[2:]])}\n"
writeTopConf(new_topology, hilb_frame_data, f'{brush_length}_full_hilbert_brushes', path, MDstep, box, energy)

print(f"Job completed on {current_time}")
