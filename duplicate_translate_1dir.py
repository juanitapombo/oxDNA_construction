import numpy as np
import copy

path='/working_directory/'

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
    if isinstance(ids[0], list):
        ids_flat = [int(item) for sublist in ids for item in sublist]
        coords = [frame_data[i].split()[:3] for i in ids_flat]
    else:
        coords = [frame_data[i].split()[:3] for i in ids]
    return np.array(coords, dtype=float)

def read_selection_ids(sel_file):
    with open(sel_file) as f:
        ids = f.readlines()
    split_ids = [line.strip().split(',') for line in ids] ## edit the delimeter as needed!
    return split_ids

def read_topFile(topFile):
    """Reads the topology file and returns the number of bases, strands, and the topology data."""
    with open(topFile) as f:
        content = f.readlines()
    split_content = [line.split() for line in content]
    nbases = split_content[0][0]
    nstrands = split_content[0][1]
    topology = split_content[1:]
    return nbases, nstrands, topology
    
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

def fix_topology(top_list, counter=0):
    counter = counter         
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

def shift_topology(top_list, nstrands, nbases):
    shift = nstrands
    counter = nbases     
    for i in range(1, len(top_list)):
        if top_list[i-1][2] != '-1':
            top_list[i-1][2] = str(-1 + counter)
        
        top_list[i-1][3] = str(counter + 1)
        
        if top_list[i][0] != top_list[i-1][0]:
            top_list[i-1][3] = '-1'
            top_list[i][2] = '-1' 
     
        top_list[i-1][0] = str(int(top_list[i-1][0]) + shift)
        counter += 1

    top_list[-1][2] = str(-1 + counter)
    top_list[-1][3] = '-1'
    
    top_list[i][0] = str(int(top_list[i][0]) + shift)
    
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

def save_selection_ids(id_list, file_name):
    flat_list = [str(item) for inner_list in id_list for item in inner_list]
    comma_separated_values = ','.join(flat_list)
    
    with open(file_name, 'w') as f:
        f.write(comma_separated_values)


############# MAIN CODE #############

trjFile = f'{path}init.dat'
topFile = f'{path}init.top'
sel_file = f'{path}origami_indices.txt'

z_transl_factor = (4*25+50)/0.8518 # oxDNA units

#### no need to modify from here on ### 

content = read_trjFile(trjFile)
nbases, nstrands, top = read_topFile(topFile)
core_ids = read_selection_ids(sel_file)

idx = 0 
MDstep, box, energy, frame_data, idx = extract_frame_data(content,idx)

# Compute translated z values 
translated_frame_data = [
    ' '.join(columns[:2] + [f"{float(columns[2]) + z_transl_factor:.6f}"] + columns[3:])
    for columns in [line.split() for line in frame_data]
]

# Concatenate original and translated data
duplicated_frame_data = frame_data + translated_frame_data

translated_top = copy.deepcopy(top)
translated_top = shift_topology(translated_top, int(nstrands), int(nbases)) # why does this line change top?? 
duplicated_top = top + translated_top

# Write new duplicated dat and top files 
box_multiplier = 1.2
box = f"b = {' '.join([str(float(x) * box_multiplier) for x in box.split()[2:]])}\n"
writeTopConf(duplicated_top, duplicated_frame_data, 'two_disks_init', path, MDstep, box, energy)

core_ids_second_origami = [[str(int(item) + int(nbases)) for item in inner_list] for inner_list in core_ids]

save_selection_ids(core_ids, f'{path}1st_origami_indices.txt')
save_selection_ids(core_ids_second_origami, f'{path}2nd_origami_indices.txt')
