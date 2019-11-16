from __future__ import print_function
import argparse
import sys
from copy import copy




if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Count UMIs from sciRNA data given stdin of UMIs. UMIs per cell are printed to STDOUT. ')
    parser.add_argument('--cell', nargs='?', type=argparse.FileType('r'), default=sys.stdin, required=True, help='Text piped in from stdin for R1.') 
    args = parser.parse_args()

    umis_per_cell = dict() 
    
    for line_number,line in enumerate(args.cell):        
        line = line.strip().split()
        allowed_regions = ["exonic","intronic"]
        
        if (line[-1] in allowed_regions):
            elements = line[0].split("|")
            cell_id = "_".join(elements[2:5])
            if cell_id in umis_per_cell:
                umis_per_cell[cell_id] = umis_per_cell[cell_id] + 1
            else:
                umis_per_cell[cell_id] = 1

    sample = elements[1]    
    
    for Cell in umis_per_cell:
        print(sample,"\t",Cell,"\t",umis_per_cell[Cell])

