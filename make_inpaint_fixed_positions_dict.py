
# helper script for ProteinMPNN
# an adapted version of the helper script make_fixed_positions_dict.py
# limited redesigning of the inpainted sequence created by RFdiffusion
# Danton Ivanochko  - updated 2023/10/20

import numpy as np
import argparse
import json


def main(args):

    # get data from .trb file
    trb_data = np.load(args.trb_file, allow_pickle=True)
    unfixed_chain = list(set([x[0] for x in trb_data['con_hal_pdb_idx']]))
    print(unfixed_chain)
    position_list = [x[1] for x in trb_data['con_hal_pdb_idx']]
    position_list_filled = list(range(position_list[0], position_list[-1]+1))
    unfixed_position_list = [x for x in position_list_filled if x not in position_list]
    print(unfixed_position_list)
    
    # get data from parsed (multiple) chains .jsonl file
    with open(args.parse_chains, 'r') as json_file:
        json_list = list(json_file)
    my_dict = {}
    for json_str in json_list:
        result = json.loads(json_str)
        all_chain_list = [item[-1:] for item in list(result) if item[:9]=='seq_chain']
        fixed_position_dict = {}   
        for chain in all_chain_list:
            seq_length = len(result[f'seq_chain_{chain}'])
            all_residue_list = (np.arange(seq_length)+1).tolist()
            if chain not in unfixed_chain:
                fixed_position_dict[chain] = all_residue_list
            else:
                idx = np.argwhere(np.array(unfixed_chain) == chain)[0][0]
                fixed_position_dict[chain] = list(set(all_residue_list)-set(unfixed_position_list))
        my_dict[result['name']] = fixed_position_dict
    with open(args.output_path, 'w') as f:
        f.write(json.dumps(my_dict) + '\n')

if __name__ == "__main__":

    # Parse Arguments
    parser = argparse.ArgumentParser()
    parser.add_argument("--parse_chains", type=str, required=True, help="Path to parsed multiple chains .jsonl file")
    parser.add_argument("--trb_file", type=str, required=True, help="Path to RFdiffusion .trb file")
    parser.add_argument("--output_path", type=str, required=True, help="Path to the output dictionary .jsonl file")
    args = parser.parse_args()

    main(args)

