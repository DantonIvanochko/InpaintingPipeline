#!/bin/bash


source /home/danton/miniconda3/etc/profile.d/conda.sh
conda activate SE3nv
conda env list


# prompts to get variables for inpainting
introduction=$"\n\nINPAINTING WORKFLOW\n\nRFdiffusion\nProteinMPNN\nAF2-MMseqs2\n\n"
prompt_input_pdb_file=$"\n>Provide an input PDB file\n>e.g. ./file.pdb or /path/to/file.pdb\n\n[USER INPUT] PDB file:"
prompt_inpaint_length=$"\n>Provide the lengths of the shortest and longest sequence to inpaint into your model:"
prompt_inpaint_length_shortest=$"\n[USER INPUT] Shortest number of AA:"
prompt_inpaint_length_longest=$"\n[USER INPUT] Longest number of AA:"
promt_contig_map=$"\n>Provide a contig_map for inpainting\n>Use a ! symbol to signify the inpainted region\n>e.g. [P1-9/!/M1-99/0 H1-276/0] \n\n[USER INPUT] Contig_map:"


# prompting for inputs
echo -e $introduction
echo -e $prompt_input_pdb_file
read input_pdb_file
echo -e $prompt_inpaint_length
echo -e $prompt_inpaint_length_shortest
read input_linker_shortest
echo -e $prompt_inpaint_length_longest
read input_linker_longest
echo -e $promt_contig_map
read input_contig_map
echo -e $"\n\nInpainting ${input_pdb_file} at ${input_contig_map} with sequences from ${input_linker_shortest} to ${input_linker_longest} amino acids long. \n\n"




# COMMON VARIABLES
# get input pdb basename for labeling output directories and files
pdb_file_basename=$(basename $input_pdb_file .pdb)
IFS='!' read -r -a input_contig_map_split <<< "${input_contig_map}"


# Generate backbone conformations --> RFdiffusion
for linkerLength in $(seq $input_linker_shortest $input_linker_longest)
do
echo "making $linkerLength amino acid long linkers"
python3.9 ~/RFdiffusion/scripts/run_inference.py \
        inference.output_prefix=RFdiffusion_output_${pdb_file_basename}/${pdb_file_basename}_RFd_length_${linkerLength} \
        inference.input_pdb=${input_pdb_file}  \
        "contigmap.contigs=${input_contig_map_split[0]}${linkerLength}${input_contig_map_split[1]}" \
        inference.num_designs=10 \
        denoiser.noise_scale_ca=0 denoiser.noise_scale_frame=0 # Default values = 1. Reducing the noise scale improves the in silico design success rates (Watson et al 2023)
done

mv ./outputs ./RFdiffusion_output_${pdb_file_basename} # move RFd log output folder to working RFd out directory



# Generate sidechain identities --> ProteinMPNN

# move all output pdbs from rfdiffusion to own unique directory
for RFd_pdb_file in ./RFdiffusion_output_${pdb_file_basename}/${pdb_file_basename}_RFd_length_*.pdb; 
do
    RFd_pdb_file_basename=$(basename $RFd_pdb_file .pdb); 
    echo "Processing ${RFd_pdb_file_basename} file."; 
    if [ ! -d ./RFdiffusion_output_${pdb_file_basename}/${RFd_pdb_file_basename} ]
    then
    mkdir ./RFdiffusion_output_${pdb_file_basename}/${RFd_pdb_file_basename}
    fi
    mv ${RFd_pdb_file} ./RFdiffusion_output_${pdb_file_basename}/${RFd_pdb_file_basename};
done


# create ProteinMPNN output directory
if [ ! -d ./ProteinMPNN_output_${pdb_file_basename} ]
then
    mkdir ./ProteinMPNN_output_${pdb_file_basename}
fi



# create files --> parse multiple chains

if [ ! -d ./ProteinMPNN_output_${pdb_file_basename}/parse_multiple_chains_${pdb_file_basename} ]
then
    mkdir ./ProteinMPNN_output_${pdb_file_basename}/parse_multiple_chains_${pdb_file_basename}
fi



for path_to_RFdiffusion_output in ./RFdiffusion_output_${pdb_file_basename}/${pdb_file_basename}_RFd_length_*/
do
    echo $path_to_RFdiffusion_output
    RFdiffusion_output_pdb_file_basename=$(basename $path_to_RFdiffusion_output/*.pdb .pdb);
    echo "Processing ${RFdiffusion_output_pdb_file_basename} file."; 

    python3.9 ~/ProteinMPNN/helper_scripts/parse_multiple_chains.py \
            --input_path=${path_to_RFdiffusion_output} \
            --output_path="./ProteinMPNN_output_${pdb_file_basename}/parse_multiple_chains_${pdb_file_basename}/${RFdiffusion_output_pdb_file_basename}_parsed_chains.jsonl";

done


# create files --> fixed positions
# using homemade fixed positions helper script!

if [ ! -d ./ProteinMPNN_output_${pdb_file_basename}/fixed_positions_${pdb_file_basename} ]
then
    mkdir ./ProteinMPNN_output_${pdb_file_basename}/fixed_positions_${pdb_file_basename}
fi


for path_to_RFdiffusion_trb_output in ./RFdiffusion_output_${pdb_file_basename}/*.trb
do

    echo $path_to_RFdiffusion_output
    RFdiffusion_output_trb_file_basename=$(basename $path_to_RFdiffusion_trb_output .trb);
    echo "Processing ${RFdiffusion_output_trb_file_basename} file.";

    python3.9 ~/ProteinMPNN/helper_scripts/make_inpaint_fixed_positions_dict.py \
        --trb_file=${path_to_RFdiffusion_trb_output} \
        --parse_chains=./ProteinMPNN_output_${pdb_file_basename}/parse_multiple_chains_${pdb_file_basename}/${RFdiffusion_output_trb_file_basename}_parsed_chains.jsonl \
        --output_path=./ProteinMPNN_output_${pdb_file_basename}/fixed_positions_${pdb_file_basename}/${RFdiffusion_output_trb_file_basename}_fixed_positions.jsonl 

done



# run ProteinMPNN

for path_to_RFdiffusion_trb_output in ./RFdiffusion_output_${pdb_file_basename}/*.trb
do

    echo $path_to_RFdiffusion_trb_output
    RFdiffusion_output_file_basename=$(basename $path_to_RFdiffusion_trb_output .trb);
    echo "Processing ${RFdiffusion_output_file_basename} file.";

        for temp in $(seq 0.1 0.05 0.3)
        do

        echo "__________________________________"
        echo "  RFdiffusion model: ${RFdiffusion_output_file_basename}";
        echo "  Sampling temp: ${temp}";

        python3.9 ~/ProteinMPNN/protein_mpnn_run.py \
               --jsonl_path ./ProteinMPNN_output_${pdb_file_basename}/parse_multiple_chains_${pdb_file_basename}/${RFdiffusion_output_file_basename}_parsed_chains.jsonl \
               --fixed_positions_jsonl ./ProteinMPNN_output_${pdb_file_basename}/fixed_positions_${pdb_file_basename}/${RFdiffusion_output_file_basename}_fixed_positions.jsonl  \
               --out_folder ./ProteinMPNN_output_${pdb_file_basename}/${RFdiffusion_output_file_basename}_temp_${temp} \
               --num_seq_per_target 20 \
               --sampling_temp ${temp} \
               --backbone_noise "0.00" ;
    done;
done


# convert proteinmpnn output to csv format

for output_fa_file in ./ProteinMPNN_output_${pdb_file_basename}/${pdb_file_basename}_RFd_length_*_temp_*/seqs/*.fa;
do
    output_fa_file_basename=$(basename $output_fa_file .fa);
    echo "Converting ${output_fa_file} to csv.";
    sed 's/$/,/' ${output_fa_file} |
    paste -d " "  - - | # merge alternating lines
    sed -e 's/\s\+//g' | # remove whitespace
    sed 's/>//g' | # remove > symbol from fasta format
    sed 's/T=//g' | # remove > symbol from fasta format
    sed 's/sample=//g' | # remove > symbol from fasta format
    sed 's/score=//g' | # remove > symbol from fasta format
    sed 's/global_//g' | # remove > symbol from fasta format
    sed 's/seq_recovery=//g' | # remove > symbol from fasta format
    sed "s/^/${output_fa_file_basename},/" | # add input file to first position of each line
    tail -n +2 >> ./ProteinMPNN_output_${pdb_file_basename}/${output_fa_file_basename}_output.csv # sort on the 4th column (score) and write to csv
done
for output_csv_file in ./ProteinMPNN_output_${pdb_file_basename}/*_output.csv;
do
    echo "Sorting ${output_csv_file} csv file by score.";
    sort -t, -k4,4 -n ${output_csv_file} -o ${output_csv_file}; # sort on the 4th column (score)
    sed -i '1s/^/input,T,sample,score,global_score,seq_recovery,seq,\n/' ${output_csv_file} # add column names and write to csv
done




# build models with AlphaFold2 (MMseqs2, Colabfold_local)


# create AF2_MMseqs2 output directory
if [ ! -d ./AF2_MMseqs2_output_${pdb_file_basename} ]
then
    mkdir ./AF2_MMseqs2_output_${pdb_file_basename}
fi

# create directory for input csv files in AF2_MMseqs2 output directory
if [ ! -d ./AF2_MMseqs2_output_${pdb_file_basename}/${pdb_file_basename}_AF2_inputs ]
then
    mkdir ./AF2_MMseqs2_output_${pdb_file_basename}/${pdb_file_basename}_AF2_inputs
fi


# convert proteinmpnn output csv to AF2 (colab) input csv format
# we only want the top 10% scoring unique sequences


for input_csv_file in ./ProteinMPNN_output_${pdb_file_basename}/*_output.csv;
do
    echo "Generating input csv file for AF2 using ${input_csv_file}.";

    input_csv_file_basename=$(basename $input_csv_file .csv);

    tail -n +2 ${input_csv_file} | sed 's/,$//' | #remove column heads and all trailing commas
    sort -t, -k7,7 -u |   sort -t, -k4,4 | # return top scoring unique sequences and then sort on score
    head -10 | # take top 10 sequeqnces (10/100 = 10%)
    cut --complement -f 5,6 -d, | # remove rows 5 and 6 (global score and seq recovery)
    sed 's/,/_/1' | sed 's/,/_/1' | sed 's/,/_/1' | sed 's/\//:/g' | # convert first 3 commans to underscores and convert slash to colon
    sed '1s/^/id,sequence\n/' > ./AF2_MMseqs2_output_${pdb_file_basename}/${pdb_file_basename}_AF2_inputs/${input_csv_file_basename}_AF2_input.csv; # convert slash to colon and write input csv for colabfold_batch AF2

done


# run AF2-MMseqs2 (colabfold_batch)

for input_AF2_csv_file in ./AF2_MMseqs2_output_${pdb_file_basename}/${pdb_file_basename}_AF2_inputs/*_AF2_input.csv; 
do

    input_AF2_csv_basename=$(basename $input_AF2_csv_file _input.csv);
    
    models_dir=./AF2_MMseqs2_output_${pdb_file_basename}/${input_AF2_csv_basename}/;
    
    colabfold_batch --templates --amber --use-gpu-relax --num-relax 5  \
                    ${input_AF2_csv_file}  \
                    ${models_dir} # output directory for data;
                    
done





