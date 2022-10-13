
import os, sys
import shutil 
import argparse
import subprocess as sp 
import linecache
import pandas as pd


DNAWORKS_EXE = os.path.join(os.path.dirname(os.path.abspath(__file__)), "dnaworks")


def arguments():

    parser = argparse.ArgumentParser()
    parser.add_argument("-f", type=str, default='dna_file.fas', 
                        help='DNA fasta file.')
    parser.add_argument("-o", type=str, default='output', 
                        help='output directory.')
    parser.add_argument("-l", type=int, default=60, 
                        help='the length of the oligos. default is 60.')
    parser.add_argument("--tm_low", type=int, default=62, 
                        help="the lower bound of annealing Tm value. Default is 62.")
    parser.add_argument("--tm_high", type=int, default=-1, 
                        help="the lower bound of annealing Tm value. Default is 70.")
    parser.add_argument("--codon", type=float, default=30, 
                        help="the codon frequency threshold. Default is 30.")
    
    args = parser.parse_args()

    if len(sys.argv) <= 2:
        parser.print_help()    
        sys.exit(0) 
    
    return args


def read_fasta(fas_fpath):
    seqs = {}
    tag, seq = "", ""
    with open(fas_fpath) as lines:
        for l in [x for x in lines if len(x)]:
            if l[0] == ">":
                if len(tag) and len(seq):
                    seqs[tag] = seq
                tag = l[1:].split()[0]
                seq = ""
            else:
                seq = l.strip("\n").strip()
        seqs[tag] = seq
    return seqs 


def parse_output(out_fpath):
    all_oligos = []

    with open(out_fpath) as lines:
        lines = [x for x in lines]
        for (i,l) in enumerate(lines):
            if "oligonucleotides need to be synthesized" in l: 
                
                num_oligos = int(l.split()[0])
                #print(i, l, i+3, i+3+num_oligos, num_oligos)
                oligos = []

                for l_idx in range(i+3, i+3+num_oligos):
                    line = linecache.getline(out_fpath, l_idx)
                    _seq = line.split()[1]
                    _len = line.split()[2] 
                    oligos.append([_seq, _len])

                
                all_oligos.append(oligos)
            else:
                pass
                #print(i, l)

    return all_oligos


def seq_split(seq, chunk_size=50):

    seq_list = [] 
    for i in range(int(len(seq) / chunk_size)):
        _seq = seq[i*chunk_size:(i+1)*chunk_size]
        seq_list.append(_seq)
    
    seq_list.append(seq[(i+1)*chunk_size:])

    return seq_list


def prepare_inp(tm_low=62, tm_high=70, nt_len=60, codon_freq=30, seq=""):

    seqs_list = seq_split(seq)
    _seq = "\n".join(seqs_list)

    if tm_high <= tm_low:
        inp_content = f"""title "dna"
email "astrozheng@gmail.com"
timelimit 3600
melting low {tm_low} tolerance  
length {nt_len}
frequency threshold {codon_freq}
concentration oligo 1E-7 sodium 0.05 magnesium 0.002
solutions 1

nucleotide
{_seq}
//
"""
    else:
        inp_content = f"""title "dna"
email "astrozheng@gmail.com"
timelimit 3600
melting low {tm_low} high {tm_high} tolerance  
length {nt_len}
frequency threshold {codon_freq}
concentration oligo 1E-7 sodium 0.05 magnesium 0.002
solutions 1

nucleotide
{_seq}
//
"""
    return inp_content


if __name__ == "__main__":

    args = arguments()

    seqs_dict = read_fasta(args.f)
    curr_dpath = os.getcwd()

    results = []

    for (name, seq) in seqs_dict.items():
        print(f"[INFO]: processing tasks {name}")
        os.makedirs(os.path.join(args.o, name), exist_ok=True)
        c = prepare_inp(args.tm_low, args.tm_high, args.l, args.codon, seq) 

        # prepare input 
        inp_fpath = os.path.join(args.o, name, "input")
        with open(inp_fpath, 'w') as tf:
            tf.write(c)

        # run task 
        out_fpath = os.path.join(args.o, name, "LOGFILE.txt")
        if not os.path.exists(out_fpath):
            os.chdir(os.path.join(args.o, name))
            
            cmd = f"{DNAWORKS_EXE} input"
            print("Running cmd: ", cmd)
            job = sp.Popen(cmd, shell=True)
            job.communicate()
            os.chdir(curr_dpath)

        # parse output 
        _data = parse_output(out_fpath)[0]
        #print(name, _data)

        # make data  
        for i in range(len(_data)):
            _idx = name + "_" + "0" * (4 - len(str(i+1))) + str(i+1) 
            assert len(_data[i][0]) == int(_data[i][1])
            _seq = [_idx, _data[i][0], "DMT OFF", int(_data[i][1])]
            results.append(_seq) 

    results = pd.DataFrame(results, columns=['id', 'seq', 'DMT_status', 'length'])
    results.to_csv(os.path.join(args.o, "results.csv"), header=True, index=False)