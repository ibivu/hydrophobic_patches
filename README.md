# NetSurfP-2.0: Prediction of Protein Structural Features

## Prerequisites

In order to install NetSurfP-2.0, you will need the MMseqs and HHblits installed
either globally, so the "mmseqs" and "hhblits" global commands are available.

To download and install HHblits, visit the download page: https://github.com/soedinglab/hh-suite

To download and install MMseqs2, visit the download page: https://github.com/soedinglab/MMseqs2

## Installation

NetSurfP-2.0 is compatible with python version 3. It is highly advised to use 
a python virtual environment to avoid version conflicts:

You can install NetSurfP-2.0 using pip:

```sh
pip install netsurfp2/
```

Or alternatively, without root:

```sh
cd netsurfp2
python setup.py develop --user
```

You can check your installation with

```sh
netsurfp2 --help
```

If you get `command not found`, you might have issues with your `$PATH`.
Instead you can try `python -m netsurfp2` instead of `netsurfp2`.

## Usage

After installing, NetSurfP-2.0 can be used both from commandline and
be intergrated into your own personal python script. How to use either
is described below. 

### Using NetSurfP-2.0 from the commandline

#### Quickstart

```sh
netsurfp2 --csv example_out/test.csv --hhdb path/to/uniclust30_2017_04 hhblits models/hhsuite.pb example.fasta example_out/
```

The package will install the executable NetSurfP-2.0. This
executable allows the prediction surface accessibility, secondary structure, 
disorder, phi/psi dihedral angles and protein-protein binding propensity 
of amino acids

*NOTE*:

> If you have more than 10 protein sequences, it is highly recommended to use 
> mmseqs option for computational time gain. Be aware that MMseqs2 requires a 
> lot of memory.

The basic usage is the following:

```
NetSurfP-2.0 sequence.fasta

where sequence.fasta is a fasta format file containing the
amino acid sequence of an antigen or multiple antigens.

The additional following options are available:

usage: netsurfp2 [-h] [--npz NPZ] [--json JSON] [--csv CSV] [--hhdb HHDB]
             [--mmdb MMDB] [--bs BS]
             {mmseqs,hhblits} model inp out

NetSurfP-2.0 commandline frontend

positional arguments:
  {mmseqs,hhblits}  Search method
  model             Input model
  inp               Input Fasta file
  out               Output directory

optional arguments:
  -h, --help        show this help message and exit
  --npz NPZ         Export results as numpy format (default: None)
  --json JSON       Export results as JSON format (default: None)
  --csv CSV         Export results as CSV format (default: None)
  --hhdb HHDB       HHBlits Database (default: /data/Databases/hhsuite/uniclus
                    t30_2017_04/uniclust30_2017_04)
  --mmdb MMDB       MMseqs Database (default:
                    /data/Databases/mmseqs/uniclust90_2017_04)
  --bs BS           Batch size (default: 25)
```

### Integrating NetSurfP-2.0 into a python script

After installation, it is possible to import NetSurfP-2.0 into a python
script. To try it, open a python terminal in the environment where NetSurfP-2.0
is available. Then write the following;

```python
import netsurfp2 as nsp2

with open("netsurfp2/example.fasta") as f:
    protlist = parse_fasta(f)

outdir = "example_directory"
searcher  = preprocess.HHblits("models/hhblits.pb") #Specifies hhblits as model
profiles  = searcher(protlist, outdir) #Retrieves sequence profiles
nsp_model = nsp2.model.TfGraphModel.load_graph("hhblits") = nsp2.utils.init_rf() #Unloads the HHblits NetSurfP-2.0 model
results   = nsp_model.predict(profiles, outdir, batch_size=25) #Results in a list of dictionaries, keys: id,desc,seq,n,rsa,asa,phi,psi,disorder,interface,q3,q8
```