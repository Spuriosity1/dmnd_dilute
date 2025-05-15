# INSTALLING


## Step 1: latticelab, a required dependency
```bash
git clone https://github.com/Spuriosity1/latticelab.git && cd latticelab
meson setup build 
ninja -C build install
```
If you can't touch `/usr/local`, instead do
```bash
meson setup build -Dprefix="/your/install/prefix"
```

## Step 2: install `dmnd_dilute`

```bash
git clone https://github.com/Spuriosity1/dmnd_dilute.git && cd dmnd_dilute
meson setup build
# if this fails to find latticelab, rune 
export PKG_CONFIG_PATH="/your/install/prefix/lib/pkgconfig"
ninja -C build
```

This will produce an executable in `build` called `dmnd_dilute`. To get information on how to call it, run
```bash
$ ./dmnd_dilute -h
```
which will print
```
Usage: dmndlat [--help] [--version] --output_dir VAR [--verbosity VAR] [--force] [--neig
hbours VAR...] [--delete_spins VAR...] [--dilution_prob VAR] [--seed VAR] [--save_lattic
e] Z1 Z2 Z3

Positional arguments:
  Z1                   First lattice vector in primitive units (three integers)  [nargs:
 3] 
  Z2                   Second lattice vector in primitive units (three integers) [nargs:
 3] 
  Z3                   Third lattice vector in primitive units (three integers) [nargs: 
3] 

Optional arguments:
  -h, --help           shows help message and exits 
  -v, --version        prints version information and exits 
  -o, --output_dir     Path to output [required]
  -v, --verbosity      [nargs=0..1] [default: 0]
  -f, --force          Overwrites output files 
  -n, --neighbours     [nargs: 0 or more] [default: {2}]
  -d, --delete_spins   Specific spin indexes to delete. [nargs: 1 or more] [default: {}]
  -p, --dilution_prob  Probability of deleting spin i 
  --seed               64-bit int to seed the RNG 
  --save_lattice       Flag to save the full lattice file 
```

# EXAMPLE USAGE
```bash
$ mkdir -p ../tmp
$ build/dmnd_dilute 4 0 0 0 4 0 0 0 4 -p 0.1 -o ../tmp --neighbours 2 4 --save_lattice
$ python3 ../lattice_indexing_lib/scripts/visualise.py "../tmp/Z1=4,0,0;Z2=0,4,0;Z3=0,0,4;nn=2,4;p=0.1000;seed=0;.lat.json" links plaqs
```

Yes, it puts semicolons in filenames, and yes, I do regret it.


