# Generation, kMC simulation, and analysis of a finite sized finned mfi zeolite
---

#### Features:
- Create MFI framework lattice for finned or smooth nanoparticle composed of intersections sites (I), straight channel sites (S), and zig-zag channel sites (Z)
- Create a rate constants' matrix for hopping of adsorbing molecule from a site to it's neighboring sites
- Run kMC simulation by choosing initial position of the adsorbing molecule
- Collect trajectory and path length of the hopping molecule
- Obtain path length and residence times
---

#### Software required
- Python 2.7
---

#### Usage on linux
- No installation required ; Ready to use python scripts
- Save all *\*.py* scripts in a folder, say ~/kMC_zeo_scripts/  

`
cd ~/kMC_zeo_scripts/  
`  
`
git clone https://github.com/palmergroupUH/kMC_zeolite.git
`


- Go to the folder  

`
cd ~/kMC_zeo_scripts/
`


- Allow executable permissions on the *\*.py* files  

`
chmod +x *.py
`


- Allow your environment, say bash, to find the executables for use: (For bash)  

`
export PATH="~/kMC_zeo_scripts/:$PATH"
`


You can now use the kMC scripts in any folder of your choice. This folder will be the default working directory and output directory

---

#### Resources required
* In this version, the bottleneck of the resources is the generation of the full-scale MFI site lattice

    - A new version is in the pipeline that would reduce the full-scale lattice to a single unit cell

* The vitual and disk memory required for the generation of MFI site lattice (Create_lattice_August2019_nptype_v4.py) scales cubically with number of unit cells along a dimension

* If the size of zeolite lattice is 500 nm to 600 nm, the virtual memory requirement is ~ 100G+, while the disk space required is ~ 60G
---

## How to use? ##
### 1. Create a 3D matrix of the system composed of

   * locations where mfi unit cells are present
   * locations of external "gas phase" where a modified unit cell is present
   * locations of the interface between zeolite and external "gas phase" 
   
    Go to desired working directory  

    `
    cd <path-to-working-directory>
    `
    
    For a finned zeolite with,  
        * Fin Width = L1 (nm),  
        * Fin Height = L2 (nm),  
        * Pitch of the square arrangement of fins = L3 (nm)
        * Smooth cubic zeolite core of edge = L4 (nm),  

        `
        fins_fixed_gas_v5.py L1 L2 L3 L4 <prefix name>
        `

    For a smooth zeolite with,  
        * Smooth cubic zeolite core of edge = L1 (nm),  

        `
        smooth_fixed_gas.py L1 <prefix name>
        `

    To compare a smooth and finned lattice of with same edge length of cubic zeolite core:
        * Discretization of floating lengths into integral number of unit cells can cause loss or gain of one unit cell layer in the zeolite lattice
        * To overcome that, it is recommended to extract the exact number of unit cells in each dimension that make the smooth core of the finned zeolite and use that in the file smooth_500_case.py to obtain the exactly same zeolite lattice for the corresponding smooth zeolite

### 2. Construct a 3D full-scale MFI-type lattice of sites and neighbors
   * It is recommended to run this command as a batch job (assuming SLURM based batch scheduling) with sufficient job time
   * For benchmark: 100 nm case require ~ 5 min (Depending on machine). Larger sizes can take up to days.
   
    Go to directory where above files were generated

   `
   cd <path-to-working-directory>
   `

   `
   sbatch -D <path-to-working-directory> --job-name=<> --error=<> --output=<> -n 1 -N 1 -t <> --mem=<> --wrap="Create_lattice_August2019_nptype_v4.py <prefix name>"
   `

### 3. Create rate matrix for rate constants betwen possible sites
   * By default, the rates are generated as specified by Forester and Smith et al. for the case of benzene diffusion in silicalite-1 at 300K.
   * An arbitrary value of rate constant is chosen to represent 
       * HighDesorption: rate constant of desorption higher than internal diffusion rate constant(s)
       * LowDesorption: rate constant of desorption lower than internal diffusion rate constant(s)

    Choose your setting. The value used in work is "HighDesorption"

    `
    cd <path-to-working-directory>
    `  
    `
    rates.py HighDesorption
    `
    
### 4. Select possible initial positions of the adsorbing molecule in the system, for independent simulations
   * By default, only the external intersection-type sites (I) of the zeolite lattice in the system (that has an entrance/exit connection to an interface site) is a valid initial position of the adsorbing molecule.
   * Create a list of all such possible and unique external sites of the zeolite lattice and randomize their ordering.
   * This serves as the list of starting configuration of the adsorbing molecule, where each configuration is for an independent simulation

   `
   cd <path-to-working-directory>
   `  
   `
   ip_allsites_v1_mmap.py <prefix name>
   `

### 5. Run kMC simulation
   * By default, the current version assumed that the user is interested in successful adsorption into the zeolite lattice and termination of simulation when the adsorbing molecule exits from the zeolite lattice of the system
   * It is recommended to run atleast maximum of ( $10^5$ , total number of unique initial positions of the adsorbing molecule) simulations for reasonable statistics
   * The code runs the simulations serially, and therefore it is recommended to run as a batch job. Assuming SLURM,

   `
   cd <path-to-working-directory>
   `  
   
   `
   sbatch -D <path-to-working-directory> --job-name=<> --error=<> --output=<> -n 1 -N 1 -t <> --mem=<> --wrap="KMC_fast_serial_memmap.py <prefix name> <first simulation number> <last simulation number> <Desorption type: HighDesorption or LowDesorption> <pdb file name for limited visualization>"
   `
---

## Understanding the code
### Generation of lattice


```python

```


```python

```


```python

```
