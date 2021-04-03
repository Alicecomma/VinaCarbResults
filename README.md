# VinaCarbResults
Windows (.NET/C#) console application to combine `.pdbqt` results files from a specific folder setup into more convenient collections of formats, taking into account whether oxygens are present inside of a pre-defined sphere.

## Use
Can be called from cmd.exe: `VinaCarbResults [0] [1] [2]`
- `[0]`: Directory (group of poses) to enter (`/*`)
- `[1]`: Enzyme results to check (given `*.enzyme`, check `/dir/pos/*.pdbqt`)
- `[2]`: Pose (element of directory) to enter (`/dir/*`)

### Directory setup
Set up your directories as follows:

    VinaCarbResults.exe
    [1].enzyme
    /[1]
    /[0]
    /[0]/[2]
    /[0]/[2]/[1].log
    /[0]/[2]/[1].pdbqt

For example (for PDB: '[3KLK](https://www.rcsb.org/structure/3KLK)', directory '[ENumbers](https://en.wikipedia.org/wiki/E_number)', pose '[E100](https://en.wikipedia.org/wiki/Curcumin)') after an AutoDock Vina run:

    3klk.enzyme
    /3klk/
    /ENumbers
    /ENumbers/E100
    /ENumbers/E100/3klk.log
    /ENumbers/E100/3klk.pdbqt

### [1].enzyme

    --radius 1.600 19.829 31.778 74.862

The above defines a sphere `[float radius] [float x] [float y] [float z]` that is checked to contain oxygens (`.endsWith("O") || .endsWith("OA")`)

## Further reading

Please see the Wiki (to the right)!
