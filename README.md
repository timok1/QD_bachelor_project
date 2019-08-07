# Quantum Dot builder

This script either builds quantum dots from scratch, or attaches ligands to
existing nanocrystals. Crystals are either a zincblende or diamond cube structure.
For a global description of the building process see the dot_building_process
pdf file.

## Usage

Copy all files to a single folder, change the directory to /scripts/ and then
run

```bash
python QD_builder.py
```

and enter the requested parameters. When the script asks for filenames of
ligands, crystals or where to save the dots, don't add the file extension. All
distances are in Angstrom.

The created quantum dots will be saved in the Created_QD folder. Here, they can
be saved to their own folder or in the main folder.

### Input parameters

* **New crystal or existing**: Either build a new crystal core, or use an
already prepared one, allowing for more customisation.
* **Lattice Constant**: The lattice constant of the unit cells in Ångström.
* **Diameter**: The diameter of the crystal that will be created, in unit cells.
This means the distance from the centre of the crystal
to the cutoff planes ((111), (110), (100)) will be half this number times the
lattice constant. This number doesn't have to be an integer. If it's not, the
cutoffs will go through the unit cells.
* **Series of ligands, or one by one?**: It is possible to either build every
QD one by one, in which case multiple ligands can be added to a single QD.
When choosing a series, the script will make QD's with all possible combinations
of chosen ligands and extensions. For instance if you choose *C10H21*, *C9H19* and
*C8H17* for the bases and *COOH* and *Na* for the extensions, six QD's will be created.
* **Filename of ligand**: Name of the file of the ligand you want to add to the
crystal. Only type the name, not the file extension(.xyz).
* **Filename of extension**: Same thing as with the ligand, only this time it
adds this extension to the ligand to create a new ligand.
* **Coverage**: Determines how many ligands will be added to the QD. This is the
fraction of surface atoms that will have this ligand attached.
* **Cap**: Single atom that will cap all remaining sites. Typically hydrogen.
* **Buffer**: The distance between ligands will be at least the sum of
*buffer* and the bonding length between the two relevant atoms. This buffer can
be set for each QD individually.




### Adding your own ligands and crystals
To add your own nanocrystals add the crystals (.xyz files) to the
input_crystals folder. These don't need to be centred at the origin. The script
currently only supports zincblende and diamond cubic structures.

It is possible to add your own ligands and crystals. In order to add your own
ligands, create an .xyz file and place it in the Ligands folder. The ligand
should be aligned along the z-axis, and the atom that connects to the quantum
dot should be placed at the origin. You can do this in for example Avogadro
(https://avogadro.cc/). Connect the ligand to an atom (like Silicon) and
optimize geometry (ctrl+alt+o). Then use the align molecules tool to align it to
the z-axis. First select the just added atom, then click on the atom in the
ligand that should connect to the crystal. Align to the z-axis. Then remove the
Silicon atom, and just select the connecting atom in the ligand. Align again,
and the ligand is positioned correctly. Note: if you want to use the replacement
option, also make sure all ligands are in the same plane as well. Otherwise the
rotation value the script uses to replace ligands will be wrong and all ligands
will be rotated in the wrong direction. Also, when you tell the script to
'extend' a ligand, the extension will replace the last atom in the file. Make
sure whether this is the right one in a text editor.


Important:
This programme needs to know the bonding lengths between different kind of
elements to work. These values are stored in bonding_distances.csv. Not every
element has been entered in this file. When you try to use elements that are
not yet supported the programme will notify you, and you can manually enter them
into the csv file.
