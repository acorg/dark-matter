## Files in this directory

`New_KS.267.obs` and `New_KS.267.seq` are the original database files that
come with the GOR4 code.

`new-gor4.obs` and `new-gor4.seq` are new database files for GOR4 which were
made according to the criteria outlined in [GOR Method for Predicting Protein Secondary Structure from Amino Acid Sequence](http://www.ulb.ac.be/di/map/tlenaert/Home_Tom_Lenaerts/INFO-F-208_files/1996%20Garnier.pdf).
For a detailed description, see [Updating the GOR4 databases, 28/6/2016](https://notebooks.antigenic-cartography.org/barbara/pages/features/updating-gor4-databases.html).
In short, all sequences in the ss.txt file, available from [here](http://www.rcsb.org/pdb/static.do?p=download/http/index.html) were filtered by length (<50 AA) and resolution (< 2.5 Angstrom).
The remaining sequences were clustered using CD-HIT, using default parameters and with a sequence identity cutoff of 30%.
Sequences containing 'X', 'U', 'B', 'J', 'O' and 'Z' were removed.
The secondary structure assignements in the ss.txt files differentiate between 8 secondary structure types, whereas GOR4 only predicts 3 (helix, strand and coil). Therefore, the 8 secondary structure types had to be converted into 3.
Sequences longer than 1200 AA were removed.
The remaining sequences were written out in the format that GOR4 requires.
