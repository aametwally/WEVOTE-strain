# WEVOTE2
Extension of WEVOTE algorithm to identify  microbes on strain level
This is currently a work in progress - current plans for developing the software are as follows:



### Installation.

The first step to installing WEVOTEstrain is to install WEVOTE: https://github.com/aametwally/WEVOTE. 
Once WEVOTE is complete, simply clone WEVOTE-strain into a directory. Once in the main directory run the following commands:
```
cd src/
make
```
This will create the executable file wvstrain.

### Selecting a species

In order the tell WEVOTE strain which species you would like to select, please navigate to https://www.ncbi.nlm.nih.gov/taxonomy and search for the species. Navigating to the species' page will provide you with the taxonomicID. For example, Shigella dysenteriae has a taxonomicID of 622.

### Running

#### Preprocessing

The first point at which the extension can be run is at the completion of WEVOTE. 

At the end of WEVOTE there is an output folder named by the user which we will refer to as MY_WEVOTE_OUTPUT. In this folder locate the file named $MY_WEVOTE_OUTPUT_WEVOTE_Details.txt. This file will contain the readID along with the individual taxonomicIDs assigned to the read from each method. Finally in the last column will be the final taxonomicID assigned to the read by WEVOTE.

In order to run WEVOTEstrain, a new FASTA file must be generated from the old one, using the WEVOTE_Details.txt file as a guide to select individual reads that match the species of interest. This will be done by  post_WEVOTE_script.sh, located in the scripts folder. The usage of the script is:
```
./post_WEVOTE_script.sh $TAXID $MY_WEVOTE_OUTPUT_WEVOTE_Details.txt $ORIGINAL_FASTA $OUTFILE_FASTA

```

The output of this script will be a file named $OUTFILE_FASTA that contains reads only for the species-of-interest.

#### Executing WEVOTE strain

Once the fasta file is generated from the preprocessing, the script simply needs to be fed as an argument to the executable wvstrain

```
./wvstrain $OUTFILE_FASTA -o $WVSTRAIN_OUTPUT
```

The output of the file is a list of each read and its corresponding strain-level taxonomicID.





