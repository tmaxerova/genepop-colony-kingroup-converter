# genepop-colony-kingroup-converter
Very amateurish but hopefully functional R scripts to convert our data between GENEPOP, COLONY and KINGROUP txt formats, as well as to generate diagnostic tables with proportion of heterozygosity for every SNP and individual in dataset. Questions about the code should be directed to Tereza Maxerová (ter.maxerova@seznam.cz), questions about the project to Michael Mikát (michael.mikat@gmail.com).
## General info
GENEPOP, COLONY and KINGROUP are types of standardised txt formats used by programs analysing relatedness of individuals in a group (e.g. hymenopteran nest). You can find sample files in formats GENEPOP, COLONY (or COLONY99) and KINGROUP (or KINGROUPdot) in a folder called "Sample_Files". You can see that for COLONY and KINGROUP, there are two versions of format. The basic version of COLONY (or just "COL") does not differentiate between male and female individuals, while in version COLONY99 all second alleles in males are substituted by "-99", as in hypemopterans the males are supposed to be haploid. Similarly in KINGROUP, the basic version does not differentiate between males and females, but in KINGROUPdot, all second alleles in males are substituted by "." .

As different programs use different formats as input or output, it is sometimes necessary to convert one type to another. Among the five scripts in this repository, you can find converters:
* from GENEPOP to COLONY/COLONY99
* from GENEPOP to KINGROUP/KINGROUPdot (where COLONY/COLONY99 are generated as by product as well)
* from COLONY to KINGROUP/KINGROUPdot
* from COLONY99 to KINGROUP/KINGROUPdot, respectively.

For the last one, there are two versions - the "SexAuto" version converts automatically all "-99"s in COLONY99 to "."s in KINGROUPdot, so you do not need to specify the sex the individuals. In the "SexAsVector" version, you can specify the male samples (where second alleles should be converted to "." in KINGROUPdot) manually, using a sex vector. 

Furthermore, it is often handy to get a simple overview of heterozygosity of SNPs and individuals in dataset. The diagnostics is produced automatically by each of the five scripts in a form of two small csv files (can be viewed e.g. in Excel), where you can see a proportion of heterozygous loci for every SNP and individual, respectively.

Each script requires one input file, then several parameters which need to be specified at the top of the R script, and produces several output files in a folder called "ResultsFolder_DATE_TIME". Below is the overview of the input files, required parameters, and the output files for each of the five scripts. At the end of the document, there is also an explained list of all parameters.

### GEN-to-COL-COL99-diag.R
Input file:
* txt GENEPOP file

Output files:
* txt COLONY file
* txt COLONY99 file
* csv COLONY table
* csv COLONY99 table
* csv table with heterozygosity diagnostics of samples
* csv table with heterozygosity diagnostics of SNPs

Required parameters:
* path_to_input_genepop_file
* output_folder_location
* sex_vector
* n_initial_empty_rows
* n_snps
* n_rubbish_rows

### GEN-to-COL-COL99-KIN-KINdot-diag.R
Input file:
* txt GENEPOP file

Output files:
* txt COLONY file
* txt COLONY99 file
* csv COLONY table
* csv COLONY99 table
* txt KINGROUP file
* txt KINGROUPdot file
* csv table with heterozygosity diagnostics of samples
* csv table with heterozygosity diagnostics of SNPs

Required parameters:
* path_to_input_genepop_file
* output_folder_location
* sex_vector
* n_initial_empty_rows
* n_snps
* n_rubbish_rows

### COL-to-KIN-KINdot-diag.R
Input file:
* txt COLONY file

Output files:
* txt KINGROUP file
* txt KINGROUPdot file
* csv table with heterozygosity diagnostics of samples
* csv table with heterozygosity diagnostics of SNPs

Required parameters:
* path_to_input_colony_file
* output_folder_location
* sex_vector

### COL99-to-KIN-KINdot_SexAuto.R
Input file:
* txt COLONY99 file

Output files:
* txt KINGROUP file
* txt KINGROUPdot file
* csv table with heterozygosity diagnostics of samples
* csv table with heterozygosity diagnostics of SNPs

Required parameters:
* path_to_input_colony99_file
* output_folder_location

### COL99-to-KIN-KINdot_SexAsVector.R
Input file:
* txt COLONY99 file

Output files:
* txt KINGROUP file
* txt KINGROUPdot file
* csv table with heterozygosity diagnostics of samples
* csv table with heterozygosity diagnostics of SNPs

Required parameters:
* path_to_input_colony99_file
* output_folder_location
* sex_vector

### Explained list of all parameters:
* path_to_input_XXXXX_file - path to the location of the input file in specific txt format (XXXXX stands for GENEPOP, COLONY or COLONY99) in your computer (e.g. "C:\\Users\\Maxerová\\Desktop\\Misa\\genepop_example.txt")
* output_folder_location - path to the place where output folder (named ResultsFolder_DATE_TIME) will be saved (e.g. "C:\\Users\\Maxerová\\Desktop\\Misa\\")
* sex_vector - a vector containing information about sex of all the individuals in your dataset, from first to last sample, f for female, m for male (e.g. for dataset with 6 individuals, it could be c(f,f,m,m,f)
* n_initial_empty_rows - number of empty rows at the beginning of the genepop file (usually 1)
* n_snps - number of SNPs in your dataset (e.g. 6)
* n_rubbish_rows - number of empty or "rubbish" (containing some uninformative text) rows between SNPs names and the actual data (usually 1 - a single row containing "Pop")
