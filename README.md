# ncbi_genome_downloader
A tool to download large numbers of genomes from NCBI. 

USAGE:
  (1) Cell A1 of the csv file must contain the text '&&head'. This allows the program to remove the
      file encoding text when it parses the csv file.
  (2) Latin species names should begin on row 2. Column A contains species names, while column B contains 
      gene names. See included example file titled 'species.csv'.
  (3) The program will output two files into the current working directory (the directory the python
      file is contained in). One is a fasta file containing the genes, titled 'fastaOutput.txt'. The
      other is a csv file containing all the gene records, titled 'record.csv'.
  (4) Users can obtain an API key by signing up for an NCBI account, then clicking on their account
      name in the upper right-hand corner, and scrolling down. More information can be found here:
      https://ncbiinsights.ncbi.nlm.nih.gov/2017/11/02/new-api-keys-for-the-e-utilities/

NOTES:
  (1) It is slow. The server can only accept 3 requests per second, so I've put half second pauses
      after every request to ensure continuity. If someone is willing to negotiate with NCBI for an API
      with a higher request rate, we can use multi-threading and speed things up a bit.
  (2) I wrote this on a mac and never tested it on windows. There may be file path irregularities or other
      unforeseen consequences.

WARNINGS:
  (1) Double check number of exons. The records are not consistent.
  (2) Check missed results by hand. Any number of various syntax differences can cause a false return.
  (3) If the protein length has a fraction, check your genome by hand. Something is wrong.

FUTURE PLANS:
  (1) Find out if there's a way to add Ensembl support.
  (2) Add a quiet mode flag.
  (3) Find out why exon counts are not consistent; find a different way to get the information.
  (4) Deal with unusual characters returns in Get_GenBank(). Currently returns that cannot be type 
      coerced to int return a blank join result.

