import csv
import os
import time
import requests
import xml.etree.ElementTree as ET


class org:
    def __init__(self, latinName, geneName):
        self.found = True
        self.latinName = latinName
        self.commonName = ""
        self.geneName = geneName
        self.geneID = ""
        self.description = ""
        self.accession = ""
        self.version = ""
        self.geneType = ""
        self.label = ""
        self.coverageDepth = ""
        self.startGene = ""
        self.endGene = ""
        self.strand = ""
        self.strandNum = ""
        self.exonCount = ""
        self.geneLength = ""
        self.proteinLength = ""
        self.exonLengths = []
        self.CDS = ""
        self.CDSlength = ""

    def printOrg(self):
        print("\nLatin name:      " + self.latinName)
        print("Common name:     " + self.commonName)
        print("Gene Name:       " + self.geneName)
        print("GeneID:          " + self.geneID)
        print("Description:     " + self.description)
        print("Gene Type:       " + self.geneType)
        print("Accession:       " + self.accession + "." + self.version)
        print("Coverage Depth:  " + self.coverageDepth)
        print("Exons:           " + self.exonCount)
        print("Exon Lengths:    " + str(self.exonLengths))
        print("Gene Label:      " + self.label)
        print("Strand:          " + self.strand)
        print("Gene start:      " + self.startGene)
        print("Gene end:        " + self.endGene)
        print("Gene length:     " + self.geneLength)
        print("Protein length   " + self.proteinLength)
        print("CDS length:      " + self.CDSlength)
        print("CDS:\n" + self.CDS)
        print("\n---------------------------------------------------------------------\n")


def Get_Input():
    csvLoc = input("Location of the csv file:\n")
    apiKey = input("API key:\n")
    csvLoc = csvLoc.strip()
    apiKey = apiKey.strip()
    return csvLoc, apiKey


def Import_CSV(csvLoc):
    # The first line the species file must be "--head". Otherwise python includes the file encoding. This gives us means to remove that
    species = {}
    with open('species.csv', newline='') as csvfile:
        thing = csv.reader(csvfile, delimiter=",")
        for i in thing:
            species.update({i[0] + " - " + i[1]: org(i[0], i[1].upper())})
        species.pop('\ufeff&&head - ')
    return species


def Get_Gene_Record(species):
    # Get the gene ID - it's buried in an xml file under <id>

    for organism in species:
        latinName = species.get(organism).latinName.replace(" ", "+")

        x = requests.get(
            "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=gene&term=" + latinName + "[ORGN]+AND+" + species.get(
                organism).geneName + "[GENE]" + "&retmode=xml" + "&api_key=" + apiKey)

        if x.text == "":
            print("Species " + species.get(organism).latinName + " not found.")
            continue

        # To not upset the NCBI server
        time.sleep(.5)

        # Get geneID
        rootEl = ET.fromstring(x.text)  # returns an Element object (not an ElementTree)

        try:
            species.get(organism).geneID = rootEl.find('.//Id').text  # .// searches all children at all levels
        except AttributeError:
            print("Species Record " + species.get(organism).geneName + " for " + species.get(
                organism).latinName + " not found.")
            species.get(organism).found = False
            continue

        # Get the gene record using the geneID
        x = requests.get('https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=gene&id=' + species.get(
            organism).geneID + "&retmode=xml" + "&api_key=" + apiKey)

        # To not upset the NCBI server even more
        time.sleep(.5)

        # new tree
        rootEl = ET.fromstring(x.text)

        # extract the common name
        try:
            species.get(organism).commonName = rootEl.find('./Entrezgene/Entrezgene_source//Org-ref_common').text
        except AttributeError:
            print("Gene " + species.get(organism).geneName + " for " + species.get(organism).latinName + " not found.")
            species.get(organism).found = False
            continue
        # extract the gene type
        species.get(organism).geneType = rootEl.find('./Entrezgene/Entrezgene_type').attrib['value']

        # extract accession number
        species.get(organism).accession = rootEl.find(
            './Entrezgene/Entrezgene_locus/Gene-commentary/Gene-commentary_accession').text

        # Version
        species.get(organism).version = rootEl.find(
            './Entrezgene/Entrezgene_locus/Gene-commentary/Gene-commentary_version').text

        # Extract the description
        species.get(organism).description = rootEl.find('./Entrezgene/Entrezgene_prot/Prot-ref/Prot-ref_desc').text

        # Extract exon count
        species.get(organism).exonCount = rootEl.find('./Entrezgene/Entrezgene_properties//Gene-commentary_text').text

        # Extract label
        species.get(organism).label = rootEl.find(
            './Entrezgene/Entrezgene_locus/Gene-commentary/Gene-commentary_label').text

        # Extract seq stop/start & length
        species.get(organism).startGene = rootEl.find('./Entrezgene//Gene-commentary_seqs//Seq-interval_from').text
        species.get(organism).endGene = rootEl.find('./Entrezgene//Gene-commentary_seqs//Seq-interval_to').text
        species.get(organism).geneLength = str(
            int(species.get(organism).endGene) - int(species.get(organism).startGene))

        # Get strand
        species.get(organism).strand = rootEl.find('./Entrezgene//Gene-commentary_seqs//Na-strand').attrib['value']

        # Stand conversion between search entrez search term and english (Gonna need this in the next search)
        if species.get(organism).strand == 'plus':
            species.get(organism).strandNum = '1'
        elif species.get(organism).strand == 'minus':
            species.get(organism).strandNum = '2'

    return (species)


def Get_CDS(species):
    def Process_Join(join):
        # take the join data out and make it pretty so we can put it in our excel file
        join = join.replace('..', ',')
        join = join.split(',')
        for i in join:
            try:
                join[join.index(i)] = int(i)
            except ValueError:
                join = []
                break

        exon = []
        count = 0
        while count < len(join):
            exon.append(join[count + 1] - join[count])
            count += 2

        return exon

    # get the giant fasta file
    for organism in species:
        if species.get(organism).found:
            theGreatFastaFile = requests.get(
                'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=' + species.get(
                    organism).accession + '&strand=' + species.get(organism).strandNum + '&seq_start=' + species.get(
                    organism).startGene + '&seq_stop' + species.get(
                    organism).endGene + '&rettype=fasta_cds_na&api_key=' + apiKey)
            # don't piss the server off (yikes)
            time.sleep(.5)
            theGreatFastaFile = theGreatFastaFile.text
            # cut it up
            begin = theGreatFastaFile.find("[db_xref=GeneID:" + species.get(organism).geneID + "]")
            theGreatFastaFile = theGreatFastaFile[begin:len(theGreatFastaFile)]

            # pause to get join data out
            beginJoin = theGreatFastaFile.find('join(', begin)
            join = theGreatFastaFile[beginJoin + len('join('): theGreatFastaFile.find(')', beginJoin)]
            species.get(organism).exonLengths = Process_Join(join)

            # continue chopping
            begin = theGreatFastaFile.find("\n", begin)
            end = theGreatFastaFile.find(">")  # The next record
            theGreatFastaFile = theGreatFastaFile[begin:end]

            CDS = theGreatFastaFile.replace("\n", "")
            CDS = CDS.strip()

            species.get(organism).CDS = CDS
    return species


def Process_CDS(species):
    for organism in species:
        if species.get(organism).found:

            # everything starts with a start codon
            begin = species.get(organism).CDS.find("ATG")

            # and ends with a stop codon
            TAA = species.get(organism).CDS.rfind("TAA")
            TAG = species.get(organism).CDS.rfind("TAG")
            TGA = species.get(organism).CDS.rfind("TGA")
            end = max(TAA, TAG, TGA) + 3

            # cut it out
            species.get(organism).CDS = species.get(organism).CDS[begin: end]

            # How long
            species.get(organism).CDSlength = str(len(species.get(organism).CDS))

            # How many amino acids
            species.get(organism).proteinLength = str(len(species.get(organism).CDS) / 3)

            # Reformat to 60 char per line to match Ensembl data
            b = 0
            e = 60
            formattedCDS = ""
            while b <= len(species.get(organism).CDS):
                formattedCDS = formattedCDS + species.get(organism).CDS[b:e] + "\n"
                b = b + 60
                e = e + 60

            # add the beginnings and ends
            species.get(organism).CDS = "> " + species.get(organism).geneName + "_" + species.get(
                organism).latinName.replace(' ', '_') + "_" + species.get(
                organism).accession + "\n" + formattedCDS + "\n"
    return species


def Get_GenBank(species):
    for organism in species:
        if species.get(organism).found:
            genBank = requests.get(
                'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=' + species.get(
                    organism).accession + '&rettype=gb&retmode=xml&api_key=' + apiKey)
            # don't anger the server gods
            time.sleep(.5)
            rootEl = ET.fromstring(genBank.text)
            uglyThing = rootEl.find('.//GBSeq_comment').text
            species.get(organism).coverageDepth = uglyThing[uglyThing.find('Genome Coverage :: ') + len(
                'Genome Coverage :: '): uglyThing.find(';', uglyThing.find('Genome Coverage :: ')) - 1]

    return species


def Create_Fasta(species):
    cwd = __file__[0:- len('ncbi_genomes.py')]

    # start clean
    if os.path.exists(cwd + "fastaOutput.txt"):
        os.remove(cwd + "fastaOutput.txt")

    fasta = open(cwd + "fastaOutput.txt", "a")

    for organism in species:
        if species.get(organism).found:
            data = species.get(organism).CDS
            fasta.write(data)
            species.get(organism).printOrg()

    fasta.close()


def Create_CSV(species):
    cwd = __file__[0:- len('ncbi_genomes.py')]

    # start clean
    if os.path.exists(cwd + "record.csv"):
        os.remove(cwd + "record.csv")

    f = open(cwd + 'record.csv', 'w', newline='')
    csvfile = csv.writer(f, dialect='excel')

    # header row
    csvfile.writerow([
        'Species',
        'Common Name',
        'Download Source',
        'Gene',
        'Gene Description (Ensembl)',
        'Transcript Access # (Ensembl)',
        'Gene Access # (Ensembl)',
        'Coverage Depth (NCBI)',
        'Gene Description (NCBI)',
        'Gene Type',
        'Accession.Version',
        'GeneID',
        'Exons',
        'Gene Length (bp)',
        'Exon Lengths',
        'CDS Length (bp)',
        'Protein Length (aa)',
        'Location',
        'Start',
        'Stop',
        'Strand'])

    for organism in species:
        if species.get(organism).found:
            csvfile.writerow([
                species.get(organism).latinName,
                species.get(organism).commonName,
                'NCBI',
                species.get(organism).geneName,
                "N/A",
                "N/A",
                "N/A",
                species.get(organism).coverageDepth,
                species.get(organism).description,
                species.get(organism).geneType,
                species.get(organism).accession + '.' + species.get(organism).version,
                species.get(organism).geneID,
                species.get(organism).exonCount,
                species.get(organism).geneLength,
                species.get(organism).exonLengths,
                species.get(organism).CDSlength,
                species.get(organism).proteinLength,
                species.get(organism).label,
                species.get(organism).startGene,
                species.get(organism).endGene,
                species.get(organism).strand])
        else:
            csvfile.writerow([
                species.get(organism).latinName, "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "",
                "", ""])
    f.close()


def Main():
    global apiKey
    csvLoc, apiKey = Get_Input()
    species = Import_CSV(csvLoc)
    species = Get_Gene_Record(species)
    species = Get_CDS(species)
    species = Process_CDS(species)
    species = Get_GenBank(species)
    Create_Fasta(species)
    Create_CSV(species)


Main()
