"""
Links to relevant information for interacting with ClinVar and SNP databases:
- https://www.ncbi.nlm.nih.gov/clinvar/docs/maintenance_use/ 
- https://www.ncbi.nlm.nih.gov/books/NBK179288/ 
- https://www.ncbi.nlm.nih.gov/clinvar/docs/identifiers/ 
"""


import os
import requests
import xml.etree.ElementTree as ET
from config import email  # Import email from config.py
from Bio import Entrez
import re
import csv


# Set email for NCBI
Entrez.email = email


# Set email
Entrez.email = email

# Function to fetch data from Entrez
def fetch_entrez_data(db, gene):
    try:
        search_handle = Entrez.esearch(db=db, term=f"{gene}[gene]", retmax=500)
        search_record = Entrez.read(search_handle)
        search_handle.close()
        return search_record["IdList"]
    except Exception as e:
        print(f"Error fetching data for {gene} from {db}: {e}")
        return None


# Function to fetch ClinVar summaries for a list of IDs
def fetch_clinvar_summaries(ids):
    """
    Fetches ClinVar summaries for a list of IDs.
    :param ids: List of ClinVar IDs.
    :return: XML data as a string or None if an error occurs.
    """
    try:
        fetch_handle = Entrez.esummary(db="clinvar", id=",".join(ids), rettype="xml")
        summaries = fetch_handle.read()
        fetch_handle.close()
        return summaries
    except Exception as e:
        print(f"Error fetching ClinVar summaries: {e}")
        return None

# Function to parse ClinVar summaries
def parse_clinvar_summaries(xml_data, gene):
    snps = []
    root = ET.fromstring(xml_data)
    for docsum in root.findall(".//DocumentSummary"):
        db_id_elem = docsum.find(".//variation_xref[db_source='dbSNP']/db_id")
        if db_id_elem is None:
            continue  # Skip if dbSNP ID is not present
        
        snp = {
            'gene': gene,
            'id': docsum.findtext("accession", default="N/A"),
            'variation_name': docsum.findtext(".//variation_name", default="N/A"),
            'variant_type': docsum.findtext(".//variant_type", default="N/A"),
            'description': docsum.findtext(".//germline_classification/description", default="N/A"),
            'Phenotype': docsum.findtext(".//trait_set/trait/trait_name", default="N/A"),
            'db_id': db_id_elem.text
        }
        snps.append(snp)
    return snps

# Function to fetch dbSNP rs number for a ClinVar variant
def fetch_dbsnp_rs_for_clinvar(db_id):
    try:
        url = f"https://www.ncbi.nlm.nih.gov/snp/{db_id}"
        response = requests.get(url)
        response.raise_for_status()

        # Extract rs number from the HTML
        rs_numbers = re.findall(r'/snp/(rs\d+)', response.text)
        rs_numbers = list(set(rs_numbers))  # Remove duplicates
        return rs_numbers
    except Exception as e:
        print(f"Error fetching dbSNP rs for ClinVar variant {db_id}: {e}")
        return []


# Function to merge ClinVar and dbSNP data
def merge_clinvar_dbsnp(clinvar_snps, dbsnp_rs_numbers):
    merged_snps = []
    for clinvar_snp in clinvar_snps:
        for rs_number in dbsnp_rs_numbers:
            merged_snp = {
                'gene': clinvar_snp['gene'],
                'variation_name': clinvar_snp['variation_name'],
                'variant_type': clinvar_snp['variant_type'],
                'description': clinvar_snp['description'],
                'db_id': clinvar_snp['db_id'],
                'Phenotype': clinvar_snp['Phenotype'],
                'rs_id': rs_number
            }
            merged_snps.append(merged_snp)
    return merged_snps

# Function to write SNPs to file
def write_snps_to_file(gene, snps, output_dir):
    file_path = os.path.join(output_dir, f"{gene}_snps.csv")
    file_exists = os.path.isfile(file_path)
    with open(file_path, 'a', newline='') as file:
        writer = csv.writer(file)
        if not file_exists:
            writer.writerow(["Gene", "Variant", "Variant Type", "Clin_Sig", "Phenotype", "dbSNP ID"])
        for snp in snps:
            writer.writerow([snp['gene'], snp['variation_name'], snp['variant_type'], snp['description'], snp['Phenotype'], snp['rs_id']])

# Main function
def main():
    # Paths to gene list files
    gene_list_files = ["results_files/other_genes.txt", "results_files/STRING_interacting_genes.txt"]

    # Directory to save the output files
    output_dir = "snp_output"
    os.makedirs(output_dir, exist_ok=True)

    for gene_list_file in gene_list_files:
        with open(gene_list_file, 'r') as file:
            genes = [line.strip() for line in file.readlines()]

        for gene in genes:
            print(f"Fetching data for gene: {gene}")

            clinvar_ids = fetch_entrez_data("clinvar", gene)

            if not clinvar_ids:
                print(f"No ClinVar IDs found for gene: {gene}")
                continue

            clinvar_summaries = fetch_clinvar_summaries(clinvar_ids)

            merged_snps = []
            if clinvar_summaries:
                clinvar_snps = parse_clinvar_summaries(clinvar_summaries, gene)
                for clinvar_snp in clinvar_snps:
                    if clinvar_snp['db_id']:
                        dbsnp_rs_numbers = fetch_dbsnp_rs_for_clinvar(clinvar_snp['db_id'])
                        if dbsnp_rs_numbers:
                            merged_snps.extend(merge_clinvar_dbsnp([clinvar_snp], dbsnp_rs_numbers))
                        else:
                            print(f"No dbSNP rs numbers found for ClinVar ID {clinvar_snp['db_id']}")
                    else:
                        print(f"No dbSNP ID found for ClinVar variant {clinvar_snp['id']}")


            if merged_snps:
                write_snps_to_file(gene, merged_snps, output_dir)
                print(f"Data written for merged SNPs of gene: {gene}")
            else:
                print(f"No merged SNPs found for gene: {gene}")

if __name__ == "__main__":
    main()
