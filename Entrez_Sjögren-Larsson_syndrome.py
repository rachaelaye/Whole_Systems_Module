import os
from Bio import Entrez
from Bio.Entrez import efetch, esearch, esummary, read
from config import email
from xml.etree import ElementTree as ET
import requests
from bs4 import BeautifulSoup

# Set email for NCBI Entrez
Entrez.email = email

# Directory to save the output files
output_dir = "sjogren_larsson_syndrome_data"
os.makedirs(output_dir, exist_ok=True)

# Function to fetch the OMIM database ID for a given condition
def fetch_omim_id(condition):
    search_term = f"{condition}[Title]"
    search_handle = esearch(db="omim", term=search_term, retmax=1)
    search_results = read(search_handle)
    search_handle.close()
    if search_results["IdList"]:
        return search_results["IdList"][0]
    return None

# Function to fetch PubMed papers linked to an OMIM entry
def fetch_pubmed_papers_linked_to_omim(omim_id):
    url = f"https://www.omim.org/entry/{omim_id}"
    headers = {
        "User-Agent": "Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/91.0.4472.124 Safari/537.36"
    }
    response = requests.get(url, headers=headers)
    response.raise_for_status()
    
    soup = BeautifulSoup(response.text, 'html.parser')
    references_section = soup.find('strong', text="REFERENCES")

    if not references_section:
        return []
    
    ol_tag = references_section.find_next('ol')

    if not ol_tag:
        return []
    
    references_content = ""
    for tag in ol_tag.find_all_next(['ol', 'div']):
        references_content += str(tag)
        if tag.name == 'div' and '<br></br>' in str(tag):
            break

    references_soup = BeautifulSoup(references_content, 'html.parser')

    references = []
    for p_tag in references_soup.find_all('p'):
        reference_text = p_tag.get_text(separator=" ", strip=True)
        reference_text = " ".join(reference_text.split())  # Replace extra spaces with single space

        # Stop copying over text when encountering [PubMed:
        if '[PubMed:' in reference_text:
            reference_text = reference_text.split('[PubMed:')[0].strip()
            
        pubmed_id = None
        pubmed_link = p_tag.find('a', href=lambda href: href and 'pubmed' in href)
        if pubmed_link:
            pubmed_id = pubmed_link.get_text(strip=True)
        
        references.append({
            'text': reference_text,
            'pubmed_id': pubmed_id
        })

    return references

# Function to write references to a file
def write_references_to_file(file_path, references):
    with open(file_path, 'w') as file:
        for reference in references:
            file.write(f"Paper: {reference['text']}\n")
            if reference['pubmed_id']:
                file.write(f"PubMed ID: {reference['pubmed_id']}\n")
            file.write("\n")

# Function to fetch data from Entrez
def fetch_entrez_data(db, id):
    try:
        search_handle = esearch(db=db, term=id, retmax=500)
        search_record = read(search_handle)
        search_handle.close()
        return search_record["IdList"]
    except Exception as e:
        print(f"Error fetching data for {id} from {db}: {e}")
        return None

# Function to fetch all PubMed papers mentioning the condition
def fetch_all_pubmed_papers(pubmed_id):
    # Get list of PubMed IDs
    pubmed_ids = pubmed_id

    # Fetch details for each PubMed ID
    papers = []
    for pmid in pubmed_ids:
        handle = efetch(db="pubmed", id=pmid, retmode="xml")
        article_xml = handle.read()
        handle.close()

        # Parse XML
        root = ET.fromstring(article_xml)

        # Extract relevant information
        article = {}
        article["PMID"] = pmid
        article["Title"] = root.find(".//ArticleTitle").text
        article["Journal"] = root.find(".//Title").text

        papers.append(article)

    return papers

# Function to fetch ClinVar entries related to the condition
def fetch_clinvar_entries(clinvar_ids):
    # Get list of ClinVar IDs
    clinvar_ids = clinvar_ids

    # Fetch details for each PubMed ID
    entries = []
    for clinvar in clinvar_ids:
        handle = esummary(db="clinvar", id=clinvar, retmode="xml")
        clinvar_xml = handle.read()
        handle.close()

        # Parse XML
        root = ET.fromstring(clinvar_xml)
     
        # Extract relevant information
        clinvar_entries = {}
        clinvar_entries["ClinVar ID"] = clinvar
        clinvar_entries["variation_name"] = root.find(".//variation_name").text
        clinvar_entries["description"] = root.find(".//description").text

        entries.append(clinvar_entries)

    return entries

# Main function to gather all required data
def main():
    condition = "Sjögren-Larsson syndrome"

    # Fetch OMIM ID
    omim_id = fetch_omim_id(condition)
    if omim_id:
        with open(os.path.join(output_dir, "omim_id.txt"), 'w') as file:
            file.write(f"OMIM ID for {condition}: {omim_id}\n")
        print(f"OMIM ID for {condition}: {omim_id}")
    else:
        print(f"OMIM ID for {condition} not found.")
        return

    # Fetch PubMed papers linked to OMIM entry
    pubmed_linked_papers = fetch_pubmed_papers_linked_to_omim(omim_id)
    if pubmed_linked_papers:
        write_references_to_file(os.path.join(output_dir, "pubmed_linked_to_omim.txt"), pubmed_linked_papers)
        print(f"PubMed Papers linked to OMIM entry {omim_id} written to file.")
    else:
        print(f"No PubMed papers found in OMIM entry {omim_id}.")

    # Get all PubMed IDs for PubMed papers linked to condition
    pubmed_ids = fetch_entrez_data("PubMed", condition)


    # Fetch all PubMed papers mentioning the condition
    all_pubmed_papers = fetch_all_pubmed_papers(pubmed_ids)
    with open(os.path.join(output_dir, "all_pubmed_papers.txt"), 'w') as file:
        for papers in all_pubmed_papers:
            file.write(f"PMID: {papers['PMID']}\n")
            file.write(f"Title: {papers['Title']}\n")
            file.write(f"Journal: {papers['Journal']}\n")
            file.write("\n")
    print(f"Fetched {len(all_pubmed_papers)} PubMed papers mentioning {condition}.")

    # Get all PubMed IDs for PubMed papers linked to condition
    clinvar_ids = fetch_entrez_data("ClinVar", condition)

    # Fetch ClinVar entries related to the condition
    clinvar_entries = fetch_clinvar_entries(clinvar_ids)
    with open(os.path.join(output_dir, "clinvar_entries.txt"), 'w') as file:
        for entry_id in clinvar_entries:
            file.write(f"ClinVar ID: {entry_id['ClinVar ID']}\n")
            file.write(f"Variant: {entry_id['variation_name']}\n")
            file.write(f"Clinical Significance: {entry_id['description']}\n")
            file.write("\n")
    print(f"Fetched {len(clinvar_entries)} ClinVar entries related to {condition}.")

if __name__ == "__main__":
    main()
