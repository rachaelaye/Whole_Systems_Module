"""
C1 Select a gene and perform a thorough literature search

Write a script to obtain a list of papers from PubMed
Create a python script using Bio.Entrez
Gene: STK36
Bio.Entrez: https://biopython.org/docs/1.75/api/Bio.Entrez.html#:~:text=All%20the%20functions%20that%20send,failures%20or%20HTTP%205XX%20codes). 
"""

from Bio import Entrez
from config import email  # Import email from config.py
from xml.etree import ElementTree as ET

# Function to fetch PubMed papers related to the gene "SUFU"
def fetch_pubmed_papers(gene):
    # Set email
    Entrez.email = email

    # Search PubMed for the gene
    search_term = f"{gene}[Title/Abstract]"
    handle = Entrez.esearch(db="pubmed", term=search_term)
    record = Entrez.read(handle)
    handle.close()

    # Get list of PubMed IDs
    pubmed_ids = record["IdList"]

    # Fetch details for each PubMed ID
    papers = []
    for pmid in pubmed_ids:
        handle = Entrez.efetch(db="pubmed", id=pmid, retmode="xml")
        article_xml = handle.read()
        handle.close()

        # Parse XML
        root = ET.fromstring(article_xml)
        
        # Extract relevant information
        article = {}
        article["PMID"] = pmid
        article["Title"] = root.find(".//ArticleTitle").text
        article["Journal"] = root.find(".//Title").text
        article["Authors"] = [author.find(".//LastName").text + " " + author.find(".//ForeName").text for author in root.findall(".//Author")]
        article["Year"] = root.find(".//PubDate//Year").text
        
        # Get page range
        page_element = root.find(".//Pagination//MedlinePgn")
        if page_element is not None:
            pages = page_element.text.split("-")
            if len(pages) == 2:
                article["Pages"] = f"{pages[0]}-{pages[1]}"
            else:
                article["Pages"] = pages[0]
        else:
            article["Pages"] = "N/A"

        papers.append(article)

    return papers

if __name__ == "__main__":
    gene = "STK36"

    papers = fetch_pubmed_papers(gene)

    # Save the results to a file
    with open(f"{gene}_pubmed_papers.txt", "w") as file:
        for paper in papers:
            file.write(f"PMID: {paper['PMID']}\n")
            file.write(f"Title: {paper['Title']}\n")
            file.write(f"Authors: {', '.join(paper['Authors'])}\n")
            file.write(f"Journal: {paper['Journal']}\n")
            file.write(f"Year: {paper['Year']}\n")
            file.write(f"Pages: {paper['Pages']}\n\n")
