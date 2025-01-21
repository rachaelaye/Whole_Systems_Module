"""
C3: 

Entrez Links Descriptions: https://www.ncbi.nlm.nih.gov/entrez/query/static/entrezlinks.html 
Entrez Help: https://www.ncbi.nlm.nih.gov/books/NBK3837/#EntrezHelp.The_Entrez_Databases

"""

from Bio import Entrez
import csv
# Import email from config.py
from config import email  

# Define the search term and database
search_term = "STK36"
db = "gene"

# Set the email address to comply with NCBI's terms of use
Entrez.email = email

try:
    # Step 1: Perform the search to get the STK36 gene ID
    handle = Entrez.esearch(db=db, term=f"{search_term}[Gene Name]", retmode="xml")
    results = Entrez.read(handle)
    handle.close()

    # Extract the STK36 gene ID
    stk36_gene_id = results["IdList"][0] if results["IdList"] else None

    if stk36_gene_id:
        # Step 2: Use elink to find interactions with the pcassay_pcassay_gene_interaction_list linkname
        interacting_gene_ids = []

        handle = Entrez.elink(dbfrom="gene", db="gene", id=stk36_gene_id, linkname="pcassay_pcassay_gene_interaction_list")
        link_results = Entrez.read(handle)
        handle.close()

        for linkset in link_results:
            if "LinkSetDb" in linkset and linkset["LinkSetDb"]:
                for link in linkset["LinkSetDb"]:
                    for link_id in link["Link"]:
                        interacting_gene_ids.append(link_id["Id"])

        # Step 3: Fetch gene symbols for the interacting gene IDs
        if interacting_gene_ids:
            handle = Entrez.esummary(db="gene", id=",".join(interacting_gene_ids), retmode="xml")
            summary_results = Entrez.read(handle)
            handle.close()

            # Extract gene symbols
            interacting_gene_symbols = []
            for docsum in summary_results["DocumentSummarySet"]["DocumentSummary"]:
                interacting_gene_symbols.append(docsum["Name"])

            # Print the interacting genes
            print(f"Interacting genes for {search_term}:")
            for gene_symbol in interacting_gene_symbols:
                print(gene_symbol)
        else:
            print(f"No interacting genes found for {search_term}.")
    else:
        print(f"No gene found for {search_term}.")

except Exception as e:
    print(f"Error occurred: {e}")


"""
Creata a list of genes from the exported TSV from STRING Database
"""
# Step 1: Define the filename of your TSV file
input_filename = 'string_interactions.tsv'  # Replace with your actual file name
output_filename = 'STRING_interacting_genes.txt'

# Step 2: Initialize an empty set to store unique node values
unique_nodes = set()

# Step 3: Read the TSV file and extract node1 and node2 values
with open(input_filename, mode='r', newline='', encoding='utf-8') as file:
    reader = csv.DictReader(file, delimiter='\t')
    for row in reader:
        node1 = row['#node1']
        node2 = row['node2']
        unique_nodes.add(node1)
        unique_nodes.add(node2)

# Step 4: Convert the set to a sorted list (if needed)
unique_nodes_list = sorted(list(unique_nodes))

# Step 5: Save the unique node values to a new text file
with open(output_filename, 'w', encoding='utf-8') as outfile:
    for node in unique_nodes_list:
        outfile.write(node + '\n')

print(f"Unique nodes saved to '{output_filename}'")


"""
Compare the genes from Bio.Entrez and STRING
"""
def read_list_from_file(filename):
    """Reads a list of strings from a file."""
    with open(filename, 'r', encoding='utf-8') as file:
        return [line.strip() for line in file]

# Step 1: Define the filenames
file1 = 'STRING_interacting_genes.txt'
file2 = 'interacting_genes_list.txt' # List from Bio.Entrez database
output_file = 'comparison_results.txt'

# Step 2: Read lists from both files
list1 = read_list_from_file(file1)
list2 = read_list_from_file(file2)

# Step 3: Compare the lists
common_elements = set(list1).intersection(set(list2))
unique_to_list1 = set(list1) - set(list2)
unique_to_list2 = set(list2) - set(list1)

# Step 4: Save the comparison results to a file
with open(output_file, 'w', encoding='utf-8') as file:
    file.write(f"Common elements:\n")
    for element in sorted(common_elements):
        file.write(f"{element}\n")
    file.write(f"\nElements unique to {file1}:\n")
    for element in sorted(unique_to_list1):
        file.write(f"{element}\n")
    file.write(f"\nElements unique to {file2}:\n")
    for element in sorted(unique_to_list2):
        file.write(f"{element}\n")

print(f"Comparison results saved to '{output_file}'")


"""

from Bio import Entrez
from config import email  # Import email from config.py

# Define the search term and database
search_term = "STK36"
db = "gene"

# Set the email address to comply with NCBI's terms of use
Entrez.email = email

# Step 1: Perform the search to get the STK36 gene ID
handle = Entrez.esearch(db=db, term=f"{search_term}[Gene Name]", retmode="xml")
results = Entrez.read(handle)
handle.close()

# Extract the STK36 gene ID
stk36_gene_id = results["IdList"][0] if results["IdList"] else None
print(f"STK36 Gene ID: {stk36_gene_id}")

if stk36_gene_id:
    # Step 2: Use elink to find all types of links for STK36
    handle = Entrez.elink(dbfrom="gene", db="gene", id=stk36_gene_id)
    link_results = Entrez.read(handle)
    handle.close()

    # Debug print to check the structure of link_results and available links
    print("Link Results:", link_results)

    # Extract all available linknames
    available_linknames = set()
    for linkset in link_results:
        if "LinkSetDb" in linkset and linkset["LinkSetDb"]:
            for link in linkset["LinkSetDb"]:
                available_linknames.add(link["LinkName"])

    print("Available Linknames:", available_linknames)

else:
    print(f"No gene found for {search_term}.")

"""