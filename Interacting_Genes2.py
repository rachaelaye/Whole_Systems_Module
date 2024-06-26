"""
C3: 

Extract from webpage.

The interactions data displayed on the NCBI Gene webpage for a specific gene (https://www.ncbi.nlm.nih.gov/gene/27148#interactions), does not get stored in a database that is accessible or queryable in the same manner as NCBI's primary databases like PubMed, Gene, or Protein.

Here's how the interactions data typically works on NCBI's Gene webpage:
Display of Interactions: The interactions data you see on the NCBI Gene webpage is curated and displayed by NCBI based on various sources of biological data. These sources can include experimental data, computational predictions, literature annotations, and curated databases.
Sources of Interaction Data: NCBI integrates interaction data from various databases and resources, including but not limited to BioGRID, IntAct, STRING, and curated literature. These interactions may involve physical interactions, genetic interactions, pathway associations, and other types of molecular interactions.
Data Presentation: NCBI aggregates and presents this interaction data in a structured format on the Gene webpage under the "Interactions" section. This includes information about interacting genes, proteins, and sometimes additional details such as experimental methods or confidence scores.
Dynamic Updates: The interaction data displayed on the NCBI Gene webpage is regularly updated as new information becomes available. NCBI aims to provide the most current and comprehensive data possible, reflecting ongoing research and database updates.
Accessing Interaction Data: While you can view and explore interaction data interactively on the NCBI Gene webpage, this specific presentation is not queryable in the same way as raw data in NCBI's primary databases. Instead, NCBI provides APIs (like Entrez utilities) and downloads (like FTP repositories) for accessing bulk data from sources that contribute to interaction data.

"""

import csv
import requests
from bs4 import BeautifulSoup

def fetch_interaction_data(base_url, params):
    """
    Fetches interaction data from the specified URL with given parameters.
    Returns BeautifulSoup object of the fetched HTML.
    """
    try:
        response = requests.get(base_url, params=params)
        response.raise_for_status()  # Raise exception for bad responses
        soup = BeautifulSoup(response.content, 'html.parser')
        return soup
    except requests.exceptions.RequestException as e:
        print(f"Error fetching data: {e}")
        return None

def parse_html_for_other_genes(soup):
    """
    Parses the HTML content to extract "Other Gene" information.
    Returns a list of "Other Gene" names.
    """
    other_genes = []

    # Locate all rows (<tr>) in the table body
    rows = soup.find_all('tr')
    
    for row in rows:
        # Find all cells (columns) in the row
        cells = row.find_all('td')
        
        if len(cells) >= 3:  # Assuming "Other Gene" is in the third column
            other_gene = cells[2].text.strip()
            other_genes.append(other_gene)
    return other_genes

def save_other_genes_to_file(other_genes, filename="other_genes.txt"):
    """
    Saves the list of "Other Gene" names to a file.
    """
    with open(filename, "w") as file:
        for gene in other_genes:
            file.write(gene + "\n")
    print(f"Other genes saved to '{filename}'")

def main():
    base_url = "https://www.ncbi.nlm.nih.gov/gene/27148/data/interaction/"
    params = {
        "p$l": "Ajax",
        "p$debugoutput": "off",
        "p$st": "entrez",
        "p$saveparams": "false",
        "pageize": "25",  # Adjust as per your pagination size
        "p$pg": 1  # Starting page number
    }

    other_genes = []
    while True:
        # Fetch interaction data
        soup = fetch_interaction_data(base_url, params)
        if soup:
            # Parse HTML for other genes
            other_genes.extend(parse_html_for_other_genes(soup))

            # Check if there's a next page
            next_page = soup.find('input', {'title': 'Next Page'})
            if not next_page:
                break  # Exit loop if no next page link is found

            # Update parameters for the next page
            params['p$pg'] += 1

        else:
            print("Error: Unable to fetch or parse data.")
            break

    # Save other genes to file
    save_other_genes_to_file(other_genes)

if __name__ == "__main__":
    main()


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
file2 = 'other_genes.txt' # List from Bio.Entrez database
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
import requests
from bs4 import BeautifulSoup

def fetch_and_save_data(url, filename):
    try:
        # Fetch the interaction data as HTML using requests
        response = requests.get(url)
        response.raise_for_status()  # Raise an exception for bad responses

        # Save the raw HTML content to a file
        with open(filename, "w", encoding="utf-8") as file:
            file.write(response.text)
        
        print(f"HTML content saved to '{filename}'")
    except requests.exceptions.RequestException as e:
        raise RuntimeError(f"Error fetching data: {e}")
    except Exception as ex:
        raise RuntimeError(f"Error: {ex}")

def parse_saved_html(filename):
    try:
        # Read the saved HTML file
        with open(filename, "r", encoding="utf-8") as file:
            html_content = file.read()

        # Parse HTML content with BeautifulSoup
        soup = BeautifulSoup(html_content, 'html.parser')

        # Extract "Other Gene" information from the parsed HTML
        other_genes = []

        # Locate all rows (<tr>) in the table body
        rows = soup.find_all('tr')
        
        for row in rows:
            # Find all cells (columns) in the row
            cells = row.find_all('td')
            
            if len(cells) >= 3:  # Assuming "Other Gene" is in the third column
                other_gene = cells[2].text.strip()
                other_genes.append(other_gene)

        # Save the "Other Gene" information to a file
        with open("other_genes_from_saved_html.txt", "w") as file:
            for gene in other_genes:
                file.write(gene + "\n")

        print(f"Other genes saved to 'other_genes_from_saved_html.txt'")

    except FileNotFoundError:
        raise RuntimeError(f"Error: Saved HTML file '{filename}' not found.")
    except Exception as ex:
        raise RuntimeError(f"Error: {ex}")

if __name__ == "__main__":
    # Define the URL for the interaction data service
    data_url = "https://www.ncbi.nlm.nih.gov/gene/27148/data/interaction/?p$l=Ajax&amp;p$debugoutput=off"
    
    # Define the filename to save the HTML content
    saved_html_filename = "saved_interaction_data.html"

    try:
        # Fetch and save HTML data from the URL
        fetch_and_save_data(data_url, saved_html_filename)

        # Parse saved HTML file and extract "Other Gene" information
        parse_saved_html(saved_html_filename)

    except RuntimeError as e:
        print(f"Error: {e}")
"""