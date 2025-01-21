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
        "page": 1,
        "pageSize": "25",
    }

    other_genes = []
    while True:
        # Fetch interaction data
        soup = fetch_interaction_data(base_url, params)
        if soup:
            # Parse HTML for other genes
            other_genes.extend(parse_html_for_other_genes(soup))

            next_page_link = soup.find('a', {'class': 'ui-ncbigrid-paged-pageControl-next next page_link', 'title': 'Goto Next Page'})
            if next_page_link:
                next_page_url = next_page_link.get('href')
                if next_page_url:
                # Update parameters for the next page in the table
                    params["page"] += 1
                else:
                    break # Exit loop if no next page link is found

            else:
                print("Error: Unable to fetch or parse data.")
                break

    # Save other genes to file
    save_other_genes_to_file(other_genes)

if __name__ == "__main__":
    main()


import requests
import pandas as pd

def scrape_interactions_table(url):
    # Send a GET request to the URL
    response = requests.get(url)
    
    if response.status_code == 200:
        # Extract data from AJAX endpoint (example, adjust based on actual AJAX request)
        ajax_url = 'https://www.ncbi.nlm.nih.gov/gene/27148/data/interaction/?p$l=Ajax&p$debugoutput=off&page=3&pageSize=25'  # Replace with actual AJAX endpoint
        ajax_response = requests.get(ajax_url)
        
        if ajax_response.status_code == 200:
            # Save JSON data to a text file
            output_file = 'interactions_data.json'  # Adjust the file path as needed
            save_to_txt(ajax_response.text, output_file)
            
            return True
        else:
            print(f"Failed to fetch AJAX data. Status code: {ajax_response.status_code}")
            return False
    else:
        print(f"Failed to fetch page. Status code: {response.status_code}")
        return False

    
def save_to_txt(data, file_path):
    try:
        with open(file_path, 'w') as f:
            f.write(data)
        print(f"Data saved to {file_path}")
    except Exception as e:
        print(f"Error saving data to {file_path}: {str(e)}")


def main():
    base_url = 'https://www.ncbi.nlm.nih.gov/gene/27148/data/interaction/'
    
    print(f"Scraping data from {base_url}")
    
    if scrape_interactions_table(base_url):
        print("Data scraped and saved successfully.")
    else:
        print("Failed to scrape data.")

if __name__ == "__main__":
    main()


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