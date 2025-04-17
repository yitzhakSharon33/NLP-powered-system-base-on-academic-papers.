from Bio import Entrez
from xml.etree import ElementTree
import requests
import json



def search_pubmed(query, max_results=5):
    Entrez.email = "your-email@example.com"  # NCBI requires an email for tracking usage
    handle = Entrez.esearch(db="pubmed", term=query, retmax=max_results)
    record = Entrez.read(handle)
    handle.close()
    
    # Fetch details for the articles found
    ids = record["IdList"]
    if ids:
        handle = Entrez.esummary(db="pubmed", id=",".join(ids))
        summaries = Entrez.read(handle)
        handle.close()
        print(summaries)
        results = []
        for summary in summaries:
            print(summary)
            results.append({
                'title': summary['Title'],
                'source': summary.get('Source'),
                'pub_date': summary.get('PubDate'),
                'link': f"https://pubmed.ncbi.nlm.nih.gov/pmc/articles/{summary.get('Id')}/pdf/",
                'authorList': summary.get('AuthorList')
            })
            print(f"https://www.ncbi.nlm.nih.gov/pmc/articles/{summary.get('Id')}/pdf/")
        return results
    else:
        print("No results found.")
        return []

search_pubmed("how to conduct bibliometric analysis",5)
