import os

from bs4 import BeautifulSoup

from API.arXiv import search_arxiv
from API.googleScholar import search_google_scholar
from API.pubMed import search_pubmed
import requests
from io import BytesIO

def search_academic_databases(query, max_results=5):
    print(f"Searching for '{query}'...\n")
    
    # Search arXiv
    print("Results from arXiv:")
    arxiv_results = search_arxiv(query, max_results)
    parseResult=[]
    for result in arxiv_results:
        print(f"Title: {result['title']}\nPublished: {result['published']}\nLink: {result['link']}\n")
        parseResult.append({'title':result['title'],'Published':result['published'], 'link': result['link'],'source_data':'arXiv'})
    # Search PubMed
    print("\nResults from PubMed:")
    pubmed_results = search_pubmed(query, max_results)
    for result in pubmed_results:
        print(f"Title: {result['title']}\nSource: {result['source']}\nPubDate: {result['pub_date']}\nLink: {result['link']}\n")
        parseResult.append({'title':result['title'],'source':result['source'] ,'Published':result['pub_date'], 'link': result['link'],'source_data':'PubMed' })
    # Search Google Scholar
    print("\nResults from Google Scholar:")
    google_scholar_results = search_google_scholar(query, max_results)
    for result in google_scholar_results:
        print(f"Title: {result['title']}\nAuthor(s): {result['author']}\nYear: {result['year']}\nLink: {result['link']}\n")
        parseResult.append({'title':result['title'],'author':result['author'] ,'Published':result['year'], 'isFree':result['isFree'],'link': result['link'],'source_data':'GoogleScholar'})
    return parseResult

def download_pdf_from_link(link,title,query,source):
    directory= f"articles/{query}/{source}"
    if not os.path.exists(directory):
        os.makedirs(directory)
        print(f"Directory '{directory}' created successfully!")
    
    if title[-1]== '.':
        title=title[:-1]
    output_file = f"{directory}/{title}.pdf"
    print(f"\nstart download pdf from link  {link} , {output_file}")
    
    # Download the PDF
    response = requests.get(link)
    if response.status_code == 200:
        with open(output_file, "wb") as file:
            file.write(response.content)
        print(f"PDF downloaded successfully: {output_file}")
        return {"isSuccess":True,"content":response.content}
    else:
        print(f" ***** Failed to download PDF. Status code: {response.status_code} *****")
        return {"isSuccess":False, "content":response.content}

# query = "deep learning in healthcare"
# # get the result from arXiv,Google Scholar and PubMed.
# result =search_academic_databases(query, 3)
# #donwload the article to pdf
# for res in result:
#     if res['source_data'] == 'GoogleScholar':
#         if  res['isFree'] == bool(True):
#             download_pdf_from_link(res['link'],res['title'],query,res['source_data'])
#     else:
#         download_pdf_from_link(res['link'], res['title'], query, res['source_data'])


