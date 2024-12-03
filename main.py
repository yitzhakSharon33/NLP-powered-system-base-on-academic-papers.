from pdf_preproccessing import read_pdf_file
from search_academic_databases import search_academic_databases, download_pdf_from_link

query = "deep learning in healthcare"
# get the result from arXiv,Google Scholar and PubMed.
result =search_academic_databases(query, 3)
#donwload the article to pdf
for res in result:
    if res['source_data'] == 'GoogleScholar':
        if  res['isFree'] == bool(True):
          path = download_pdf_from_link(res['link'],res['title'],query,res['source_data'])
          test = read_pdf_file(path)
            
    else:
        path = download_pdf_from_link(res['link'], res['title'], query, res['source_data'])
        test = read_pdf_file(path)
        

        


