from pdf_preproccessing import preprocess_text, extract_section, extract_text, extract_abstract
from search_academic_databases import search_academic_databases, download_pdf_from_link
import pdfplumber

query = "how to conduct bibliometric analysis"
# get the result from arXiv,Google Scholar and PubMed.
result =search_academic_databases(query, 3)
#donwload the article to pdf
abstracts=[]
for res in result:
    try:
        if res['source_data'] == 'GoogleScholar':
            if  res['isFree'] == bool(True):
              content = download_pdf_from_link(res['link'],res['title'],query,res['source_data'])
              if content["isSuccess"] == True:
                  text = extract_text(content["content"])
                  token = preprocess_text(text)
                  # section = extract_section(content["content"])
                  section = extract_abstract(text)
                  if section != None:
                      abstracts.append({"text": section, "title": res["title"]})
                  
        else:
            content = download_pdf_from_link(res['link'], res['title'], query, res['source_data'])
            if content["isSuccess"] == True:
                text =extract_text(content["content"])
                token = preprocess_text(text)
                # section = extract_section(content["content"])
                section =extract_abstract(text)
                if section != None:
                    abstracts.append({"text":section ,"title":res["title"]})
    except:
        print("error "+res["title"])
print("finish")


