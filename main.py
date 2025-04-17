from pdf_preproccessing import preprocess_text, extract_section, extract_text, extract_abstract
from search_academic_databases import search_academic_databases, download_pdf_from_link
import pdfplumber

query = "Deep Learning rejection"
word_to_check = "3D reconstruction"

# get the result from arXiv,Google Scholar and PubMed.
result =search_academic_databases(query, 10)
print(result)
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
                      res["text"] = section
                      abstracts.append(res)
        else :
            content = download_pdf_from_link(res['link'], res['title'], query, res['source_data'])
            if content["isSuccess"] == True:
                text =extract_text(content["content"])
                token = preprocess_text(text)
                # section = extract_section(content["content"])
                section =extract_abstract(text)
                if section != None:
                    res["text"]=section
                    abstracts.append(res)
    except:
        print("error "+res["title"])
print(len(abstracts))


import spacy
from spacy.lang.en import English

# Load the English NLP model
nlp = spacy.load("en_core_web_sm")


# Word to check
wordToCheck= word_to_check.split(' ')

# Process each abstract
for idx, abstract in enumerate(abstracts):
    doc = nlp(abstract['text'])
    lemmas = [token.lemma_.lower() for token in doc]
    if len(wordToCheck)>1:
        if wordToCheck[0].lower() in lemmas and wordToCheck[1].lower() in lemmas:
            print(f"✅ '{word_to_check}' found in abstract {abstract['title']}:")
            # print(f"    {abstract}")
        else:
            print(f"❌ '{word_to_check}' NOT found in abstract {abstract['title']}.")
    else:
        if word_to_check.lower() in lemmas:
            print(f"✅ '{word_to_check}' found in abstract {abstract['title']}:")
        else:
            print(f"❌ '{word_to_check}' NOT found in abstract {abstract['title']}.")
        
