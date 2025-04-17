from docutils.nodes import author
from scholarly import scholarly
import json
def search_google_scholar(query, max_results=5):
    search_query = scholarly.search_pubs(query)
    results = []
    for _ in range(max_results):
        try:
            pub = next(search_query)
            print("search_google_scholar")
            print(json.dumps(pub))
            title =  pub['bib']['title']
            author = pub['bib']['author']
            year = pub['bib']['pub_year']
            citations = pub['citedby_url']
            num_citations = pub['num_citations']
            isFree = bool(False)
            if pub.get('eprint_url') is not None:
                link = pub['eprint_url']
                isFree = bool(True)
            else:
                link =  pub['pub_url']
        
            results.append({
                'title': title,
                'author':author,
                'year': year,
                'link':link,
                'isFree': isFree,
                'citations':citations,
                'num_citations':num_citations
            })
        except StopIteration:
            break
    return results
search_google_scholar('matching learning',2)

