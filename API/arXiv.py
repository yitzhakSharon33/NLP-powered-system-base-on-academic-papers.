import requests
from docutils.nodes import author


def search_arxiv(query, max_results=5):
    base_url = "http://export.arxiv.org/api/query?"
    params = {
        'search_query': query,
        'start': 0,
        'max_results': max_results,
        'sortBy': 'submittedDate',
        'sortOrder': 'descending'
    }

    response = requests.get(base_url, params=params)

    if response.status_code == 200:
        # Parse the XML response (arXiv returns XML format)
        import xml.etree.ElementTree as ET
        tree = ET.ElementTree(ET.fromstring(response.text))
        root = tree.getroot()
        ns = {'atom': 'http://www.w3.org/2005/Atom'}

        # Iterate over each entry in the arXiv response
        entries = root.findall('{http://www.w3.org/2005/Atom}entry')
        results = []
        print('arxiv')
        print(entries)
        for entry in entries:
            title = entry.find('{http://www.w3.org/2005/Atom}title').text
            summary = entry.find('{http://www.w3.org/2005/Atom}summary').text
            published = entry.find('{http://www.w3.org/2005/Atom}published').text
            link = entry.find('{http://www.w3.org/2005/Atom}link').attrib['href']
            authors = entry.find('{http://www.w3.org/2005/Atom}author')
            categories = entry.findall('{http://www.w3.org/2005/Atom}:category')
            journal_ref =entry.find('{http://www.w3.org/2005/Atom}:journal_ref')
            results.append({
                'title': title,
                'summary': summary,
                'published': published,
                'link': link.replace('/abs/','/pdf/'),
                # 'author':[author.find('{http://www.w3.org/2005/Atom}:name', ).text for author in authors],
                'categories': [category.attrib['term'] for category in categories],
                'journal_ref':journal_ref,
                'source_data':'arxiv'
            })
    
        return results
    else:
        print(f"Error: {response.status_code}")
        return None