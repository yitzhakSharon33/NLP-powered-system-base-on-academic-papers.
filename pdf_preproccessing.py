import re
import nltk
from nltk.tokenize import word_tokenize, sent_tokenize
from nltk.corpus import stopwords
from nltk.stem import PorterStemmer, WordNetLemmatizer
from sklearn.feature_extraction.text import TfidfVectorizer
from PyPDF2 import PdfReader
from io import BytesIO
import pdfplumber

# Download necessary NLTK data (only needs to be done once)
nltk.download('punkt')
nltk.download('punkt_tab')

nltk.download('stopwords')
nltk.download('wordnet')

# Text Preprocessing Pipeline

def text_cleaning(text):
    # Remove non-alphanumeric characters (punctuation, special symbols)
    text = re.sub(r'[^a-zA-Z\s]', '', text)
    # Lowercase the text
    text = text.lower()
    return text

def tokenization(text):
    # Tokenize the text into words
    tokens = word_tokenize(text)
    return tokens

def remove_stopwords(tokens):
    # Remove common stopwords
    stop_words = set(stopwords.words('english'))
    filtered_tokens = [word for word in tokens if word not in stop_words]
    return filtered_tokens

def stemming(tokens):
    # Apply stemming using PorterStemmer
    stemmer = PorterStemmer()
    stemmed_tokens = [stemmer.stem(word) for word in tokens]
    return stemmed_tokens

def lemmatization(tokens):
    # Apply lemmatization using WordNetLemmatizer
    lemmatizer = WordNetLemmatizer()
    lemmatized_tokens = [lemmatizer.lemmatize(word) for word in tokens]
    return lemmatized_tokens

def vectorization(corpus):
    # Vectorize text using TF-IDF
    vectorizer = TfidfVectorizer(max_features=100)
    X = vectorizer.fit_transform(corpus)
    return X, vectorizer.get_feature_names_out()


def contentPreProccsing(document):
    document = document.decode('utf-8', errors='ignore')
    # Preprocessing Steps
    cleaned_text = text_cleaning(document)
    #print(f"Cleaned Text: {cleaned_text}")
    
    tokens = tokenization(cleaned_text)
    print(f"Tokens: {len(tokens)}")
    
    filtered_tokens = remove_stopwords(tokens)
    print(f"Filtered Tokens (no stopwords): {len(filtered_tokens)}")
    
    # Uncomment one of the following lines based on your preference for stemming or lemmatization
    processed_tokens = stemming(filtered_tokens)
    # processed_tokens = lemmatization(filtered_tokens)
    
    print(f"Processed Tokens (after stemming or lemmatization): {len(processed_tokens)}")
    
    # If you'd like to perform vectorization, use the following
    corpus = [' '.join(processed_tokens)]  # Convert list of tokens back to string
    X, feature_names = vectorization(corpus)
    
    print("\nTF-IDF Features:")
    for feature, score in zip(feature_names, X.toarray()[0]):
        print(f"{feature}: {score:.4f}")


import re
import nltk
from nltk.tokenize import word_tokenize
from nltk.corpus import stopwords

# Download NLTK resources if necessary
nltk.download('punkt')
nltk.download('stopwords')

def extract_text(content):
    if type(content) is str:
        return content
    pdf_file = BytesIO(content)
    
    # Read the PDF from the in-memory file
    reader = PdfReader(pdf_file)
    text = ''
    # Extract text from each page
    for i, page in enumerate(reader.pages):
        text += page.extract_text()
    return text

def preprocess_text(text):
    
        
    # Remove non-alphabetic characters
    text = re.sub(r'[^A-Za-z\s]', '', text)
    
    # Tokenize text
    tokens = word_tokenize(text.lower())
    
    # Remove stopwords
    stop_words = set(stopwords.words('english'))
    filtered_tokens = [word for word in tokens if word not in stop_words]
    
    return filtered_tokens

def extract_section(content):
    pdf_file = BytesIO(content)
    
    # Read the PDF with pdfplumber
    with pdfplumber.open(pdf_file) as pdf:
        start_keyword = "abstract"
        end_keyword = "Introduction"
        extracting = False
        section_text = ""
        
        for page in pdf.pages:
            page_text = page.extract_text()
            if start_keyword in page_text:
                extracting = True
            if extracting:
                section_text += page_text
                section_text = section_text.split("abstract")[1]
            if extracting and end_keyword in page_text:
                break
    
    # Output the extracted section
    print("Extracted Section:")
    return (section_text)


import re

def extract_abstract(text):
    pattern = r"(?:Abstract|ABSTRACT|Summary|Introduction)\s*[:\n]([\s\S]*?)(?=\n[A-Z]{2,}|\n1\.|\nIntroduction)"
    match = re.search(pattern, text)
    if match:
        return match.group(1).strip()
    return None