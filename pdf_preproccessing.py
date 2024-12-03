import pikepdf
from PyPDF2 import PdfReader
import os

from PyPDF2.errors import PdfReadError



def read_pdf_file(path):
    if validate_file_path(path):
        print(path)
        try:
            reader = PdfReader(path)
            text = "".join([page.extract_text() for page in reader.pages])
            paragraphs = text.split('\n\n')
            return paragraphs
        except PdfReadError as e:
            print(f"Error reading PDF: {path} , {e}")
            repair_pdf(path,path)
            isRunFirstTime = bool(True)
            if isRunFirstTime :
                isRunFirstTime = bool(False)
                read_pdf_file(path)
    else:
        return None

def validate_file_path(file_path):
    if os.path.exists(file_path) and os.path.isfile(file_path):
        return True
    return False

def repair_pdf(input_path, output_path):
    try:
        pdf = pikepdf.open(input_path)
        pdf.save(output_path)
        print(f"Repaired PDF saved to: {output_path}")
    except pikepdf.PdfError as e:
        print(f"Failed to repair PDF: {e}")