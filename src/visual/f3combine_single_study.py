# Combine f3 single study figures into one figure
# Row 1: f3cor_qtdid + f3egene_upset_qtdid
# Row 2: f3neff_scatter_qtdid_annot
from utils import *
import pypdfium2 as pdfium
from PIL import Image
import io

target_qtdids = [
    "QTD000021",
    "QTD000031",
    "QTD000066",
    "QTD000067",
    "QTD000069",
    "QTD000073",
    "QTD000081",
    "QTD000115",
    "QTD000371",
    "QTD000372",
]

input_path = f"{save_path}/single_study"
output_path = f"{save_path}/single_study_combined"
os.makedirs(output_path, exist_ok=True)


def pdf_to_image(pdf_path, dpi=300):
    """Convert PDF to PIL Image"""
    pdf = pdfium.PdfDocument(pdf_path)
    page = pdf[0]
    # Render at high DPI for quality
    scale = dpi / 72  # PDF default is 72 DPI
    bitmap = page.render(scale=scale)
    pil_image = bitmap.to_pil()
    pdf.close()
    return pil_image


def combine_pdfs(qtdid):
    """Combine three PDFs into one figure for a given study"""
    # Define input PDF paths
    cor_pdf = f"{input_path}/f3cor_{qtdid}.pdf"
    egene_pdf = f"{input_path}/f3egene_upset_{qtdid}.pdf"
    neff_pdf = f"{input_path}/f3neff_scatter_{qtdid}_annot.pdf"

    # Check if all files exist
    for pdf_file in [cor_pdf, egene_pdf, neff_pdf]:
        if not os.path.exists(pdf_file):
            print(f"Warning: {pdf_file} not found, skipping {qtdid}")
            return

    # Convert PDFs to images
    cor_img = pdf_to_image(cor_pdf)
    egene_img = pdf_to_image(egene_pdf)
    neff_img = pdf_to_image(neff_pdf)

    # Get dimensions
    cor_w, cor_h = cor_img.size
    egene_w, egene_h = egene_img.size
    neff_w, neff_h = neff_img.size

    # Row 1: cor + egene (side by side)
    row1_w = cor_w + egene_w
    row1_h = max(cor_h, egene_h)

    # Row 2: neff (full width, centered)
    row2_w = neff_w
    row2_h = neff_h

    # Combined figure dimensions
    total_w = max(row1_w, row2_w)
    total_h = row1_h + row2_h

    # Create combined image with white background
    combined = Image.new("RGB", (total_w, total_h), "white")

    # Paste row 1 (centered if narrower than total width)
    row1_x_offset = (total_w - row1_w) // 2
    combined.paste(cor_img, (row1_x_offset, 0))
    combined.paste(egene_img, (row1_x_offset + cor_w, 0))

    # Paste row 2 (centered)
    row2_x_offset = (total_w - row2_w) // 2
    combined.paste(neff_img, (row2_x_offset, row1_h))

    # Save as PDF
    output_file = f"{output_path}/f3combined_{qtdid}.pdf"
    combined.save(output_file, "PDF", resolution=300)
    print(f"Saved: {output_file}")


if __name__ == "__main__":
    for qtdid in target_qtdids:
        combine_pdfs(qtdid)
        print(f"Processed {qtdid}")
