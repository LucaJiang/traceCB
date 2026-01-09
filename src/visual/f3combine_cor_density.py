# Combine correlation density figures for all studies and individual studies
# Row 1: all_studies (2x width) + 2 individual studies
# Row 2: 4 individual studies
# Row 3: 4 individual studies
from utils import *
import pypdfium2 as pdfium
from PIL import Image, ImageDraw, ImageFont

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


def add_title_to_image(img, title, font_size=60):
    """Add title text at the top of an image"""
    # Create a new image with extra space for title
    title_height = int(font_size * 1.5)
    new_img = Image.new("RGB", (img.width, img.height + title_height), "white")

    # Paste original image below title area
    new_img.paste(img, (0, title_height))

    # Draw title text
    draw = ImageDraw.Draw(new_img)
    try:
        # Try to use a nice font
        font = ImageFont.truetype(
            "/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf", font_size
        )
    except:
        # Fall back to default font
        font = ImageFont.load_default()

    # Calculate text position (centered)
    bbox = draw.textbbox((0, 0), title, font=font)
    text_width = bbox[2] - bbox[0]
    text_x = (img.width - text_width) // 2
    text_y = (title_height - font_size) // 2

    draw.text((text_x, text_y), title, fill="black", font=font)

    return new_img


def combine_cor_density_figures():
    """Combine all correlation density figures into one composite figure"""
    # Load all_studies figure
    all_studies_pdf = f"{input_path}/f3cor_density_all_studies.pdf"
    if not os.path.exists(all_studies_pdf):
        print(f"Error: {all_studies_pdf} not found")
        return

    all_studies_img = pdf_to_image(all_studies_pdf)
    # Load individual study figures
    individual_imgs = []
    for qtdid in target_qtdids:
        pdf_path = f"{input_path}/f3cor_density_{qtdid}.pdf"
        if not os.path.exists(pdf_path):
            print(f"Warning: {pdf_path} not found, skipping {qtdid}")
            continue

        img = pdf_to_image(pdf_path)
        individual_imgs.append(img)

    if len(individual_imgs) != 10:
        print(f"Warning: Expected 10 individual images, got {len(individual_imgs)}")
        return

    # Get dimensions
    # all_studies image should be roughly 2x the width of individual images
    single_w, single_h = individual_imgs[0].size
    all_w, all_h = all_studies_img.size

    # Resize all_studies to match 2x individual width if needed
    target_all_w = single_w * 2
    if abs(all_w - target_all_w) > 50:  # Allow some tolerance
        # Resize all_studies image to 2x width while maintaining aspect ratio
        scale = target_all_w / all_w
        new_all_h = int(all_h * scale)
        all_studies_img = all_studies_img.resize(
            (target_all_w, new_all_h), Image.LANCZOS
        )
        all_w, all_h = all_studies_img.size

    # Calculate row heights (use max height in each row)
    row1_h = max(all_h, single_h)
    row2_h = single_h
    row3_h = single_h

    # Total dimensions
    total_w = single_w * 4  # 4 images wide
    total_h = row1_h + row2_h + row3_h

    # Create combined image with white background
    combined = Image.new("RGB", (total_w, total_h), "white")

    # Row 1: all_studies (2x width) + 2 individual images
    y_offset = 0
    x_offset = 0
    # Center all_studies vertically in row 1
    y_center_offset = (row1_h - all_h) // 2
    combined.paste(all_studies_img, (x_offset, y_offset + y_center_offset))
    x_offset += all_w

    for i in range(2):
        y_center_offset = (row1_h - single_h) // 2
        combined.paste(individual_imgs[i], (x_offset, y_offset + y_center_offset))
        x_offset += single_w

    # Row 2: 4 individual images
    y_offset += row1_h
    x_offset = 0
    for i in range(2, 6):
        combined.paste(individual_imgs[i], (x_offset, y_offset))
        x_offset += single_w

    # Row 3: 4 individual images
    y_offset += row2_h
    x_offset = 0
    for i in range(6, 10):
        combined.paste(individual_imgs[i], (x_offset, y_offset))
        x_offset += single_w

    # Save as PDF
    output_file = f"{output_path}/f3cor_density_combined.pdf"
    combined.save(output_file, "PDF", resolution=300)
    print(f"Saved: {output_file}")


if __name__ == "__main__":
    combine_cor_density_figures()
