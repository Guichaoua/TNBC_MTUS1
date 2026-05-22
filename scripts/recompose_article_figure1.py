from pathlib import Path

from PIL import Image, ImageDraw, ImageFont


PROJECT = Path(__file__).resolve().parents[1]
SOURCE = PROJECT / "figures" / "source" / "fig1__v3.png"
OUT_DIR = PROJECT / "figures"
OUT_DIR.mkdir(parents=True, exist_ok=True)


def load_font(size: int, bold: bool = False) -> ImageFont.FreeTypeFont:
    candidates = [
        "/System/Library/Fonts/Supplemental/Arial Bold.ttf" if bold else "/System/Library/Fonts/Supplemental/Arial.ttf",
        "/System/Library/Fonts/Supplemental/Helvetica Bold.ttf" if bold else "/System/Library/Fonts/Supplemental/Helvetica.ttf",
        "/Library/Fonts/Arial Bold.ttf" if bold else "/Library/Fonts/Arial.ttf",
    ]
    for candidate in candidates:
        path = Path(candidate)
        if path.exists():
            return ImageFont.truetype(str(path), size)
    return ImageFont.load_default()


def scaled_rect(rect: tuple[int, int, int, int], scale: int) -> tuple[int, int, int, int]:
    return tuple(v * scale for v in rect)


def paste_scaled_crop(
    canvas: Image.Image,
    source: Image.Image,
    crop_box: tuple[int, int, int, int],
    x: int,
    y: int,
    target_width: int,
    scale: int,
) -> None:
    panel = source.crop(crop_box)
    target = (target_width * scale, round(panel.height * target_width / panel.width) * scale)
    panel = panel.resize(target, Image.Resampling.LANCZOS).convert("RGB")
    canvas.paste(panel, (x * scale, y * scale))


def draw_label(
    draw: ImageDraw.ImageDraw,
    label: str,
    title: str,
    x: int,
    y: int,
    scale: int,
    title_size: int = 34,
) -> None:
    label_font = load_font(42 * scale, bold=True)
    title_font = load_font(title_size * scale, bold=True)
    x *= scale
    y *= scale
    draw.text((x, y), f"{label}.", fill=(28, 28, 28), font=label_font)
    draw.text((x + 64 * scale, y + 5 * scale), title, fill=(28, 28, 28), font=title_font)


def draw_right_aligned_text(
    draw: ImageDraw.ImageDraw,
    text: str,
    right_x: int,
    center_y: int,
    scale: int,
    size: int = 15,
) -> None:
    font = load_font(size * scale)
    bbox = draw.textbbox((0, 0), text, font=font)
    width = bbox[2] - bbox[0]
    height = bbox[3] - bbox[1]
    draw.text(
        ((right_x * scale) - width, (center_y * scale) - height // 2),
        text,
        fill=(28, 28, 28),
        font=font,
    )


def main() -> None:
    src = Image.open(SOURCE).convert("RGBA")
    scale = 3
    canvas = src.resize((src.width * scale, src.height * scale), Image.Resampling.LANCZOS).convert("RGB")
    draw = ImageDraw.Draw(canvas)

    # White out the informal figure title and the old oversized panel headers,
    # while preserving the four original scientific panels and their positions.
    for rect in [
        (0, 0, 210, 62),
        (220, 135, 705, 205),
        (1090, 0, 1534, 58),
        (1110, 395, 1534, 458),
        (700, 760, 1534, 1148),
    ]:
        draw.rectangle(scaled_rect(rect, scale), fill="white")

    draw_label(draw, "A", "TNBC discovery", 225, 146, scale, title_size=36)
    draw_label(draw, "B", "TNBC selected", 1105, 8, scale, title_size=32)
    draw_label(draw, "C", "Breast cancer selected", 1110, 405, scale, title_size=31)
    draw_label(draw, "D", "TCGA selected", 1105, 780, scale, title_size=31)
    paste_scaled_crop(canvas, src, (740, 818, 1520, 1145), 850, 850, 660, scale)
    draw.rectangle(scaled_rect((850, 852, 1034, 1000), scale), fill="white")
    for text, center_y in [
        ("MYC_TARGETS_V1", 866),
        ("MYC_TARGETS_V2", 896),
        ("DNA_REPAIR", 926),
        ("OXIDATIVE_PHOSPHORYLATION", 956),
        ("UNFOLDED_PROTEIN_RESPONSE", 986),
    ]:
        draw_right_aligned_text(draw, text, 1030, center_y, scale)

    png_path = OUT_DIR / "Figure1_Pathway_programme_4panel.png"
    pdf_path = OUT_DIR / "Figure1_Pathway_programme_4panel.pdf"
    canvas.save(png_path, dpi=(600, 600))
    canvas.save(pdf_path, "PDF", resolution=600)
    print(png_path)
    print(pdf_path)


if __name__ == "__main__":
    main()
