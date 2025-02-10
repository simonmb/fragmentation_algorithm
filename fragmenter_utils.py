import io
import random
from rdkit.Chem import AllChem
from rdkit.Chem.Draw import rdMolDraw2D
from PIL import Image, ImageDraw, ImageFont

def get_text_size(draw, text, font):
    """Gets the size of the text when drawn using the given font.

    Args:
        draw (PIL.ImageDraw.ImageDraw): The drawing context.
        text (str): The text whose size needs to be measured.
        font (PIL.ImageFont.ImageFont): The font used for the text.

    Returns:
        tuple: A tuple containing the width and height of the text.
    """
    try:
        return draw.textsize(text, font=font)
    except AttributeError:
        bbox = draw.textbbox((0, 0), text, font=font)
        return (bbox[2] - bbox[0], bbox[3] - bbox[1])

def draw_mol_with_highlights_and_legend(mol, highlight_dict, group_names=None, group_colors=None,
                                        img_size=(600, 600),
                                        legend_width=None,
                                        circle_radius=None,
                                        line_width=None,
                                        highlight_alpha=120):
    """Draws a molecule with highlighted atom groups and a legend.

    Args:
        mol (rdkit.Chem.rdchem.Mol): The RDKit molecule object to be drawn.
        highlight_dict (dict): A dictionary mapping group names to lists of atom indices.
        group_names (dict, optional): A mapping of group keys to custom names. Defaults to None.
        group_colors (dict, optional): A mapping of group keys to RGB color tuples. Defaults to None.
        img_size (tuple, optional): The size of the output image in pixels. Defaults to (600, 600).
        legend_width (int, optional): The width of the legend panel. Defaults to one-third of img_size[0].
        circle_radius (int, optional): The radius of the circles highlighting atoms. Defaults to a size-based calculation.
        line_width (int, optional): The width of the highlight lines. Defaults to a size-based calculation.
        highlight_alpha (int, optional): The transparency level of the highlight overlays (0-255). Defaults to 120.

    Returns:
        PIL.Image.Image: The final image with the molecule drawing and legend.
    """

    if not legend_width:
        legend_width = int(img_size[0]/3)
    if not circle_radius:
        circle_radius = int((img_size[0] / 300) * 15)
    if not line_width:
        line_width = int((img_size[0] / 300) * 12)

    AllChem.Compute2DCoords(mol)
    drawer = rdMolDraw2D.MolDraw2DCairo(img_size[0], img_size[1])
    drawer.drawOptions().addAtomIndices = True
    drawer.drawOptions().bondLineWidth = (img_size[0] / 300) * 2
    drawer.DrawMolecule(mol)
    drawer.FinishDrawing()
    png_data = drawer.GetDrawingText()
    base_img = Image.open(io.BytesIO(png_data)).convert("RGBA")
    
    atom_coords = {}
    for i in range(mol.GetNumAtoms()):
        pt = drawer.GetDrawCoords(i)
        atom_coords[i] = pt
    
    if group_colors is None:
        palette = [
            (31, 119, 180),  # blue
            (255, 127, 14),  # orange
            (44, 160, 44),   # green
            (214, 39, 40),   # red
            (148, 103, 189), # purple
            (140, 86, 75),   # brown
            (227, 119, 194), # pink
            (127, 127, 127), # grey
            (188, 189, 34),  # olive
            (23, 190, 207)   # cyan
        ]
        group_color_mapping = {}
        for i, key in enumerate(highlight_dict.keys()):
            group_color_mapping[key] = palette[i % len(palette)]
    else:
        group_color_mapping = group_colors

    for key, groups in highlight_dict.items():
        overlay = Image.new("RGBA", base_img.size, (255, 255, 255, 0))
        overlay_draw = ImageDraw.Draw(overlay, 'RGBA')
        color = group_color_mapping.get(key, (random.randint(0,255),
                                               random.randint(0,255),
                                               random.randint(0,255)))
        rgba_color = (color[0], color[1], color[2], highlight_alpha)
        
        for group in groups:
            # bonds
            for bond in mol.GetBonds():
                a1 = bond.GetBeginAtomIdx()
                a2 = bond.GetEndAtomIdx()
                if a1 in group and a2 in group:
                    p1 = tuple(atom_coords.get(a1))
                    p2 = tuple(atom_coords.get(a2))
                    if p1 and p2:
                        overlay_draw.line([p1, p2], fill=rgba_color, width=line_width)
            # atoms
            for atom_idx in group:
                if atom_idx in atom_coords:
                    x, y = atom_coords[atom_idx]
                    bbox = [x - circle_radius, y - circle_radius,
                            x + circle_radius, y + circle_radius]
                    overlay_draw.ellipse(bbox, fill=rgba_color)
        
        base_img = Image.alpha_composite(base_img, overlay)
    
    border_width = 5
    bordered_img = Image.new("RGBA",
                             (base_img.width + 2 * border_width, base_img.height + 2 * border_width),
                             (240, 240, 240, 255))
    bordered_img.paste(base_img, (border_width, border_width), base_img)
    
    final_width = bordered_img.width + legend_width
    final_img = Image.new("RGBA", (final_width, bordered_img.height), (255, 255, 255, 255))
    final_img.paste(bordered_img, (0, 0))
    
    legend_area = (bordered_img.width, 0, final_width, bordered_img.height)
    legend_bg_color = (245, 245, 245, 255)
    legend_draw = ImageDraw.Draw(final_img)
    legend_draw.rectangle(legend_area, fill=legend_bg_color)

    margin = 10
    try:
        title_font = ImageFont.truetype("arial.ttf", 16)
    except IOError:
        title_font = ImageFont.load_default()
    title_text = "Groups"
    title_width, title_height = get_text_size(legend_draw, title_text, font=title_font)
    title_x = bordered_img.width + (legend_width - title_width) // 2
    title_y = margin
    legend_draw.text((title_x, title_y), title_text, fill=(0, 0, 0, 255), font=title_font)
    
    current_y = title_y + title_height + margin
    box_size = 20
    spacing = 10
    try:
        item_font = ImageFont.truetype("arial.ttf", 14)
    except IOError:
        item_font = ImageFont.load_default()
    
    for key in group_color_mapping.keys():
        rect_x0 = bordered_img.width + margin
        rect_y0 = current_y
        rect_x1 = rect_x0 + box_size
        rect_y1 = rect_y0 + box_size

        try:
            legend_draw.rounded_rectangle([rect_x0, rect_y0, rect_x1, rect_y1],
                                          radius=4,
                                          fill=(group_color_mapping[key][0],
                                                group_color_mapping[key][1],
                                                group_color_mapping[key][2],
                                                255))
        except AttributeError:
            legend_draw.rectangle([rect_x0, rect_y0, rect_x1, rect_y1],
                                  fill=(group_color_mapping[key][0],
                                        group_color_mapping[key][1],
                                        group_color_mapping[key][2],
                                        255))
        group_name = key
        if group_names:
            group_name = group_names[key]
        text = f"{group_name} ({len(highlight_dict[key])})"
        text_width, text_height = get_text_size(legend_draw, text, font=item_font)
        text_x = rect_x1 + spacing
        text_y = rect_y0 + (box_size - text_height) // 2  # center the text vertically
        legend_draw.text((text_x, text_y), text, fill=(0, 0, 0, 255), font=item_font)
        current_y += box_size + spacing
    
    return final_img
