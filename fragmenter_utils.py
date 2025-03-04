import io
import random
import math

from PIL import Image, ImageDraw, ImageFont

from rdkit import Chem
from rdkit.Chem.rdchem import ChiralType
from rdkit.Chem import AllChem
from rdkit.Chem.Draw import rdMolDraw2D
from rdkit.Chem.rdchem import ChiralType


def get_text_size(draw, text, font):
    """Calculates the size of the text when rendered using the specified font.

    Args:
        draw (PIL.ImageDraw.ImageDraw): The drawing context used to measure the text.
        text (str): The text string to be measured.
        font (PIL.ImageFont.ImageFont): The font used for rendering the text.

    Returns:
        tuple: A tuple (width, height) representing the dimensions of the text.
    """
    try:
        return draw.textsize(text, font=font)
    except AttributeError:
        bbox = draw.textbbox((0, 0), text, font=font)
        return (bbox[2] - bbox[0], bbox[3] - bbox[1])


def draw_mol_with_highlights_and_legend(mol,
                                        group_dict,
                                        group_names=None,
                                        group_colors=None,
                                        show_aromatic_info=True,
                                        img_width_and_height=None):
    """Draws a 2D depiction of a molecule with highlighted groups and a legend.

    The function renders the molecule using RDKit, highlights specified groups
    of atoms and bonds with semi-transparent colors, and adds a legend showing the
    group names and counts. Optionally, aromatic bonds and atoms can be highlighted.

    Args:
        mol (rdkit.Chem.Mol): The molecule to be drawn. It must consist of a single fragment.
        group_dict (dict): A dictionary where each key maps to a list of groups (each group is
            a list of atom indices) to be highlighted.
        group_names (dict, optional): A mapping of group keys to their display names for the legend.
            Defaults to None.
        group_colors (dict, optional): A mapping of group keys to RGB color tuples for highlighting.
            If None, a default palette is used. Defaults to None.
        show_aromatic_info (bool, optional): If True, aromatic bonds and atoms are highlighted.
            Defaults to True.
        img_width_and_height (int, optional): The width and height (in pixels) for the image.
            If not provided, a suitable size is determined automatically. Defaults to None.

    Returns:
        PIL.Image.Image: An image of the molecule with highlights and a legend.

    Raises:
        ValueError: If the molecule consists of more than one fragment.
    """
    if len(Chem.GetMolFrags(mol)) != 1:
        raise ValueError('This function can only handle molecules with one fragment.')

    def get_bond_length_in_pixel(mol, img_width_and_height):
        drawer = rdMolDraw2D.MolDraw2DCairo(img_width_and_height, img_width_and_height)
        drawer.DrawMolecule(mol)
        drawer.FinishDrawing()
        bonds = list(mol.GetBonds())
        if not bonds:
            return 80
        else:
            bond = bonds[0]
            p1 = drawer.GetDrawCoords(bond.GetBeginAtomIdx())
            p2 = drawer.GetDrawCoords(bond.GetEndAtomIdx())
            dx = p1.x - p2.x
            dy = p1.y - p2.y
            bond_length_in_pixel = math.sqrt(dx * dx + dy * dy)
            return bond_length_in_pixel

    highlight_alpha = 140

    # 1. Draw the molecule normally
    AllChem.Compute2DCoords(mol)

    img_width_and_height_start = img_width_and_height
    if img_width_and_height is None:
        img_width_and_height = 600
        img_width_and_height_start = img_width_and_height

        bond_length_in_pixel = get_bond_length_in_pixel(mol, img_width_and_height)
        pixel_threshold = 70
        while bond_length_in_pixel < pixel_threshold:
            img_width_and_height = int(1.2 * img_width_and_height)
            bond_length_in_pixel = get_bond_length_in_pixel(mol, img_width_and_height)
    else:
        bond_length_in_pixel = get_bond_length_in_pixel(mol, img_width_and_height)

    size_multiplier = max(img_width_and_height / img_width_and_height_start, 1)

    drawer = rdMolDraw2D.MolDraw2DCairo(img_width_and_height, img_width_and_height)
    drawer.drawOptions().addAtomIndices = True
    drawer.drawOptions().bondLineWidth = size_multiplier * 3
    drawer.DrawMolecule(mol)
    drawer.FinishDrawing()
    png_data = drawer.GetDrawingText()
    base_img = Image.open(io.BytesIO(png_data)).convert("RGBA")

    atom_coords = {}
    for i in range(mol.GetNumAtoms()):
        pt = drawer.GetDrawCoords(i)
        atom_coords[i] = pt

    legend_width = int(img_width_and_height / 3)
    circle_radius = int(bond_length_in_pixel / 4)
    line_width = int(bond_length_in_pixel / 5)

    # 2. Draw overlay for groups and aromatic bonds  
    has_aromatic_info = False
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
        for i, key in enumerate(group_dict.keys()):
            group_color_mapping[key] = palette[i % len(palette)]
    else:
        group_color_mapping = group_colors

    overlay = Image.new("RGBA", base_img.size, (255, 255, 255, 0))
    overlay_draw = ImageDraw.Draw(overlay, 'RGBA')
    for key, groups in group_dict.items():
        color = group_color_mapping.get(key, (random.randint(0, 255),
                                               random.randint(0, 255),
                                               random.randint(0, 255)))
        rgba_color = (color[0], color[1], color[2], highlight_alpha)
        
        for group in groups:
            # Highlight bonds within the group
            for bond in mol.GetBonds():
                a1 = bond.GetBeginAtomIdx()
                a2 = bond.GetEndAtomIdx()
                if a1 in group and a2 in group:
                    p1 = tuple(atom_coords.get(a1))
                    p2 = tuple(atom_coords.get(a2))
                    if p1 and p2:
                        overlay_draw.line([p1, p2], fill=rgba_color, width=line_width)
            # Highlight atoms in the group
            for atom_idx in group:
                if atom_idx in atom_coords:
                    x, y = atom_coords[atom_idx]
                    bbox = [x - circle_radius, y - circle_radius,
                            x + circle_radius, y + circle_radius]
                    overlay_draw.ellipse(bbox, fill=rgba_color)
    
    line_color_aromatic = (98, 190, 235, 255)
    line_width_aromatic = int(max((line_width * size_multiplier) / 5, 4))
    if show_aromatic_info:
        for bond in mol.GetBonds():
            if bond.GetBondType() == Chem.rdchem.BondType.AROMATIC:
                has_aromatic_info = True
                p1 = tuple(atom_coords.get(bond.GetBeginAtomIdx()))
                p2 = tuple(atom_coords.get(bond.GetEndAtomIdx()))
                if p1 and p2:
                    overlay_draw.line([p1, p2], fill=line_color_aromatic, width=line_width_aromatic)
        for atom in mol.GetAtoms():
            if atom.GetIsAromatic():
                has_aromatic_info = True
                x, y = drawer.GetDrawCoords(atom.GetIdx())
                bbox = [x - circle_radius, y - circle_radius,
                        x + circle_radius, y + circle_radius]
                overlay_draw.ellipse(bbox, outline=line_color_aromatic, width=line_width_aromatic)

    base_img = Image.alpha_composite(base_img, overlay)

    # 3. Draw legend
    border_width = 5
    bordered_img = Image.new("RGBA",
                             (base_img.width + 2 * border_width, base_img.height + 2 * border_width),
                             (240, 240, 240, 255))
    bordered_img.paste(base_img, (border_width, border_width), base_img)
    
    final_width = bordered_img.width + legend_width
    final_img = Image.new("RGBA", (final_width, bordered_img.height), (255, 255, 255, 255))
    final_img.paste(bordered_img, (0, 0))
    
    legend_area = (bordered_img.width, 0, final_width, bordered_img.height)
    legend_bg_gray_level = 245
    legend_bg_color = (legend_bg_gray_level, legend_bg_gray_level, legend_bg_gray_level, 255)
    legend_draw = ImageDraw.Draw(final_img)
    legend_draw.rectangle(legend_area, fill=legend_bg_color)

    margin = 10 * size_multiplier
    try:
        font = ImageFont.truetype("arial.ttf", int(16 * size_multiplier))
    except IOError:
        font = ImageFont.load_default()
    title_text = "Description"
    title_width, title_height = get_text_size(legend_draw, title_text, font=font)
    title_x = int(bordered_img.width + (legend_width - title_width) / 2)
    title_y = margin
    legend_draw.text((title_x, title_y), title_text, fill=(0, 0, 0, 255), font=font)
    
    current_y = title_y + title_height + margin
    box_size = 20 * size_multiplier
    spacing = 10 * size_multiplier
    
    for key in group_color_mapping:
        rect_x0 = bordered_img.width + margin
        rect_y0 = current_y
        rect_x1 = rect_x0 + box_size
        rect_y1 = rect_y0 + box_size

        original_color = list(group_color_mapping[key])
        original_color.append(highlight_alpha)
        original_color = tuple(original_color)

        # Mix with white
        r, b, g, a = original_color
        ratio_actual = a / 255
        actual_color = [((ratio_actual * v) + (1 - ratio_actual) * 255) for v in original_color[:3]]
        actual_color = tuple([int(min(max(v, 0), 255)) for v in actual_color])

        try:
            legend_draw.rounded_rectangle([rect_x0, rect_y0, rect_x1, rect_y1],
                                          radius=4,
                                          fill=actual_color)
        except AttributeError:
            legend_draw.rectangle([rect_x0, rect_y0, rect_x1, rect_y1],
                                  fill=actual_color)
        group_name = key
        if group_names:
            group_name = group_names[key]
        text = f"{group_name} ({len(group_dict[key])})"
        _, text_height = get_text_size(legend_draw, text, font=font)
        text_x = rect_x1 + spacing
        text_y = int(rect_y0 + (box_size - text_height) / 2)  # Center the text vertically
        legend_draw.text((text_x, text_y), text, fill=(0, 0, 0, 255), font=font)
        current_y += box_size + spacing

    if has_aromatic_info:
        line_x0 = bordered_img.width + margin
        line_x1 = line_x0 + box_size
        text = "arom. atom/bond"
        _, text_height = get_text_size(legend_draw, text, font=font)
        text_x = line_x1 + spacing
        text_y = int(current_y + (box_size - text_height) / 2)  # Center vertically
        line_y = int(text_y + box_size / 2)
        legend_draw.line([(line_x0, line_y), (line_x1, line_y)], fill=line_color_aromatic, width=line_width_aromatic)
        legend_draw.text((text_x, text_y), text, fill=(0, 0, 0, 255), font=font)

    return final_img


def get_table_with_atom_properties_relevant_to_SMARTS(mol):
    """Generates a table of atom properties relevant to SMARTS pattern matching.

    The function extracts various properties for each atom in the molecule such as
    atomic symbol, charge, aromaticity, degree, hydrogen counts, ring information,
    and chirality. These properties can assist in understanding how SMARTS patterns
    interact with molecular structures.

    Args:
        mol (rdkit.Chem.Mol): The molecule for which the atom properties are to be tabulated.

    Returns:
        tuple: A tuple containing four elements:
            - headers1 (list): Primary header labels.
            - headers2 (list): Secondary header labels providing additional context.
            - data (list of lists): A list of rows with atom properties.
            - formatted_rows (list): A list of strings, each representing a formatted row.
    """
    # Reference: https://www.daylight.com/dayhtml/doc/theory/theory.smarts.html

    headers1 = [
        "idx", "Sym[AN]", "Charge", "Arom", "Degree",
        "Tot. Hs", "Impl. Hs", "R. Count",
        "R. Size", "Val.", "Conn.", "R. Conn.",
        "Chir.", "CIP"
    ]
    headers2 = [
        "", "", "", "/Aliph", "(D<n>)",
        "(H<n>)", "(h<n>)", "(R<n>)",
        "(r<n>)", "(v<n>)", "(X<n>)", "(x<n>)",
        "", ""
    ]
    
    data = []
    ri = mol.GetRingInfo()
    for atom in mol.GetAtoms():
        has_chirality = atom.GetChiralTag() != ChiralType.CHI_UNSPECIFIED
        in_rings_of_size = []
        for n in range(3, 11):
            if atom.IsInRingSize(n):
                in_rings_of_size.append(n)
        if not in_rings_of_size and atom.IsInRing():
            in_rings_of_size = ['>10']

        ring_connectivity = sum(1 for bond in atom.GetBonds() if ri.NumBondRings(bond.GetIdx()) > 0)
        
        idx = atom.GetIdx()
        row = [
            f"{idx}",
            f"{atom.GetSymbol()} [{atom.GetAtomicNum()}]",
            str(atom.GetFormalCharge()),
            "a" if atom.GetIsAromatic() else "A",
            f"{atom.GetDegree()}",
            f"{atom.GetTotalNumHs()}",
            f"{atom.GetImplicitValence()}",
            str(ri.NumAtomRings(idx)),
            str(ri.AtomRingSizes(idx)) if atom.IsInRing() else "",
            f"{atom.GetTotalValence()}",
            f"{atom.GetTotalDegree()}",
            str(ring_connectivity),
            "✔" if has_chirality else "✘",
            f"{atom.GetProp('_CIPCode') if atom.HasProp('_CIPCode') else ''}"
        ]
        data.append(row)
    
    col_widths = []
    for i in range(len(headers1)):
        max_width = max(len(headers1[i]), len(headers2[i]))
        for row in data:
            if len(row[i]) > max_width:
                max_width = len(row[i])
        col_widths.append(max_width)
    
    fmt = " | ".join(f"{{:<{w}}}" for w in col_widths)
    sep = "-+-".join("-" * w for w in col_widths)
    
    formatted_rows = []
    formatted_rows.append(fmt.format(*headers1))
    formatted_rows.append(fmt.format(*headers2))
    formatted_rows.append(sep)
    for row in data:
        formatted_rows.append(fmt.format(*row))

    return headers1, headers2, data, formatted_rows
