from typing import Dict, Tuple
from architecture import model
from attr import dataclass
import torch
from api.model.dto import (
    Atom2D, Bond2D, PredictRequest, PredictResponseData, PredictSingleRequest, PredictSingleResponseData,
    PredictMultipleRequest, PredictMultipleResponseData, FragmentRequest, FragmentResponseData,
    PredictCheckRequest, PredictCheckResponseData, InferAllRequest, InferAllResponseData,
    DownloadReportRequest, DownloadReportResponseData, PredictedBond, EvaluatedFragmentBond
)
import logging
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.Draw import rdMolDraw2D
import hashlib
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)
from deepbde.architecture.inference_util import multi_predict, single_predict
from io import BytesIO
from reportlab.lib.pagesizes import letter
from reportlab.pdfgen import canvas
import cairosvg

@dataclass
class MoleculeInfo:
    """
    Información completa de una molécula generada a partir de SMILES.
    - molecule_id: ID único (hash SHA256)
    - smiles_canonical: SMILES canónico con H explícitos
    - canvas: Tamaño del lienzo para visualización
    - atoms: Diccionario con información 2D de átomos
    - bonds: Diccionario con información 2D de enlaces
    - image_svg: Representación SVG enriquecida de la molécula
    """
    mol: Chem.Mol
    molecule_id: str
    smiles_canonical: str
    canvas: Dict[str, int]
    atoms: Dict[str, Atom2D]
    bonds: Dict[str, Bond2D]
    image_svg: str

def generate_id_smiles(smiles: str) -> Tuple[Chem.Mol, str, str]:
    """Generates a unique ID and canonical SMILES representation for a molecule.
    Args:
        smiles (str): The SMILES representation of the molecule.
    Returns:
        Tuple[Chem.Mol, str, str]: Molecule object, unique ID, and canonical SMILES.
    """
    mol = None
    try:
        mol = Chem.MolFromSmiles(smiles)
    except Exception as e:
        raise ValueError(f"Invalid SMILES: {smiles}") from e
    if mol is None:
        raise ValueError(f"Invalid SMILES: {smiles}")
    Chem.SanitizeMol(mol)
    mol = Chem.AddHs(mol)
    smiles_canonical = Chem.MolToSmiles(mol, canonical=True, allHsExplicit=True, kekuleSmiles=True, isomericSmiles=True)
    mol = Chem.MolFromSmiles(smiles_canonical)  
    Chem.SanitizeMol(mol)
    mol = Chem.AddHs(mol)
    smiles_canonical = Chem.MolToSmiles(mol, canonical=True, allHsExplicit=True, kekuleSmiles=True, isomericSmiles=True)
    mol_id = hashlib.sha256(smiles_canonical.encode()).hexdigest()[:16]
    return mol, mol_id, smiles_canonical

def verify_smiles(smiles: str, mol_id: str) -> bool:
    """Verifies if the provided molecule ID matches the generated ID from SMILES.
    
    Args:
        smiles (str): The SMILES representation of the molecule.
        mol_id (str): Unique ID to verify.
    Returns:
        bool: True if IDs match, raises ValueError otherwise.
    """
    real_mol_id = generate_id_smiles(smiles)[1]
    if mol_id != real_mol_id:
        raise ValueError(f"SMILES ID mismatch: expected {real_mol_id}, got {mol_id}")
    return True

def is_valid_for_bde(mol: MoleculeInfo, bond_idx: int) -> bool:
    """Checks if a bond at the given index is a single bond.
    
    Args:
        mol (MoleculeInfo): The molecule information object.
        bond_idx (int): The index of the bond to check.
    Returns:
        bool: True if the bond is single, False otherwise.
    """
    if bond_idx < 0 or bond_idx >= mol.GetNumBonds(): # pyright: ignore[reportAttributeAccessIssue]
        raise ValueError(f"Bond index {bond_idx} out of range for molecule with {mol.GetNumBonds()} bonds") # pyright: ignore[reportAttributeAccessIssue]
    bond = mol.GetBondWithIdx(bond_idx) # type: ignore
    return bond.GetBondType() ==  (Chem.BondType.SINGLE and not bond.IsInRing()) # pyright: ignore[reportAttributeAccessIssue]

def calculate_canvas_size(mol: Chem.Mol, padding: int = 50, min_size: int = 300) -> Dict[str, int]:
    """Calculates the canvas size based on the molecule's bounding box.
    
    Args:
        mol (Chem.Mol): The RDKit molecule object with 2D coordinates.
        padding (int): Padding to add around the molecule.
        min_size (int): Minimum width and height of the canvas.
    Returns:
        Dict[str, int]: Dictionary with 'width' and 'height' keys.
    """
    conf = mol.GetConformer()
    positions = [conf.GetAtomPosition(i) for i in range(mol.GetNumAtoms())]# type: ignore[attr-defined]
    min_x, max_x = min(pos.x for pos in positions), max(pos.x for pos in positions)# type: ignore[attr-defined]
    min_y, max_y = min(pos.y for pos in positions), max(pos.y for pos in positions)# type: ignore[attr-defined]
    width = int(max_x - min_x) + padding# type: ignore[attr-defined]
    height = int(max_y - min_y) + padding# type: ignore[attr-defined]
    return {"width": max(width, min_size), "height": max(height, min_size)}

def generate_molecule_svg(mol: Chem.Mol, canvas: Dict[str, int]) -> Tuple[str, rdMolDraw2D.MolDraw2DSVG]:
    """Generates an SVG representation of the molecule with atom and bond indices in the default color.
    Args:
        mol (Chem.Mol): The RDKit molecule object with 2D coordinates.
        canvas (Dict[str, int]): Dictionary with 'width' and 'height' of the canvas.
    Returns:
        Tuple[str, rdMolDraw2D.MolDraw2DSVG]: SVG string of the molecule drawing and the drawer object.
    """
    drawer = rdMolDraw2D.MolDraw2DSVG(canvas["width"], canvas["height"])
    opts = drawer.drawOptions()
    opts.addStereoAnnotation = True
    opts.explicitMethyl = True
    opts.includeAtomTags = True
    opts.addAtomIndices = True  # Enable atom indices in the default color (typically black)
    opts.addBondIndices = True  # Enable bond indices in the default color (typically black)
    drawer.DrawMolecule(mol)# type: ignore[attr-defined]
    drawer.FinishDrawing()
    raw_svg = drawer.GetDrawingText()
    return ''.join(line.strip() for line in raw_svg.splitlines()), drawer
ignore_missing_imports = True
def get_atoms_info(mol: Chem.Mol, drawer: rdMolDraw2D.MolDraw2DSVG) -> Dict[str, Atom2D]:
    """Extracts atom information from the molecule.
    
    Args:
        mol (Chem.Mol): The RDKit molecule object.
        drawer (rdMolDraw2D.MolDraw2DSVG): The RDKit SVG drawer object.
    Returns:
        Dict[str, Atom2D]: Dictionary mapping atom indices to Atom2D objects.
    """
    atoms = {}
    for atom in mol.GetAtoms():
        idx = atom.GetIdx()
        pt = drawer.GetDrawCoords(idx)
        symbol = atom.GetSymbol()
        atom_smiles = Chem.MolToSmiles(Chem.MolFromSmiles(f'[{symbol}]')) if symbol else None
        atoms[str(idx)] = Atom2D(
            x=float(pt.x),
            y=float(pt.y),
            symbol=symbol,
            smiles=atom_smiles
        )
    return atoms

def get_bonds_info(mol: Chem.Mol, drawer: rdMolDraw2D.MolDraw2DSVG) -> Dict[str, Bond2D]:
    """Extracts bond information from the molecule.
    Args:
        mol (Chem.Mol): The RDKit molecule object.
        drawer (rdMolDraw2D.MolDraw2DSVG): The RDKit SVG drawer object.
    Returns:
        Dict[str, Bond2D]: Dictionary mapping bond indices to Bond2D objects.
    """
    bonds = {}
    for bond in mol.GetBonds():
        b_idx = bond.GetIdx()
        start_idx = bond.GetBeginAtomIdx()
        end_idx = bond.GetEndAtomIdx()
        start_pt = drawer.GetDrawCoords(start_idx)
        end_pt = drawer.GetDrawCoords(end_idx)
        atom1 = mol.GetAtomWithIdx(start_idx)
        atom2 = mol.GetAtomWithIdx(end_idx)
        bond_type = bond.GetBondType()
        bond_type_map = {
            Chem.BondType.SINGLE: "single",
            Chem.BondType.DOUBLE: "double",
            Chem.BondType.TRIPLE: "triple",
            Chem.BondType.AROMATIC: "aromatic"
        }
        bond_type_str = bond_type_map.get(bond_type, str(bond_type).lower())
        bond_atoms = f"{atom1.GetSymbol()}-{atom2.GetSymbol()}"
        bonds[str(b_idx)] = Bond2D(
            start=start_idx,
            end=end_idx,
            start_coords={"x": float(start_pt.x), "y": float(start_pt.y)},
            end_coords={"x": float(end_pt.x), "y": float(end_pt.y)},
            bond_atoms=bond_atoms,
            bond_type=bond_type_str if bond_type_str in {"single", "double", "triple", "aromatic"} else "single"  # type: ignore
        )
    return bonds

def get_all_info_molecule(smiles: str) -> MoleculeInfo:
    """
    Extracts all relevant information from a molecule given its SMILES representation.
    Args:
        smiles (str): The SMILES representation of the molecule.
    Returns:
        MoleculeInfo: An object containing all relevant information about the molecule.
    """
    mol, mol_id, smiles_canonical = generate_id_smiles(smiles)
    AllChem.Compute2DCoords(mol) # pyright: ignore[reportAttributeAccessIssue]
    canvas = calculate_canvas_size(mol)
    image_svg, drawer = generate_molecule_svg(mol, canvas)
    atoms = get_atoms_info(mol, drawer)
    bonds = get_bonds_info(mol, drawer)
    assert len(atoms) == mol.GetNumAtoms(), f"Atom count mismatch: expected {mol.GetNumAtoms()}, got {len(atoms)}"
    assert len(bonds) == mol.GetNumBonds(), f"Bond count mismatch: expected {mol.GetNumBonds()}, got {len(bonds)}"
    return MoleculeInfo(
        mol=mol,
        molecule_id=mol_id,
        smiles_canonical=smiles_canonical,
        canvas=canvas,
        atoms=atoms,
        bonds=bonds,
        image_svg=image_svg
    )
    
def predict_controller(request: PredictRequest) -> PredictResponseData:
    """Controller for predicting the bond dissociation energy (BDE) of a molecule.
    Args:
        request (PredictRequest): Request object containing the SMILES string.
    Returns:
        PredictResponseData: Response containing molecule data and SVG.
    """
    all_info = get_all_info_molecule(request.smiles)
    return PredictResponseData(
        smiles_canonical=all_info.smiles_canonical,
        image_svg=all_info.image_svg,
        canvas=all_info.canvas,
        atoms=all_info.atoms,
        bonds=all_info.bonds,
        molecule_id=all_info.molecule_id
    )


def get_bde_for_bond_indices(mol_info: MoleculeInfo, idx: int ) -> float | None:
    """Predicts the bond dissociation energy (BDE) for a specific bond index in a molecule.
    Args:
        smiles (str): The SMILES representation of the molecule.
        idx (int): The index of the bond to predict.
    Returns:
        float: The predicted BDE for the specified bond.
    """
    if not is_valid_for_bde(mol_info.mol, idx):
        return None
    bde = single_predict(mol_info.smiles_canonical, idx)
    return bde.item() if isinstance(bde, torch.Tensor) else bde
    
def predict_single_controller(request: PredictSingleRequest) -> PredictSingleResponseData:
    """
    Controller para la predicción de BDE de un único enlace de una molécula.
    Args:
        request (PredictSingleRequest): Contiene el SMILES, el ID de la molécula y el índice del enlace.
    Returns:
        PredictSingleResponseData: Contiene el SMILES canónico y la predicción del enlace.
    """
    all_info = get_all_info_molecule(request.smiles)
    verify_smiles(all_info.smiles_canonical, request.molecule_id)
    mol = all_info.mol
    if not is_valid_for_bde(mol, request.bond_idx):
        raise ValueError(f"Bond index {request.bond_idx} is not a single bond or is in a ring.")
    bond = mol.GetBondWithIdx(request.bond_idx)
    begin_idx = bond.GetBeginAtomIdx()
    end_idx = bond.GetEndAtomIdx()
    atom1 = mol.GetAtomWithIdx(begin_idx)
    atom2 = mol.GetAtomWithIdx(end_idx)
    bde = get_bde_for_bond_indices(all_info, request.bond_idx)
    predicted_bond = PredictedBond(
        idx=request.bond_idx,
        bde=bde,
        begin_atom_idx=begin_idx,
        end_atom_idx=end_idx,
        bond_atoms=f"{atom1.GetSymbol()}-{atom2.GetSymbol()}",
        bond_type=bond.GetBondType().name.lower()
    )
    return PredictSingleResponseData(
        smiles_canonical=all_info.smiles_canonical,
        bond=predicted_bond
    )

def predict_multiple_controller(request: PredictMultipleRequest) -> PredictMultipleResponseData:
    """
    Controller para la predicción de BDE de múltiples enlaces de una molécula.
    Args:
        request (PredictMultipleRequest):
            - smiles (str): SMILES de la molécula.
            - bond_indices (List[int]): Lista de índices de enlaces a predecir.
    Returns:
        PredictMultipleResponseData: Respuesta con la lista de enlaces y sus BDEs predichos.
    """
    all_info = get_all_info_molecule(request.smiles)
    verify_smiles(all_info.smiles_canonical, request.molecule_id)

    predicted_bonds = []
    for idx in request.bond_indices:
        if not is_valid_for_bde(all_info.mol, idx):
            raise ValueError(f"Bond index {idx} is not a single bond or is in a ring.")

        bde = get_bde_for_bond_indices(all_info, idx)
        bond = all_info.mol.GetBondWithIdx(idx)
        begin_idx = bond.GetBeginAtomIdx()
        end_idx = bond.GetEndAtomIdx()
        atom1 = all_info.mol.GetAtomWithIdx(begin_idx)
        atom2 = all_info.mol.GetAtomWithIdx(end_idx)

        predicted_bond = PredictedBond(
            idx=idx,
            bde=bde,
            begin_atom_idx=begin_idx,
            end_atom_idx=end_idx,
            bond_atoms=f"{atom1.GetSymbol()}-{atom2.GetSymbol()}",
            bond_type=bond.GetBondType().name.lower()
        )
        predicted_bonds.append(predicted_bond)

    return PredictMultipleResponseData(
        smiles=request.smiles,
        molecule_id=request.molecule_id,
        bonds=predicted_bonds
    )
    
def infer_all_controller(request: InferAllRequest) -> InferAllResponseData:
    """
    Controller para inferir la BDE de todos los enlaces posibles de una molécula.
    Args:
        request (InferAllRequest): Contiene el SMILES de la molécula.
    Returns:
        InferAllResponseData: Respuesta con la lista de enlaces y sus BDEs predichos o null si no se pueden predecir.
    """
    all_info = get_all_info_molecule(request.smiles)
    mol = all_info.mol

    predicted_bonds = []
    for bond in mol.GetBonds():
        bond_idx = bond.GetIdx()
        begin_idx = bond.GetBeginAtomIdx()
        end_idx = bond.GetEndAtomIdx()
        atom1 = mol.GetAtomWithIdx(begin_idx)
        atom2 = mol.GetAtomWithIdx(end_idx)
        if is_valid_for_bde(mol, bond_idx):
            bde = get_bde_for_bond_indices(all_info, bond_idx)
        else:
            bde = None
        predicted_bond = PredictedBond(
            idx=bond_idx,
            bde=bde,
            begin_atom_idx=begin_idx,
            end_atom_idx=end_idx,
            bond_atoms=f"{atom1.GetSymbol()}-{atom2.GetSymbol()}",
            bond_type=bond.GetBondType().name.lower()
        )
        predicted_bonds.append(predicted_bond)
    return InferAllResponseData(
        smiles_canonical=all_info.smiles_canonical,
        molecule_id=all_info.molecule_id,
        bonds=predicted_bonds
    )

def fragment_controller(request: FragmentRequest) -> FragmentResponseData:
    # Esta función aún no está implementada. Se debe completar según los requisitos del proyecto.
    return FragmentResponseData(
        smiles_canonical="",
        molecule_id="",
        bonds=[]
    )
    
    
    
def get_fragments_from_bond(mol: Chem.Mol, bond_idx: int) -> list[str]:
    """Gets the two fragments generated by breaking a bond in a molecule.

    Args:
        mol (Chem.Mol): The RDKit molecule object.
        bond_idx (int): The index of the bond to break.

    Returns:
        list[str]: A list of SMILES strings for the two fragments.
    """
    # Create a copy of the molecule to modify
    mol_copy = Chem.RWMol(mol)
    bond = mol_copy.GetBondWithIdx(bond_idx)
    begin_idx = bond.GetBeginAtomIdx()
    end_idx = bond.GetEndAtomIdx()

    # Remove the bond to split the molecule into fragments
    mol_copy.RemoveBond(begin_idx, end_idx)

    # Get the fragments as SMILES
    fragments = Chem.GetMolFrags(mol_copy, asMols=True)
    return [Chem.MolToSmiles(frag, canonical=True) for frag in fragments]

def  report_txt(smile:str) -> str:
    """Generates a report in text format for the molecule.
    
    Args:
        smile (str): The SMILES representation of the molecule.
    Returns:
        str: Text report of the molecule.
    """
    # El txt debe tener el smile introducido el id, todas las bde 
    # y que enlaces se selecciono y que moleculas da como resultado al romper el enlace.
    all_info = get_all_info_molecule(smile)
    report_lines = [
        f"Molecule ID: {all_info.molecule_id}",
        f"Canonical SMILES: {all_info.smiles_canonical}",
        "Bonds and BDEs:"
    ]
    for bond in all_info.mol.GetBonds():
        bond_idx = bond.GetIdx()
        begin_idx = bond.GetBeginAtomIdx()
        end_idx = bond.GetEndAtomIdx()
        atom1 = all_info.mol.GetAtomWithIdx(begin_idx)
        atom2 = all_info.mol.GetAtomWithIdx(end_idx)
        bde = get_bde_for_bond_indices(all_info, bond_idx)
        report_lines.append(
            f"Bond {bond_idx}: {atom1.GetSymbol()}-{atom2.GetSymbol()} (BDE: {bde})"
        )
    report_lines.append("Fragmentation Results:")
    for bond in all_info.mol.GetBonds():
        bond_idx = bond.GetIdx()
        bde = get_bde_for_bond_indices(all_info, bond_idx)
        if bde is not None:
            fragments = get_fragments_from_bond(all_info.mol, bond_idx)
            report_lines.append(
                f"Fragment from Bond {bond_idx}: (BDE: {bde})"
            )
            report_lines.append(f"Generated Fragments: {fragments[0]} and {fragments[1]}")
    report_lines="\n".join(report_lines)
    logger.info("Generated report:\n%s", report_lines)    
    return report_lines


def report_pdf(smile:str) -> str:
    """Generates a report in PDF format for the molecule.
    
    Args:
        smile (str): The SMILES representation of the molecule.
    Returns:
        str: Base64-encoded PDF report of the molecule.
    """
    all_info = get_all_info_molecule(smile)

    # Create a buffer to store the PDF
    buffer = BytesIO()
    pdf_canvas = canvas.Canvas(buffer, pagesize=letter)

    # Convert the SVG to PNG using cairosvg
    svg_data = all_info.image_svg
    png_buffer = BytesIO()
    cairosvg.svg2png(bytestring=svg_data.encode("utf-8"), write_to=png_buffer)
    png_buffer.seek(0)

    # Add the PNG image of the main molecule to the PDF
    pdf_canvas.drawImage(png_buffer, 50, 700, width=200, height=200)  # Position the PNG image

    # Add the molecule data next to the image
    text_x = 300  # Position for the text
    text_y = 800
    pdf_canvas.setFont("Helvetica", 10)
    pdf_canvas.drawString(text_x, text_y, f"Molecule ID: {all_info.molecule_id}")
    text_y -= 15
    pdf_canvas.drawString(text_x, text_y, f"Canonical SMILES: {all_info.smiles_canonical}")
    text_y -= 15
    pdf_canvas.drawString(text_x, text_y, "Bonds and BDEs:")
    text_y -= 15

    for bond in all_info.mol.GetBonds():
        bond_idx = bond.GetIdx()
        begin_idx = bond.GetBeginAtomIdx()
        end_idx = bond.GetEndAtomIdx()
        atom1 = all_info.mol.GetAtomWithIdx(begin_idx)
        atom2 = all_info.mol.GetAtomWithIdx(end_idx)
        bde = get_bde_for_bond_indices(all_info, bond_idx)
        pdf_canvas.drawString(
            text_x, text_y, f"Bond {bond_idx}: {atom1.GetSymbol()}-{atom2.GetSymbol()} (BDE: {bde})"
        )
        text_y -= 15

    pdf_canvas.drawString(50, text_y, "Fragmentation Results:")
    text_y -= 15

    for bond in all_info.mol.GetBonds():
        bond_idx = bond.GetIdx()
        bde = get_bde_for_bond_indices(all_info, bond_idx)
        if bde is not None:
            fragments = get_fragments_from_bond(all_info.mol, bond_idx)

            # Add a horizontal line to separate sections
            pdf_canvas.line(50, text_y, 550, text_y)
            text_y -= 10

            # Add BDE information
            pdf_canvas.drawString(50, text_y, f"Bond {bond_idx} (BDE: {bde}):")
            text_y -= 15

            # Generate and add images for each fragment
            for i, fragment_smiles in enumerate(fragments):
                fragment_mol = Chem.MolFromSmiles(fragment_smiles)
                fragment_svg, _ = generate_molecule_svg(fragment_mol, all_info.canvas)
                fragment_png_buffer = BytesIO()
                cairosvg.svg2png(bytestring=fragment_svg.encode("utf-8"), write_to=fragment_png_buffer)
                fragment_png_buffer.seek(0)

                # Position fragment images below the BDE info
                pdf_canvas.drawImage(fragment_png_buffer, 50 + (i * 250), text_y - 200, width=200, height=200)

            # Adjust text_y for the next section
            text_y -= 220

    # Finalize the PDF
    pdf_canvas.save()
    buffer.seek(0)

    # Return the PDF as a base64-encoded string
    pdf_content = buffer.read()
    buffer.close()
    return pdf_content.decode("latin1")


def download_report_controller(request: DownloadReportRequest) -> DownloadReportResponseData:
    """Generates and downloads a report for the molecule. with all information of bde and fragments.
    Args:
        request (DownloadReportRequest): _description_

    Returns:
        DownloadReportResponseData: _description_
    """
    # Esta función aún no está implementada. Se debe completar según los requisitos del proyecto.
    # Aquí se debe generar el informe y codificarlo en base64 o un txt.
    if request.format not in {"pdf", "txt"}:
        raise ValueError("Formato no soportado: use 'pdf' o 'txt'.")
    if request.format == "txt":
        report_content = report_txt(request.smiles)
    else:
        report_content = report_pdf(request.smiles)
    report_base64 = report_content.encode('utf-8').decode('utf-8')
    return DownloadReportResponseData(
        type="pdf" if request.format == "pdf" else "txt",
        report_base64=report_base64
    )
    
    
    
    
    
    
def predict_check_controller(request: PredictCheckRequest) -> PredictCheckResponseData:
    # Esta función aún no está implementada. Se debe completar según los requisitos del proyecto.
    dummy_bond = PredictedBond(
        idx=0,
        bde=0.0,
        begin_atom_idx=0,
        end_atom_idx=1,
        bond_atoms="H-H",
        bond_type="single"
    )
    return PredictCheckResponseData(
        smiles_canonical="",
        bond=dummy_bond,
        products=[]
    )





