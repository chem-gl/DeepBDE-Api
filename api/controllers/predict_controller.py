from api.model.dto import (
    PredictRequest, PredictResponseData, PredictSingleRequest, PredictSingleResponseData,
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
from typing import Dict
def generate_id_smiles(smiles: str) -> tuple[Chem.Mol, str, str]:
    """Generates a unique ID and canonical SMILES representation for a molecule.
    Args:
        smiles (str): The SMILES representation of the molecule.
    Returns:
        str, str: A unique ID and the canonical SMILES with explicit hydrogens for the molecule.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"Invalid SMILES: {smiles}")
    Chem.SanitizeMol(mol)
    mol = Chem.AddHs(mol)
    
    smiles_canonical = Chem.MolToSmiles(mol,  canonical=True, allHsExplicit=True,kekuleSmiles =True,isomericSmiles=True)

    mol_id = hashlib.sha256(smiles_canonical.encode()).hexdigest()[:16]
    logger.info(f"Received SMILES: {smiles}")
    logger.info(f"Canonical SMILES: {smiles_canonical}")
    logger.info(f"Generated molecule ID: {mol_id}")
    return  mol, mol_id, smiles_canonical
def verify_smiles(smiles: str, mol_id: str) -> bool:
    """
    Verifies the SMILES and generates a molecule object.
    Args:
        smiles (str): The SMILES representation of the molecule.
        mol_id (str): Unique ID of the molecule.
    """
    real_mol_id = generate_id_smiles(smiles)[1]
    if mol_id != real_mol_id:
        raise ValueError(f"SMILES ID mismatch: expected {real_mol_id}, got {mol_id}")
    return True
def predict_controller(request: PredictRequest) -> PredictResponseData:
    """
    Controller for predicting the bond dissociation energy (BDE) of a molecule.
    Args:
        smiles (str): The SMILES of the molecule.
    Returns:
        PredictResponseData: The response containing the canonical SMILES, SVG image, canvas dimensions, etc.
    Nota: Si se integra con un framework async (ej. FastAPI), considerar hacer esta función async.
    """
    mol, mol_id, smiles_canonical = generate_id_smiles(request.smiles)
    AllChem.Compute2DCoords(mol)
    # --- Lógica para calcular el tamaño dinámico del lienzo ---
    # Calcular el bounding box de la molécula
    min_x, max_x = float('inf'), float('-inf')
    min_y, max_y = float('inf'), float('-inf')
    for atom in mol.GetAtoms():
        pos = mol.GetConformer().GetAtomPosition(atom.GetIdx())
        if pos.x < min_x: min_x = pos.x
        if pos.x > max_x: max_x = pos.x
        if pos.y < min_y: min_y = pos.y
        if pos.y > max_y: max_y = pos.y
    # Calcular dimensiones con un margen de 50 pixels
    padding = 50
    width = int(max_x - min_x) + padding
    height = int(max_y - min_y) + padding
    # Asegurar un tamaño mínimo para moléculas pequeñas
    width = max(width, 300)
    height = max(height, 300)
    drawer = rdMolDraw2D.MolDraw2DSVG(width, height)#6195d7c5aefdd3c0
    opts = drawer.drawOptions()
    opts.addStereoAnnotation = True
    # Forzar la visualización explícita de todos los carbonos y nitrógenos
    opts.explicitMethyl = True
    opts.includeAtomTags = True
    drawer.DrawMolecule(mol)
    drawer.FinishDrawing()
    raw_svg = drawer.GetDrawingText()
    image_svg = ''.join(line.strip() for line in raw_svg.splitlines())

    atoms: Dict[str, Dict[str, float]] = {}
    for atom in mol.GetAtoms():
        idx = atom.GetIdx()
        pt = drawer.GetDrawCoords(idx)
        atoms[str(idx)] = {"x": float(pt.x), "y": float(pt.y)}

    bonds: Dict[str, Dict[str, Dict[str, float]]] = {}
    for bond in mol.GetBonds():
        b_idx = bond.GetIdx()
        start_idx = bond.GetBeginAtomIdx()
        end_idx = bond.GetEndAtomIdx()
        bonds[str(b_idx)] = {
            "start": atoms[str(start_idx)],
            "end": atoms[str(end_idx)]
        }
    canvas = {"width": width, "height": height}
    return PredictResponseData(
        smiles_canonical=smiles_canonical,
        image_svg=image_svg,
        canvas=canvas,
        atoms=atoms,
        bonds=bonds,
        molecule_id=mol_id
    )

def predict_single_controller(request: PredictSingleRequest) -> PredictSingleResponseData:
    bond = PredictedBond(
        idx=request.bond_idx, bde=113.7,
        begin_atom_idx=0, end_atom_idx=1,
        bond_atoms="C-O"
    )
    return PredictSingleResponseData(
        smiles_canonical=request.smiles,
        bond=bond
    )

def predict_multiple_controller(request: PredictMultipleRequest) -> PredictMultipleResponseData:
    bonds = [
        PredictedBond(idx=1, bde=112.3, begin_atom_idx=0, end_atom_idx=1, bond_atoms="C-O"),
        PredictedBond(idx=2, bde=110.5, begin_atom_idx=1, end_atom_idx=2, bond_atoms="C-C")
    ]
    return PredictMultipleResponseData(
        smiles="CCO",
        molecule_id=request.molecule_id,
        bonds=bonds
    )

def fragment_controller(request: FragmentRequest) -> FragmentResponseData:
    bonds = [
        EvaluatedFragmentBond(idx=0, begin_atom=0, end_atom=1, bond_atoms="C-C", is_fragmentable=True)
    ]
    smiles_list = ["CCO", "CC", "O"] if request.export_smiles else None
    xyz_block = "...XYZ data..." if request.export_xyz else None
    return FragmentResponseData(
        smiles_canonical="CCO",
        molecule_id=request.molecule_id,
        bonds=bonds,
        smiles_list=smiles_list,
        xyz_block=xyz_block
    )

def predict_check_controller(request: PredictCheckRequest) -> PredictCheckResponseData:
    bond = PredictedBond(
        idx=request.bond_idx, bde=113.7,
        begin_atom_idx=0, end_atom_idx=1,
        bond_atoms="C-O"
    )
    return PredictCheckResponseData(
        smiles_canonical=request.smiles,
        bond=bond,
        products=request.products
    )

def infer_all_controller(request: InferAllRequest) -> InferAllResponseData:
    bonds = [
        PredictedBond(idx=0, bde=112.3, begin_atom_idx=0, end_atom_idx=1, bond_atoms="C-O"),
        PredictedBond(idx=1, bde=110.5, begin_atom_idx=1, end_atom_idx=2, bond_atoms="C-C")
    ]
    return InferAllResponseData(
        smiles_canonical=request.smiles,
        molecule_id="a1b2c3d4e5f6a7b8",
        bonds=bonds
    )

def download_report_controller(request: DownloadReportRequest) -> DownloadReportResponseData:
    return DownloadReportResponseData(
        report_base64=request.smiles + '_report_base64_string',
    )