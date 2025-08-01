from architecture import model
import torch
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
from deepbde.architecture.inference_util import multi_predict, single_predict
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

def is_single_bond(mol: Chem.Mol, bond_idx: int) -> bool:
    """
    Checks if the given bond index corresponds to a single bond in the molecule.
    Args:
        mol (Chem.Mol): The RDKit molecule object.
        bond_idx (int): The index of the bond to check.
    Returns:
        bool: True if the bond is a single bond, False otherwise.
    """
    if bond_idx < 0 or bond_idx >= mol.GetNumBonds():
        raise ValueError(f"Bond index {bond_idx} out of range for molecule with {mol.GetNumBonds()} bonds")
    bond = mol.GetBondWithIdx(bond_idx)
    return bond.GetBondType() == Chem.BondType.SINGLE

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

    from api.model.dto import Atom2D, Bond2D
    atoms: Dict[str, Atom2D] = {}
    for atom in mol.GetAtoms():
        idx = atom.GetIdx()
        pt = drawer.GetDrawCoords(idx)
        symbol = atom.GetSymbol()
        # Generate a simple atom SMILES (for most atoms, just the symbol; for special cases, use RDKit)
        atom_smiles = Chem.MolToSmiles(Chem.MolFromSmiles(f'[{symbol}]')) if symbol else None
        atoms[str(idx)] = Atom2D(
            x=float(pt.x),
            y=float(pt.y),
            symbol=symbol,
            smiles=atom_smiles
        )

    bonds: Dict[str, Bond2D] = {}
    for bond in mol.GetBonds():
        b_idx = bond.GetIdx()
        start_idx = bond.GetBeginAtomIdx()
        end_idx = bond.GetEndAtomIdx()
        start_pt = drawer.GetDrawCoords(start_idx)
        end_pt = drawer.GetDrawCoords(end_idx)
        atom1 = mol.GetAtomWithIdx(start_idx)
        atom2 = mol.GetAtomWithIdx(end_idx)
        bond_atoms = f"{atom1.GetSymbol()}-{atom2.GetSymbol()}"
        bonds[str(b_idx)] = Bond2D(
            start=start_idx,
            end=end_idx,
            start_coords={"x": float(start_pt.x), "y": float(start_pt.y)},
            end_coords={"x": float(end_pt.x), "y": float(end_pt.y)},
            bond_atoms=bond_atoms
        )
    canvas = {"width": width, "height": height}
    return PredictResponseData(
        smiles_canonical=smiles_canonical,
        image_svg=image_svg,
        canvas=canvas,
        atoms=atoms,
        bonds=bonds,
        molecule_id=mol_id
    )


# Controlador para la predicción de un enlace específico
# /api/v1/predict/single

def predict_single_controller(request: PredictSingleRequest) -> PredictSingleResponseData:
    if not verify_smiles(request.smiles, request.molecule_id):
        raise ValueError("SMILES verification failed")
    mol, _, smiles_canonical = generate_id_smiles(request.smiles)
    if request.bond_idx < 0 or request.bond_idx >= mol.GetNumBonds():
        raise ValueError(f"Bond index {request.bond_idx} out of range for molecule with {mol.GetNumBonds()} bonds")
    bond = mol.GetBondWithIdx(request.bond_idx)
    if bond is None:
        raise ValueError(f"Bond index {request.bond_idx} does not exist in the molecule")
    if not is_single_bond(mol, request.bond_idx):
        raise ValueError(f"Bond index {request.bond_idx} is not a single bond")
    begin_atom_idx = bond.GetBeginAtomIdx()
    end_atom_idx = bond.GetEndAtomIdx()
    atom1 = mol.GetAtomWithIdx(begin_atom_idx)
    atom2 = mol.GetAtomWithIdx(end_atom_idx)
    bond_atoms = f"{atom1.GetSymbol()}-{atom2.GetSymbol()}"
    import torch
    with torch.no_grad():
        pred = single_predict(smiles_canonical, request.bond_idx)
        if hasattr(pred, 'numpy'):
            pred = pred.numpy()
        bde = float(pred) if not hasattr(pred, '__len__') or len(pred) == 1 else float(pred[0])
    predicted_bond = PredictedBond(
        idx=request.bond_idx,
        bde=bde,
        begin_atom_idx=begin_atom_idx,
        end_atom_idx=end_atom_idx,
        bond_atoms=bond_atoms
    )
    return PredictSingleResponseData(
        smiles_canonical=smiles_canonical,
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
    # Verificar SMILES y obtener molécula
    if not verify_smiles(request.smiles, request.molecule_id):
        raise ValueError("SMILES verification failed")
    mol, _, smiles_canonical = generate_id_smiles(request.smiles)
    bond_indices = request.bond_indices
    errors = []
    for idx in bond_indices:
        if idx < 0 or idx >= mol.GetNumBonds():
            errors.append(f"Bond index {idx} out of range for molecule with {mol.GetNumBonds()} bonds")
        else:
            bond = mol.GetBondWithIdx(idx)
            if bond is None:
                errors.append(f"Bond index {idx} does not exist in the molecule")
            elif not is_single_bond(mol, idx):
                errors.append(f"Bond index {idx} is not a single bond")
    if errors:
        logger.error(f"Índices de enlace inválidos: {errors}")
        raise ValueError(f"Invalid bond indices: {errors}")

    try:
        bde_array = multi_predict(smiles_canonical, bond_indices)
        bde_list = bde_array.numpy().tolist() if hasattr(bde_array, 'numpy') else list(bde_array)
    except Exception as e:
        logger.error(f"Error en multi_predict: {e}")
        raise ValueError(f"Error in multi_predict: {e}")
    
    bonds = []
    for idx, bde in zip(bond_indices, bde_list):
        bond = mol.GetBondWithIdx(idx)
        begin_atom_idx = bond.GetBeginAtomIdx()
        end_atom_idx = bond.GetEndAtomIdx()
        begin_atom = mol.GetAtomWithIdx(begin_atom_idx)
        end_atom = mol.GetAtomWithIdx(end_atom_idx)
        bond_atoms = f"{begin_atom.GetSymbol()}-{end_atom.GetSymbol()}"
        # Logging para depuración
        logger.info(f"Predicción BDE para idx={idx}: valor={bde} tipo={type(bde)}")
        # Si bde es una lista, tomar el primer elemento
        if isinstance(bde, list) and len(bde) > 0:
            bde = bde[0]
        try:
            bde_float = float(bde)
        except Exception as ex:
            logger.error(f"No se pudo convertir bde a float: valor={bde} tipo={type(bde)} error={ex}")
            raise ValueError(f"BDE value for bond {idx} is not a valid number: {bde}")
        bonds.append(PredictedBond(
            idx=idx,
            bde=bde_float,
            begin_atom_idx=begin_atom_idx,
            end_atom_idx=end_atom_idx,
            bond_atoms=bond_atoms
        ))

    return PredictMultipleResponseData(
        smiles=smiles_canonical,
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