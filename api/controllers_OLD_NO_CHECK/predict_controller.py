 
import torch
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.Draw import rdMolDraw2D
import hashlib
from typing import Dict

from deepbde.architecture.inference_util import single_predict
def predict_bde_controller(request: PredictRequest) -> PredictResponse:
    """
    Controlador para la predicción general de enlaces de una molécula.
    Parámetros:
        request (PredictRequest):
            - smiles (str): SMILES de la molécula.
    Devuelve:
        PredictResponse:
            - smiles_canonical (str): SMILES canónico con hidrógenos explícitos.
            - image_svg (str): Imagen SVG enriquecida con información atómica.
            - canvas (dict): Dimensiones del lienzo.
            - atoms (dict): Coordenadas de los átomos.
            - bonds (dict): Coordenadas de los enlaces.
            - mol_id (str): ID único de la molécula.
    """
    mol = Chem.MolFromSmiles(request.smiles)
    if mol is None:
        raise ValueError(f"SMILES inválido: {request.smiles}")
    try:
        Chem.SanitizeMol(mol)
    except (Chem.MolSanitizeException, Chem.KekulizeException) as e:
        raise ValueError(f"Error al sanitizar molécula: {e}")
    mol = Chem.AddHs(mol)
    smiles_canonical = Chem.MolToSmiles(mol, canonical=True, allHsExplicit=True)
    mol_id = hashlib.sha256(smiles_canonical.encode()).hexdigest()[:16]

    # Generar coordenadas 2D
    try:
        from rdkit.Chem import rdCoordGen
        rdCoordGen.AddCoords(mol)
    except ImportError:
        AllChem.Compute2DCoords(mol)

    # Determinar dimensiones del lienzo
    n_atoms = mol.GetNumAtoms()
    dimension = max(300, min(100 + 20 * n_atoms, 800))
    width, height = dimension, dimension

    # Dibujar molécula en SVG
    drawer = rdMolDraw2D.MolDraw2DSVG(width, height)
    opts = drawer.drawOptions()
    opts.addStereoAnnotation = True
    drawer.DrawMolecule(mol)
    drawer.FinishDrawing()
    raw_svg = drawer.GetDrawingText()
    image_svg = ''.join(line.strip() for line in raw_svg.splitlines())

    # Coordenadas de átomos
    atoms: Dict[str, Dict[str, float]] = {}
    for atom in mol.GetAtoms():
        idx = atom.GetIdx()
        pt = drawer.GetDrawCoords(idx)
        atoms[str(idx)] = {"x": float(pt.x), "y": float(pt.y)}

    # Coordenadas de enlaces
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
    return PredictResponse(
        smiles_canonical=smiles_canonical,
        image_svg=image_svg,
        canvas=canvas,
        atoms=atoms,
        bonds=bonds,
        mol_id=mol_id
    )




def predict_single_bde_controller(request: PredictSingleRequest) -> PredictSingleResponse:
    """
    Controlador para la predicción de BDE de un enlace específico.

    Parámetros:
        request (PredictSingleRequest):
            - smiles (str): SMILES de la molécula.
            - bond_idx (int): Índice del enlace a predecir.
    Devuelve:
        PredictSingleResponse:
            - bde (float): Valor de BDE predicho para el enlace.
    """
    mol = Chem.MolFromSmiles(request.smiles)
    if mol is None:
        raise ValueError(f"SMILES inválido: {request.smiles}")
    try:
        Chem.SanitizeMol(mol)
    except (Chem.MolSanitizeException, Chem.KekulizeException) as e:
        raise ValueError(f"Error al sanitizar molécula: {e}")
    smiles_canonical = Chem.MolToSmiles(mol, canonical=True)
    with torch.no_grad():
        pred = single_predict(smiles_canonical, request.bond_idx)
        if hasattr(pred, 'numpy'):
            pred = pred.numpy()
        bde = float(pred) if not hasattr(pred, '__len__') or len(pred) == 1 else float(pred[0])
    return PredictSingleResponse(bde=bde)

def predict_multiple_bde_controller(request: PredictMultipleRequest) -> PredictMultipleResponse:
    """
    Controlador para la predicción de BDE de múltiples enlaces de una molécula.
    Parámetros:
        request (PredictMultipleRequest):
            - smiles (str): SMILES de la molécula.
            - bond_indices (List[int]): Lista de índices de enlaces a predecir.
    Devuelve:
        PredictMultipleResponse:
            - bonds (List[PredictMultipleBond]): Lista de enlaces con índice y BDE predicho.
    """
    mol = Chem.MolFromSmiles(request.smiles)
    if mol is None:
        raise ValueError(f"SMILES inválido: {request.smiles}")
    try:
        Chem.SanitizeMol(mol)
    except (Chem.MolSanitizeException, Chem.KekulizeException) as e:
        raise ValueError(f"Error al sanitizar molécula: {e}")
    smiles_canonical = Chem.MolToSmiles(mol, canonical=True)
    num_bonds = mol.GetNumBonds()
    valid_indices = [idx for idx in request.bond_indices if 0 <= idx < num_bonds]
    if not valid_indices:
        return PredictMultipleResponse(bonds=[])
    bonds = []
    try:
        for idx in valid_indices:
            with torch.no_grad():
                pred = single_predict(smiles_canonical, idx)
                if hasattr(pred, 'numpy'):
                    pred = pred.numpy()
                bde = float(pred) if not hasattr(pred, '__len__') or len(pred) == 1 else float(pred[0])
                bonds.append(PredictBond(idx=idx, bde=bde))
        if not bonds or len(bonds) != len(valid_indices):
            raise RuntimeError("The model cannot process the prediction for the requested smile or bonds.")
        return PredictMultipleResponse(bonds=bonds)
    except Exception:
        raise RuntimeError("The model cannot process the prediction for the requested  smile or bonds.")
