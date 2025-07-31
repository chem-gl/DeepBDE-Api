from api.model.dto import (
    PredictRequest, PredictResponse,
    PredictSingleRequest, PredictSingleResponse,
    PredictMultipleRequest, PredictMultipleResponse, PredictMultipleBond
) 
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
            - smiles_canonical (str): SMILES canónico de la molécula ingresada.
            - image_svg (str): Imagen SVG enriquecida con la información de átomos y enlaces.
            - canvas (dict): Metadatos del lienzo (width, height).
            - atoms (dict): Posiciones de los átomos.
            - bonds (dict): Posiciones de los enlaces (inicio y fin de cada bond).
            - mol_id (str): ID único asociado al SMILES canónico.
    """
    mol = Chem.MolFromSmiles(request.smiles)
    if mol is None:
        raise ValueError(f"SMILES inválido: {request.smiles}")
    try:
        Chem.SanitizeMol(mol)
    except (Chem.MolSanitizeException, Chem.KekulizeException) as e:
        raise ValueError(f"Error al sanitizar molécula: {e}")
    smiles_canonical = Chem.MolToSmiles(mol, canonical=True)
    mol_id = hashlib.sha256(smiles_canonical.encode()).hexdigest()[:16]
    try:
        from rdkit.Chem import rdCoordGen
        rdCoordGen.AddCoords(mol)
    except ImportError:
        AllChem.Compute2DCoords(mol)
    n_atoms = mol.GetNumAtoms()
    dimension = max(300, min(100 + 20 * n_atoms, 800))
    width, height = dimension, dimension
    drawer = rdMolDraw2D.MolDraw2DSVG(width, height)
    opts = drawer.drawOptions()
    opts.addStereoAnnotation = True
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
    # Validar y canonizar SMILES
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
    return PredictSingleResponse(bde=432.1)  # Mocked value for demonstration

def predict_multiple_bde_controller(request: PredictMultipleRequest) -> PredictMultipleResponse:
    """
    Controlador para la predicción de BDE de múltiples enlaces de una molécula.

    Parámetros:
        request (PredictMultipleRequest):
            - smiles (str): SMILES de la molécula.
            - bond_idxs (List[int]): Lista de índices de enlaces a predecir.
    Devuelve:
        PredictMultipleResponse:
            - bonds (List[PredictMultipleBond]): Lista de enlaces con índice y BDE predicho.
    """
    # MOCK
    return PredictMultipleResponse(bonds=[
        PredictMultipleBond(idx=1, bde=112.3),
        PredictMultipleBond(idx=2, bde=110.5)
    ])
