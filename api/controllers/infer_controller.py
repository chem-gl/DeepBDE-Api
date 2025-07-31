from api.model.dto import InferAllRequest, InferAllResponse, InferAllBond
import numpy as np

 

def infer_all_controller(request: InferAllRequest) -> InferAllResponse:
    """
    Realiza la inferencia de BDE para todos los enlaces válidos usando DeepBDE y RDKit.
    """
    smiles = request.smiles
    bond_idxs, preds = predict_all(smiles)
    if hasattr(preds, 'numpy'):
        preds = preds.numpy()
    bonds = [InferAllBond(idx=int(idx), bde=float(bde)) for idx, bde in zip(bond_idxs, preds)]
    return InferAllResponse(bonds=bonds)
