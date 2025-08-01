from api.model.dto import PredictCheckRequest, PredictCheckResponseData, BondDetails
def predict_check_controller(request: PredictCheckRequest) -> PredictCheckResponseData:
    # Mock: always returns a valid response
    return PredictCheckResponseData(
        smiles_canonical=request.smiles,
        bond=BondDetails(
            idx=request.bond_idx,
            bde=113.7,
            begin_atom_idx=0,
            end_atom_idx=1,
            bond_atoms="C-O"
        ),
        products=request.products
    )
