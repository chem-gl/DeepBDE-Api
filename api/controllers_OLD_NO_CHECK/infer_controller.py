from api.model.dto import BondDetails, InferAllRequest, InferAllResponseData


def infer_all_controller(request: InferAllRequest) -> InferAllResponseData:
    # Mock: returns two bonds
    bonds = [
        BondDetails(idx=0, bde=110.2, begin_atom_idx=0, end_atom_idx=1, bond_atoms="C-O"),
        BondDetails(idx=1, bde=112.3, begin_atom_idx=1, end_atom_idx=2, bond_atoms="C-C")
    ]
    return InferAllResponseData(
        smiles_canonical=request.smiles,
        mol_id="a1b2c3d4e5f6a7b8",
        bonds=bonds
    )
