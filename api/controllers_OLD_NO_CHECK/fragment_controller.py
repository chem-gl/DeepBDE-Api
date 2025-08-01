from typing import List, Dict, Any

from api.model.dto import (
    FragmentRequest, FragmentResponseData, FragmentBond
)

def fragment_smiles_controller(request: FragmentRequest) -> FragmentResponseData:
    # Mock: returns a valid FragmentResponseData
    bonds = [
        FragmentBond(idx=0, begin_atom=0, end_atom=1, bond_atoms="C-O", is_fragmentable=True),
        FragmentBond(idx=1, begin_atom=1, end_atom=2, bond_atoms="C-C", is_fragmentable=True)
    ]
    return FragmentResponseData(
        smiles_canonical=request.smiles,
        mol_id=request.mol_id,
        bonds=bonds,
        smiles_list=[request.smiles, "CC", "O"] if request.export_smiles else None,
        xyz_block=None
    )

def fragment_xyz_controller(request: FragmentRequest) -> FragmentResponseData:
    # Mock: returns a valid FragmentResponseData with xyz_block
    return FragmentResponseData(
        smiles_canonical=request.smiles,
        mol_id=request.mol_id,
        bonds=[],
        smiles_list=None,
        xyz_block="3\n\nC 0.000 0.000 0.000\nC 1.200 0.000 0.000\nO 2.400 0.000 0.000\n"
    )
