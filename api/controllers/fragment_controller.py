from api.model.dto import (
    FragmentSmilesRequest, FragmentSmilesResponse,
    FragmentXYZRequest, FragmentXYZResponse
)
from typing import List, Dict, Any

def fragment_smiles_controller(request: FragmentSmilesRequest) -> FragmentSmilesResponse:
    return FragmentSmilesResponse(
        parent=request.smiles,
        fragments=["[CH3]", "[OH]"],
        all_fragments=[
            {"bond_idx": 0, "fragments": ["CC", "O"]},
            {"bond_idx": 1, "fragments": ["C", "CO"]}
        ]
    )

def fragment_xyz_controller(request: FragmentXYZRequest) -> FragmentXYZResponse:
    xyz = "3\n\nC 0.000 0.000 0.000\nC 1.200 0.000 0.000\nO 2.400 0.000 0.000\n\n1\n\nO 0.000 0.000 0.000"
    return FragmentXYZResponse(xyz=xyz)
