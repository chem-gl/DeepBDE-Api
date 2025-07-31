from dataclasses import dataclass
from typing import List, Optional, Dict, Any

@dataclass
class PredictRequest:
    smiles: str

@dataclass
class PredictResponseBond:
    idx: int
    atoms: List[int]
    bde: float

@dataclass
class PredictResponse:
    image: str
    bonds: List[PredictResponseBond]

@dataclass
class PredictSingleRequest:
    smiles: str
    bond_idx: int

@dataclass
class PredictSingleResponse:
    bde: float

@dataclass
class PredictMultipleRequest:
    smiles: str
    bond_indices: List[int]

@dataclass
class PredictMultipleBond:
    idx: int
    bde: float

@dataclass
class PredictMultipleResponse:
    bonds: List[PredictMultipleBond]

@dataclass
class FragmentSmilesRequest:
    smiles: str
    all_bonds: Optional[bool] = False

@dataclass
class FragmentSmilesResponse:
    parent: str
    fragments: List[str]
    all_fragments: List[Dict[str, Any]]

@dataclass
class FragmentXYZRequest:
    smiles: str
    all_bonds: Optional[bool] = False

@dataclass
class FragmentXYZResponse:
    xyz: str

@dataclass
class PredictCheckRequest:
    smiles: str
    bond_idx: int
    products: List[str]

@dataclass
class PredictCheckResponse:
    bde: float

@dataclass
class InferAllRequest:
    smiles: str

@dataclass
class InferAllBond:
    idx: int
    bde: float

@dataclass
class InferAllResponse:
    bonds: List[InferAllBond]

@dataclass
class DownloadReportRequest:
    smiles: str
    bond_idx: int
    format: str

@dataclass
class DownloadReportResponse:
    data: str
