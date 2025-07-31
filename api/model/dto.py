
# DTOs (Data Transfer Objects) para la API de DeepBDE
# Cada clase representa la estructura de datos de entrada o salida de los endpoints.
from dataclasses import dataclass
from typing import List, Optional, Dict, Any


@dataclass
class PredictRequest:
    """
    Request para /predict/ (imagen 2D y enlaces):
    - smiles: SMILES de la molécula
    """
    smiles: str


@dataclass
class PredictResponseBond:
    """
    Representa un enlace de la molécula:
    - idx: índice del enlace (RDKit)
    - atoms: índices de los átomos que conecta
    - bde: energía de disociación (opcional, solo en endpoints de predicción)
    """
    idx: int
    atoms: List[int]
    bde: Optional[float] = None


@dataclass
class PredictResponse:
    """
    Respuesta de /predict/:
    - image: imagen 2D en base64
    - bonds: lista de enlaces (sin BDE)
    """
    image: str
    bonds: List[PredictResponseBond]


@dataclass
class PredictSingleRequest:
    """
    Request para /predict/single/:
    - smiles: SMILES de la molécula
    - bond_idx: índice del enlace a predecir
    """
    smiles: str
    bond_idx: int


@dataclass
class PredictSingleResponse:
    """
    Respuesta de /predict/single/:
    - bde: energía de disociación predicha
    """
    bde: float


@dataclass
class PredictMultipleRequest:
    """
    Request para /predict/multiple/:
    - smiles: SMILES de la molécula
    - bond_indices: lista de índices de enlaces a predecir
    """
    smiles: str
    bond_indices: List[int]


@dataclass
class PredictMultipleBond:
    """
    Representa un enlace y su BDE predicha (para /predict/multiple/):
    - idx: índice del enlace
    - bde: energía de disociación
    """
    idx: int
    bde: float


@dataclass
class PredictMultipleResponse:
    """
    Respuesta de /predict/multiple/:
    - bonds: lista de enlaces y sus BDEs
    """
    bonds: List[PredictMultipleBond]


@dataclass
class FragmentSmilesRequest:
    """
    Request para /fragment/smiles/:
    - smiles: SMILES de la molécula
    - all_bonds: si True, devuelve todos los fragmentos posibles
    """
    smiles: str
    all_bonds: Optional[bool] = False


@dataclass
class FragmentSmilesResponse:
    """
    Respuesta de /fragment/smiles/:
    - parent: SMILES original
    - fragments: lista de fragmentos principales
    - all_fragments: lista de todos los fragmentos por enlace
    """
    parent: str
    fragments: List[str]
    all_fragments: List[Dict[str, Any]]


@dataclass
class FragmentXYZRequest:
    """
    Request para /fragment/xyz/:
    - smiles: SMILES de la molécula
    - all_bonds: si True, devuelve todos los fragmentos posibles
    """
    smiles: str
    all_bonds: Optional[bool] = False


@dataclass
class FragmentXYZResponse:
    """
    Respuesta de /fragment/xyz/:
    - xyz: string con coordenadas XYZ
    """
    xyz: str


@dataclass
class PredictCheckRequest:
    """
    Request para /predict/check/:
    - smiles: SMILES de la molécula
    - bond_idx: índice del enlace
    - products: lista de productos esperados
    """
    smiles: str
    bond_idx: int
    products: List[str]


@dataclass
class PredictCheckResponse:
    """
    Respuesta de /predict/check/:
    - bde: energía de disociación predicha
    """
    bde: float


@dataclass
class InferAllRequest:
    """
    Request para /infer/all/:
    - smiles: SMILES de la molécula
    """
    smiles: str


@dataclass
class InferAllBond:
    """
    Representa un enlace y su BDE predicha (para /infer/all/):
    - idx: índice del enlace
    - bde: energía de disociación
    """
    idx: int
    bde: float


@dataclass
class InferAllResponse:
    """
    Respuesta de /infer/all/:
    - bonds: lista de enlaces y sus BDEs
    """
    bonds: List[InferAllBond]


@dataclass
class DownloadReportRequest:
    """
    Request para /download_report/:
    - smiles: SMILES de la molécula
    - bond_idx: índice del enlace
    - format: formato de reporte (ej: 'pdf')
    """
    smiles: str
    bond_idx: int
    format: str


@dataclass
class DownloadReportResponse:
    """
    Respuesta de /download_report/:
    - data: string con el reporte generado (ej: PDF en base64)
    """
    data: str
