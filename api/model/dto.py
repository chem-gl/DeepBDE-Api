from typing import List, Optional, Dict, Any, Literal, TypeVar, Generic
from enum import Enum
from pydantic import BaseModel, Field, field_validator, model_validator, ValidationError

# ---------- Enumeración de códigos de error ----------
class ErrorCode(str, Enum):
    """Enumeración de códigos de error estándar de la API."""
    SMILES_INVALID = "SMILES_INVALID"
    BOND_INDEX_OUT_OF_RANGE = "BOND_INDEX_OUT_OF_RANGE"
    UNSUPPORTED_FORMAT = "UNSUPPORTED_FORMAT"
    MISSING_EXPORT_OPTION = "MISSING_EXPORT_OPTION"

# ---------- Clases comunes ----------
class ErrorDetail(BaseModel):
    """
    Detalle de error estándar para respuestas de la API.
    - code: Código de error único (ErrorCode).
    - message: Descripción amigable del error.
    """
    code: ErrorCode  # Código de error único, basado en la enumeración ErrorCode.
    message: str  # Descripción amigable del error.

class PredictedBond(BaseModel):
    """
    Detalle de un enlace químico y su BDE predicha.
    - idx: Índice del enlace en la molécula.
    - bde: Energía de disociación predicha para el enlace.
    - begin_atom_idx: Índice del átomo inicial del enlace.
    - end_atom_idx: Índice del átomo final del enlace.
    - bond_atoms: Representación legible del enlace (ejemplo: 'C-O').
    """
    idx: int  # Índice del enlace en la molécula.
    bde: float  # Energía de disociación predicha para el enlace.
    begin_atom_idx: int  # Índice del átomo inicial del enlace.
    end_atom_idx: int  # Índice del átomo final del enlace.
    bond_atoms: str  # Representación legible del enlace (ejemplo: 'C-O').

class EvaluatedFragmentBond(BaseModel):
    """
    Representa un enlace evaluado para fragmentación.
    - idx: Índice del enlace en la molécula.
    - begin_atom: Índice del átomo inicial del enlace.
    - end_atom: Índice del átomo final del enlace.
    - bond_atoms: Representación legible del enlace (ejemplo: 'C-H').
    - is_fragmentable: Indica si el enlace puede fragmentarse (enlace simple, etc.).
    """
    idx: int  # Índice del enlace en la molécula.
    begin_atom: int  # Índice del átomo inicial del enlace.
    end_atom: int  # Índice del átomo final del enlace.
    bond_atoms: str  # Representación legible del enlace (ejemplo: 'C-H').
    is_fragmentable: bool  # Indica si el enlace puede fragmentarse (enlace simple, etc.).

# ---------- Plantilla de respuesta API ----------
T = TypeVar('T')

class APIResponse(BaseModel, Generic[T]):
    """
    Respuesta unificada para todos los endpoints.
    - status: 'success' o 'error'
    - data: Contenido de la respuesta
    - error: Detalle del error (opcional)
    """
    status: Literal["success", "error"] = "success"
    data: Optional[T] = None
    error: Optional[ErrorDetail] = None

# ---------- Modelos de endpoints ----------

# --- /predict/ ---
class PredictRequest(BaseModel):
    """
    Entrada para /predict/.
    - smiles: SMILES de la molécula
    """
    smiles: str

class PredictResponseData(BaseModel):
    """
    Salida de /predict/.
    - smiles_canonical: SMILES canónico con hidrógenos explícitos
    - image_svg: Imagen SVG enriquecida
    - canvas: Metadatos del lienzo
    - atoms: Posiciones de átomos
    - bonds: Posiciones de enlaces
    - molecule_id: ID único (hash SHA256)
    """
    smiles_canonical: str
    image_svg: str
    canvas: Dict[str, int]
    atoms: Dict[str, Dict[str, float]]
    bonds: Dict[str, Dict[str, Dict[str, float]]]
    molecule_id: str

# --- /predict/single/ ---
class PredictSingleRequest(BaseModel):
    """
    Entrada para /predict/single/.
    - smiles: SMILES de la molécula
    - bond_idx: Índice del enlace (>= 0)
    """
    smiles: str
    bond_idx: int

    @field_validator('bond_idx')
    @classmethod
    def _validate_bond_idx(cls, v: int) -> int:
        if v < 0:
            raise ValueError("bond_idx debe ser >= 0")
        return v

class PredictSingleResponseData(BaseModel):
    """
    Salida de /predict/single/.
    - smiles_canonical: SMILES canónico
    - bond: Detalle del enlace predicho
    """
    smiles_canonical: str
    bond: PredictedBond

# --- /predict/multiple/ ---
class PredictMultipleRequest(BaseModel):
    """
    Entrada para /predict/multiple/.
    - molecule_id: ID de la molécula
    - bond_indices: Lista de índices de enlaces
    - smiles: (opcional)
    """
    molecule_id: str
    bond_indices: List[int] = Field(..., min_items=1)
    smiles: Optional[str] = None

class PredictMultipleResponseData(BaseModel):
    """
    Salida de /predict/multiple/.
    - smiles: SMILES canónico
    - molecule_id: ID de la molécula
    - bonds: Lista de detalles de enlaces
    """
    smiles: str
    molecule_id: str
    bonds: List[PredictedBond]

# --- /fragment/ ---
class FragmentRequest(BaseModel):
    """
    Entrada para /fragment/.
    - smiles: SMILES con hidrógenos explícitos
    - molecule_id: ID de la molécula
    - bond_idx: Índice del enlace (opcional)
    - export_smiles: Si True, devuelve lista de SMILES
    - export_xyz: Si True, devuelve bloque XYZ
    """
    smiles: str
    molecule_id: str
    bond_idx: Optional[int] = None
    export_smiles: bool = False
    export_xyz: bool = False

    @model_validator(mode='before')
    @classmethod
    def _validate_export_options(cls, data: Any) -> Any:
        if not data.get('export_smiles') and not data.get('export_xyz'):
            raise ValueError("Debe seleccionarse al menos export_smiles o export_xyz")
        return data

class FragmentResponseData(BaseModel):
    """
    Salida de /fragment/.
    - smiles_canonical: SMILES canónico
    - molecule_id: ID de la molécula
    - bonds: Lista de enlaces evaluados
    - smiles_list: Lista de SMILES (opcional)
    - xyz_block: Cadena XYZ (opcional)
    """
    smiles_canonical: str
    molecule_id: str
    bonds: List[EvaluatedFragmentBond]
    smiles_list: Optional[List[str]] = None
    xyz_block: Optional[str] = None

# --- /predict/check/ ---
class PredictCheckRequest(BaseModel):
    """
    Entrada para /predict/check/.
    - smiles: SMILES original
    - bond_idx: Índice del enlace
    - products: Lista de SMILES esperados
    """
    smiles: str
    bond_idx: int
    products: List[str] = Field(..., min_items=1)

class PredictCheckResponseData(BaseModel):
    """
    Salida de /predict/check/.
    - smiles_canonical: SMILES canónico
    - bond: Detalle del enlace analizado
    - products: Lista de SMILES resultantes
    """
    smiles_canonical: str
    bond: PredictedBond
    products: List[str]

# --- /infer/all/ ---
class InferAllRequest(BaseModel):
    """
    Entrada para /infer/all/.
    - smiles: SMILES de la molécula
    """
    smiles: str

class InferAllResponseData(BaseModel):
    """
    Salida de /infer/all/.
    - smiles_canonical: SMILES canónico
    - molecule_id: ID de la molécula
    - bonds: Lista de detalles de enlaces
    """
    smiles_canonical: str
    molecule_id: str
    bonds: List[PredictedBond]

# --- /download_report/ ---
class DownloadReportRequest(BaseModel):
    """
    Entrada para /download_report/.
    - smiles: SMILES de la molécula
    - bond_idx: Índice del enlace
    - format: 'pdf' o 'html'
    """
    smiles: str
    bond_idx: int
    format: str

    @field_validator('format')
    @classmethod
    def _validate_format(cls, v: str) -> str:
        if v not in {"pdf", "html"}:
            raise ValueError("Formato no soportado: use 'pdf' o 'html'.")
        return v

class DownloadReportResponseData(BaseModel):
    """
    Salida de /download_report/.
    - report_base64: Reporte en base64
    """
    report_base64: str  # Base64 del reporte archivo
