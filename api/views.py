from rest_framework.response import Response
from drf_spectacular.utils import extend_schema, OpenApiExample, OpenApiResponse
from rest_framework.views import APIView
from pydantic import ValidationError

from .model.dto import (
    PredictRequest, PredictSingleRequest, PredictMultipleRequest,
    FragmentRequest, InferAllRequest, DownloadReportRequest, PredictCheckRequest,
    PredictResponseData, PredictSingleResponseData, PredictMultipleResponseData,
    FragmentResponseData, InferAllResponseData, DownloadReportResponseData, PredictCheckResponseData,
    PredictedBond, EvaluatedFragmentBond,
    ErrorDetail, ErrorCode, APIResponse
)

# Utilidad para una respuesta de error consistente
def _handle_validation_error(e: ValidationError, status_code=400):
    """Crea una respuesta de error estructurada a partir de un ValidationError."""
    first = e.errors()[0]
    loc = first['loc'][-1]
    code = ErrorCode.BOND_INDEX_OUT_OF_RANGE if loc == 'bond_idx' else ErrorCode.SMILES_INVALID
    detail = ErrorDetail(code=code, message=first['msg'])
    payload = APIResponse[None](status="error", error=detail).model_dump()
    return Response(payload, status=status_code)

class PredictView(APIView):
    """
    Endpoint para obtener la imagen 2D de la molécula y la lista de enlaces con sus índices y átomos involucrados (sin predecir BDEs).
    """
    @extend_schema(
        description="""
        Devuelve información enriquecida de la molécula: SMILES canónico, imagen SVG, metadatos de canvas, posiciones de átomos y enlaces, id único.
        
        Returns enriched information about the molecule: canonical SMILES, SVG image, canvas metadata, atom and bond positions, unique ID.
        """,
        request=PredictRequest,
        responses={
            200: OpenApiResponse(
                response=APIResponse[PredictResponseData],
                description="Respuesta exitosa / Successful response"
            ),
            400: OpenApiResponse(
                response=APIResponse[None],
                description="Error de validación / Validation error"
            )
        },
        examples=[
            OpenApiExample(
                'Ejemplo de entrada / Example input',
                value={"smiles": "CCO"},
                request_only=True
            ),
            OpenApiExample(
                'Ejemplo alternativo / Alternative example',
                value={"smiles": "C1=CC=CC=C1"},
                request_only=True
            )
        ]
    )
    def post(self, request):
        try:
            # Validación de entrada (con Pydantic)
            req = PredictRequest(**request.data)

            # Lógica de negocio simulada
            resp_data = PredictResponseData(
                smiles_canonical=req.smiles,
                image_svg="<svg>...</svg>",
                canvas={"width": 300, "height": 300},
                atoms={"0": {"x": 123.45, "y": 67.89}},
                bonds={"0": {"start": {"x": 123.45, "y": 67.89}, "end": {"x": 156.78, "y": 120.34}}},
                molecule_id="a1b2c3d4e5f6a7b8"
            )

            resp = APIResponse[PredictResponseData](status="success", data=resp_data)
            return Response(resp.model_dump())
        except ValidationError as e:
            return _handle_validation_error(e)

class PredictSingleView(APIView):
    """
    Endpoint para predecir la energía de disociación de un enlace específico de la molécula (SMILES y bond_idx).
    """
    @extend_schema(
        description="""
        Predice la energía de disociación para un enlace específico de la molécula.
        
        Predicts the dissociation energy for a specific bond of the molecule.
        """,
        request=PredictSingleRequest,
        responses={
            200: OpenApiResponse(
                response=APIResponse[PredictSingleResponseData],
                description="Respuesta exitosa / Successful response"
            ),
            400: OpenApiResponse(
                response=APIResponse[None],
                description="Error de validación / Validation error"
            )
        },
        examples=[
            OpenApiExample(
                "Ejemplo de entrada / Example input",
                value={"smiles": "CCO", "bond_idx": 1},
                request_only=True
            ),
            OpenApiExample(
                "Ejemplo alternativo / Alternative example",
                value={"smiles": "C1=CC=CC=C1", "bond_idx": 2},
                request_only=True
            )
        ]
    )
    def post(self, request):
        try:
            req = PredictSingleRequest(**request.data)
            bond = PredictedBond(
                idx=req.bond_idx, bde=113.7,
                begin_atom_idx=0, end_atom_idx=1,
                bond_atoms="C-O"
            )
            resp_data = PredictSingleResponseData(
                smiles_canonical=req.smiles,
                bond=bond
            )
            resp = APIResponse[PredictSingleResponseData](status="success", data=resp_data)
            return Response(resp.model_dump())
        except ValidationError as e:
            return _handle_validation_error(e)

class PredictMultipleView(APIView):
    """
    Endpoint para predecir las energías de disociación para varios enlaces de la molécula (SMILES y bond_indices).
    """
    @extend_schema(
        description="""
        Predice las energías de disociación para varios enlaces de la molécula.
        
        Predicts the dissociation energies for multiple bonds of the molecule.
        """,
        request=PredictMultipleRequest,
        responses={
            200: OpenApiResponse(
                response=APIResponse[PredictMultipleResponseData],
                description="Respuesta exitosa / Successful response"
            ),
            400: OpenApiResponse(
                response=APIResponse[None],
                description="Error de validación / Validation error"
            )
        },
        examples=[
            OpenApiExample(
                "Entrada de ejemplo / Example input",
                value={"mol_id": "a1b2c3d4e5f6a7b8", "bond_indices": [1, 2]},
                request_only=True
            ),
            OpenApiExample(
                "Entrada alternativa / Alternative input",
                value={"mol_id": "h1g2f3e4d5c6b7a8", "bond_indices": [0, 2, 3]},
                request_only=True
            )
        ]
    )
    def post(self, request):
        try:
            req = PredictMultipleRequest(**request.data)
            bonds = [
                PredictedBond(idx=1, bde=112.3, begin_atom_idx=0, end_atom_idx=1, bond_atoms="C-O"),
                PredictedBond(idx=2, bde=110.5, begin_atom_idx=1, end_atom_idx=2, bond_atoms="C-C")
            ]
            resp_data = PredictMultipleResponseData(
                smiles="CCO",
                molecule_id=req.molecule_id,
                bonds=bonds
            )
            resp = APIResponse[PredictMultipleResponseData](status="success", data=resp_data)
            return Response(resp.model_dump())
        except ValidationError as e:
            return _handle_validation_error(e)

class FragmentView(APIView):
    """
    Unified endpoint to generate molecular fragments in SMILES or XYZ format.
    """
    @extend_schema(
        description="""
        Genera fragmentos moleculares en formato SMILES o XYZ.
        
        Generates molecular fragments in SMILES or XYZ format.
        """,
        request=FragmentRequest,
        responses={
            200: OpenApiResponse(
                response=APIResponse[FragmentResponseData],
                description="Respuesta exitosa / Successful response"
            ),
            400: OpenApiResponse(
                response=APIResponse[None],
                description="Error de validación / Validation error"
            )
        },
        examples=[
            OpenApiExample(
                "Entrada de ejemplo / Example input",
                value={"smiles": "CCO", "molecule_id": "a1b2c3d4e5f6a7b8", "export_smiles": True},
                request_only=True
            ),
            OpenApiExample(
                "Entrada alternativa / Alternative input",
                value={"smiles": "C1=CC=CC=C1", "molecule_id": "h1g2f3e4d5c6b7a8", "export_xyz": True},
                request_only=True
            )
        ]
    )
    def post(self, request):
        try:
            req = FragmentRequest(**request.data)
            bonds = [
                EvaluatedFragmentBond(idx=0, begin_atom=0, end_atom=1, bond_atoms="C-C", is_fragmentable=True)
            ]
            smiles_list = ["CCO", "CC", "O"] if req.export_smiles else None
            xyz_block = "...XYZ data..." if req.export_xyz else None
            resp_data = FragmentResponseData(
                smiles_canonical="CCO",
                molecule_id=req.molecule_id,
                bonds=bonds,
                smiles_list=smiles_list,
                xyz_block=xyz_block
            )
            resp = APIResponse[FragmentResponseData](status="success", data=resp_data)
            return Response(resp.model_dump())
        except ValidationError as e:
            return _handle_validation_error(e)

class PredictCheckView(APIView):
    """
    Endpoint para verificar si los productos generados por la escisión de un enlace corresponden a los esperados y predecir la BDE.
    """
    @extend_schema(
        description="""
        Verifica productos generados por la escisión de un enlace y predice la BDE.
        
        Verifies products generated by the cleavage of a bond and predicts the BDE.
        """,
        request=PredictCheckRequest,
        responses={
            200: OpenApiResponse(
                response=APIResponse[PredictCheckResponseData],
                description="Respuesta exitosa / Successful response"
            ),
            400: OpenApiResponse(
                response=APIResponse[None],
                description="Error de validación / Validation error"
            )
        },
        examples=[
            OpenApiExample(
                "Entrada de ejemplo / Example input",
                value={"smiles": "CCO", "bond_idx": 1, "products": ["CC", "O"]},
                request_only=True
            ),
            OpenApiExample(
                "Entrada alternativa / Alternative input",
                value={"smiles": "C1=CC=CC=C1", "bond_idx": 2, "products": ["C1=CC", "C=C1"]},
                request_only=True
            )
        ]
    )
    def post(self, request):
        try:
            req = PredictCheckRequest(**request.data)
            bond = PredictedBond(
                idx=req.bond_idx, bde=113.7,
                begin_atom_idx=0, end_atom_idx=1,
                bond_atoms="C-O"
            )
            resp_data = PredictCheckResponseData(
                smiles_canonical=req.smiles,
                bond=bond,
                products=req.products
            )
            resp = APIResponse[PredictCheckResponseData](status="success", data=resp_data)
            return Response(resp.model_dump())
        except ValidationError as e:
            return _handle_validation_error(e)

class InferAllView(APIView):
    """
    Endpoint para predecir las energías de disociación para todos los enlaces simples de la molécula dada.
    """
    @extend_schema(
        description="""
        Predice las energías de disociación para todos los enlaces simples de la molécula dada.
        
        Predicts the dissociation energies for all single bonds of the given molecule.
        """,
        request=InferAllRequest,
        responses={
            200: OpenApiResponse(
                response=APIResponse[InferAllResponseData],
                description="Respuesta exitosa / Successful response"
            ),
            400: OpenApiResponse(
                response=APIResponse[None],
                description="Error de validación / Validation error"
            )
        },
        examples=[
            OpenApiExample(
                "Entrada de ejemplo / Example input",
                value={"smiles": "CCO"},
                request_only=True
            ),
            OpenApiExample(
                "Entrada alternativa / Alternative input",
                value={"smiles": "C1=CC=CC=C1"},
                request_only=True
            )
        ]
    )
    def post(self, request):
        try:
            req = InferAllRequest(**request.data)
            bonds = [
                PredictedBond(idx=0, bde=112.3, begin_atom_idx=0, end_atom_idx=1, bond_atoms="C-O"),
                PredictedBond(idx=1, bde=110.5, begin_atom_idx=1, end_atom_idx=2, bond_atoms="C-C")
            ]
            resp_data = InferAllResponseData(
                smiles_canonical=req.smiles,
                molecule_id="a1b2c3d4e5f6a7b8",
                bonds=bonds
            )
            resp = APIResponse[InferAllResponseData](status="success", data=resp_data)
            return Response(resp.model_dump())
        except ValidationError as e:
            return _handle_validation_error(e)

class DownloadReportView(APIView):
    """
    Endpoint para generar y descargar un informe PDF con los resultados de la predicción para la molécula y enlace indicados.
    """
    @extend_schema(
        description="""
        Genera y descarga un informe PDF con los resultados de la predicción para la molécula y enlace indicados.
        
        Generates and downloads a PDF report with the prediction results for the indicated molecule and bond.
        """,
        request=DownloadReportRequest,
        responses={
            200: OpenApiResponse(
                response=APIResponse[DownloadReportResponseData],
                description="Respuesta exitosa / Successful response"
            ),
            400: OpenApiResponse(
                response=APIResponse[None],
                description="Error de validación / Validation error"
            )
        },
        examples=[
            OpenApiExample(
                "Entrada de ejemplo / Example input",
                value={"smiles": "CCO", "bond_idx": 1, "format": "pdf"},
                request_only=True
            ),
            OpenApiExample(
                "Entrada alternativa / Alternative input",
                value={"smiles": "C1=CC=CC=C1", "bond_idx": 2, "format": "pdf"},
                request_only=True
            )
        ]
    )
    def post(self, request):
        try:
            DownloadReportRequest(**request.data)  # Validate input
            resp_data = DownloadReportResponseData(
                report_base64="...PDF data in base64..."
            )
            resp = APIResponse[DownloadReportResponseData](status="success", data=resp_data)
            return Response(resp.model_dump())
        except ValidationError as e:
            return _handle_validation_error(e)