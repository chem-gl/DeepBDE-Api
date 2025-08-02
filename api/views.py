import logging
from rest_framework.response import Response
from drf_spectacular.utils import extend_schema, OpenApiExample, OpenApiResponse
from rest_framework.views import APIView
from pydantic import ValidationError

from api.controllers.predict_controller import (
    predict_controller, predict_single_controller, predict_multiple_controller,
    fragment_controller, predict_check_controller, infer_all_controller,
    download_report_controller
)

from .model.dto import (
    PredictRequest, PredictSingleRequest, PredictMultipleRequest,
    FragmentRequest, InferAllRequest, DownloadReportRequest, PredictCheckRequest,
    PredictResponseData, PredictSingleResponseData, PredictMultipleResponseData,
    FragmentResponseData, InferAllResponseData, DownloadReportResponseData, PredictCheckResponseData,
    ErrorDetail, ErrorCode, APIResponse
)

# Define aliases for APIResponse with specific data types
PredictResponse = APIResponse[PredictResponseData]
PredictSingleResponse = APIResponse[PredictSingleResponseData]
PredictMultipleResponse = APIResponse[PredictMultipleResponseData]
FragmentResponse = APIResponse[FragmentResponseData]
InferAllResponse = APIResponse[InferAllResponseData]
DownloadReportResponse = APIResponse[DownloadReportResponseData]
PredictCheckResponse = APIResponse[PredictCheckResponseData]
ErrorResponse = APIResponse[None]

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
    example1 = APIResponse[PredictResponseData](
        status="success",
        data=PredictResponseData(
            smiles_canonical="CCO",
            image_svg="<svg>...</svg>",
            canvas={"width": 300, "height": 300},
            atoms={
                "0": {"x": 50.0, "y": 100.0, "symbol": "C", "smiles": "[CH3]"},
                "1": {"x": 150.0, "y": 100.0, "symbol": "C", "smiles": "[CH2]"},
                "2": {"x": 250.0, "y": 100.0, "symbol": "O", "smiles": "O"}
            },
            bonds={
                "0": {
                    "start": 0,
                    "end": 1,
                    "start_coords": {"x": 50.0, "y": 100.0},
                    "end_coords": {"x": 150.0, "y": 100.0},
                    "bond_atoms": "C-C",
                    "bond_type": "single"
                },
                "1": {
                    "start": 1,
                    "end": 2,
                    "start_coords": {"x": 150.0, "y": 100.0},
                    "end_coords": {"x": 250.0, "y": 100.0},
                    "bond_atoms": "C-O",
                    "bond_type": "single"
                }
            },
            molecule_id="a1b2c3d4e5f6a7b8"
        )
    ).model_dump()

    @extend_schema(
        description="""
        Devuelve información enriquecida de la molécula: SMILES canónico, imagen SVG, metadatos de canvas, posiciones de átomos y enlaces, id único.
        
        Returns enriched information about the molecule: canonical SMILES, SVG image, canvas metadata, atom and bond positions, unique ID.
        """,
        request=PredictRequest,

        responses={
            200: OpenApiResponse(
                response=PredictResponse,
                description="Respuesta exitosa / Successful response",
                examples=[
                    OpenApiExample(
                        name="Ejemplo exitoso",
                        value=example1,
                        response_only=True
                    )
                ]
            ),
            400: OpenApiResponse(
                response=ErrorResponse,
                description="Error de validación / Validation error",
                examples=[
                    OpenApiExample(
                        name="Error por SMILES inválido",
                        value={
                            "status": "error",
                            "detail": "Formato SMILES inválido"
                        },
                        response_only=True
                    )
                ]
            )
        },

        examples=[
            OpenApiExample(
                "Ejemplo de entrada / Example input",
                value={"smiles": "CCO"},
                request_only=True
            ),
            OpenApiExample(
                "Ejemplo alternativo / Alternative example",
                value={"smiles": "C1=CC=CC=C1"},
                request_only=True
            )
        ]
    )
    def post(self, request):
        try:
            req = PredictRequest(**request.data)
            resp_data = predict_controller(req)
            resp = APIResponse[PredictResponseData](status="success", data=resp_data)
            return Response(resp.model_dump())
        except ValidationError as e:
            return _handle_validation_error(e)

class PredictSingleView(APIView):
    """
    Endpoint para predecir la energía de disociación de un enlace específico de la molécula (SMILES y bond_idx).
    """
    example1 = APIResponse[PredictSingleResponseData](
        status="success",
        data=PredictSingleResponseData(
            smiles_canonical="CCO",
            bond={
                "idx": 1,
                "bde": 95.0,
                "begin_atom_idx": 0,
                "end_atom_idx": 1,
                "bond_atoms": "C-O",
                "bond_type": "single"
            }
        )
    ).model_dump()

    @extend_schema(
        description="""
        Predice la energía de disociación para un enlace específico de la molécula.
        
        Predicts the dissociation energy for a specific bond of the molecule.
        """,
        request=PredictSingleRequest,
        responses={
            200: OpenApiResponse(
                response=PredictSingleResponse,
                description="Respuesta exitosa / Successful response",
                examples=[
                    OpenApiExample(
                        name="Ejemplo exitoso",
                        value=example1,
                        response_only=True
                    )
                ]
            ),
            400: OpenApiResponse(
                response=ErrorResponse,
                description="Error de validación / Validation error",
                examples=[
                    OpenApiExample(
                        name="Error por índice de enlace inválido",
                        value={
                            "status": "error",
                            "detail": "Índice de enlace fuera de rango"
                        },
                        response_only=True
                    )
                ]
            )
        },
        examples=[
            OpenApiExample(
                "Ejemplo de entrada / Example input",
                value={"smiles": "CCO", "bond_idx": 1, "molecule_id": "150018eccd174140"},
                request_only=True
            ),
            OpenApiExample(
                "Ejemplo alternativo / Alternative example",
                value={"smiles": "C1=CC=CC=C1", "bond_idx": 2, "molecule_id": "4b7620ed22c55dfd"},
                request_only=True
            )
        ]
    )
    def post(self, request):
        try:
            req = PredictSingleRequest(**request.data)
            resp_data = predict_single_controller(req)
            resp = APIResponse[PredictSingleResponseData](status="success", data=resp_data)
            return Response(resp.model_dump())
        except ValidationError as e:
            return _handle_validation_error(e)

class PredictMultipleView(APIView):
    """
    Endpoint para predecir las energías de disociación para varios enlaces de la molécula (SMILES y bond_indices).
    """
    example1 = APIResponse[PredictMultipleResponseData](
        status="success",
        data=PredictMultipleResponseData(
            smiles="CCO",
            molecule_id="a1b2c3d4e5f6a7b8",
            bonds=[
                {
                    "idx": 1,
                    "bde": 95.0,
                    "begin_atom_idx": 0,
                    "end_atom_idx": 1,
                    "bond_atoms": "C-O",
                    "bond_type": "single"
                },
                {
                    "idx": 2,
                    "bde": 100.0,
                    "begin_atom_idx": 1,
                    "end_atom_idx": 2,
                    "bond_atoms": "C-C",
                    "bond_type": "single"
                }
            ]
        )
    ).model_dump()

    @extend_schema(
        description="""
        Predice las energías de disociación para varios enlaces de la molécula.
        
        Predicts the dissociation energies for multiple bonds of the molecule.
        """,
        request=PredictMultipleRequest,
        responses={
            200: OpenApiResponse(
                response=PredictMultipleResponse,
                description="Respuesta exitosa / Successful response",
                examples=[
                    OpenApiExample(
                        name="Ejemplo exitoso",
                        value=example1,
                        response_only=True
                    )
                ]
            ),
            400: OpenApiResponse(
                response=ErrorResponse,
                description="Error de validación / Validation error",
                examples=[
                    OpenApiExample(
                        name="Error por índices de enlace inválidos",
                        value={
                            "status": "error",
                            "error": {
                                "code": "BOND_INDEX_OUT_OF_RANGE",
                                "message": "Índices de enlace fuera de rango"
                            }
                        },
                        response_only=True
                    )
                ]
            )
        },
        examples=[
            OpenApiExample(
                "Entrada de ejemplo / Example input",
                value={ "smiles": "CCO", "bond_indices": [1, 2], "molecule_id": "150018eccd174140" },
                request_only=True
            ),
            OpenApiExample(
                "Entrada alternativa / Alternative input",
                value={"smiles": "C1=CC=CC=C1", "bond_indices": [0, 2, 3], "molecule_id": "h1g2f3e4d5c6b7a8"},
                request_only=True
            )
        ]
    )
    def post(self, request):
        try:
            logging.info(f"Received request for multiple predictions: {request.data}")
            req = PredictMultipleRequest(**request.data)

            resp_data = predict_multiple_controller(req)
 
            resp = APIResponse[PredictMultipleResponseData](status="success", data=resp_data)
            return Response(resp.model_dump())
        except ValidationError as e:
            return _handle_validation_error(e)

class FragmentView(APIView):
    """
    Unified endpoint to generate molecular fragments in SMILES or XYZ format.
    """
    example1 = APIResponse[FragmentResponseData](
        status="success",
        data=FragmentResponseData(
            smiles_canonical="CCO",
            molecule_id="a1b2c3d4e5f6a7b8",
            bonds=[
                {
                    "idx": 1,
                    "begin_atom": 0,
                    "end_atom": 1,
                    "bond_atoms": "C-O",
                    "is_fragmentable": True,
                    "bond_type": "single"
                }
            ],
            smiles_list=["CC", "O"],
            xyz_block=None
        )
    ).model_dump()

    @extend_schema(
        description="""
        Genera fragmentos moleculares en formato SMILES o XYZ.
        
        Generates molecular fragments in SMILES or XYZ format.
        """,
        request=FragmentRequest,
        responses={
            200: OpenApiResponse(
                response=FragmentResponse,
                description="Respuesta exitosa / Successful response",
                examples=[
                    OpenApiExample(
                        name="Ejemplo exitoso",
                        value=example1,
                        response_only=True
                    )
                ]
            ),
            400: OpenApiResponse(
                response=ErrorResponse,
                description="Error de validación / Validation error",
                examples=[
                    OpenApiExample(
                        name="Error por formato inválido",
                        value={
                            "status": "error",
                            "detail": "Formato no soportado"
                        },
                        response_only=True
                    )
                ]
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
            resp_data = fragment_controller(req)
            resp = APIResponse[FragmentResponseData](status="success", data=resp_data)
            return Response(resp.model_dump())
        except ValidationError as e:
            return _handle_validation_error(e)

class PredictCheckView(APIView):
    """
    Endpoint para verificar si los productos generados por la escisión de un enlace corresponden a los esperados y predecir la BDE.
    """
    example1 = APIResponse[PredictCheckResponseData](
        status="success",
        data=PredictCheckResponseData(
            smiles_canonical="CCO",
            bond={
                "idx": 1,
                "bde": 95.0,
                "begin_atom_idx": 0,
                "end_atom_idx": 1,
                "bond_atoms": "C-O",
                "bond_type": "single"
            },
            products=["CC", "O"]
        )
    ).model_dump()

    @extend_schema(
        description="""
        Verifica productos generados por la escisión de un enlace y predice la BDE.
        
        Verifies products generated by the cleavage of a bond and predicts the BDE.
        """,
        request=PredictCheckRequest,
        responses={
            200: OpenApiResponse(
                response=PredictCheckResponse,
                description="Respuesta exitosa / Successful response",
                examples=[
                    OpenApiExample(
                        name="Ejemplo exitoso",
                        value=example1,
                        response_only=True
                    )
                ]
            ),
            400: OpenApiResponse(
                response=ErrorResponse,
                description="Error de validación / Validation error",
                examples=[
                    OpenApiExample(
                        name="Error por productos inválidos",
                        value={
                            "status": "error",
                            "detail": "Productos no coinciden con los esperados"
                        },
                        response_only=True
                    )
                ]
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
            resp_data = predict_check_controller(req)
            resp = APIResponse[PredictCheckResponseData](status="success", data=resp_data)
            return Response(resp.model_dump())
        except ValidationError as e:
            return _handle_validation_error(e)

class InferAllView(APIView):
    """
    Endpoint para predecir las energías de disociación para todos los enlaces simples de la molécula dada.
    """
    example1 = APIResponse[InferAllResponseData](
        status="success",
        data=InferAllResponseData(
            smiles_canonical="CCO",
            molecule_id="a1b2c3d4e5f6a7b8",
            bonds=[
                {
                    "idx": 1,
                    "bde": 95.0,
                    "begin_atom_idx": 0,
                    "end_atom_idx": 1,
                    "bond_atoms": "C-O",
                    "bond_type": "single"
                },
                {
                    "idx": 2,
                    "bde": 100.0,
                    "begin_atom_idx": 1,
                    "end_atom_idx": 2,
                    "bond_atoms": "C-C",
                    "bond_type": "single"
                }
            ]
        )
    ).model_dump()

    @extend_schema(
        description="""
        Predice las energías de disociación para todos los enlaces simples de la molécula dada.
        
        Predicts the dissociation energies for all single bonds of the given molecule.
        """,
        request=InferAllRequest,
        responses={
            200: OpenApiResponse(
                response=InferAllResponse,
                description="Respuesta exitosa / Successful response",
                examples=[
                    OpenApiExample(
                        name="Ejemplo exitoso",
                        value=example1,
                        response_only=True
                    )
                ]
            ),
            400: OpenApiResponse(
                response=ErrorResponse,
                description="Error de validación / Validation error",
                examples=[
                    OpenApiExample(
                        name="Error por SMILES inválido",
                        value={
                            "status": "error",
                            "detail": "Formato SMILES inválido"
                        },
                        response_only=True
                    )
                ]
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
            resp_data = infer_all_controller(req)
            resp = APIResponse[InferAllResponseData](status="success", data=resp_data)
            return Response(resp.model_dump())
        except ValidationError as e:
            return _handle_validation_error(e)

class DownloadReportView(APIView):
    """
    Endpoint para generar y descargar un informe PDF con los resultados de la predicción para la molécula y enlace indicados.
    """
    example1 = APIResponse[DownloadReportResponseData](
        status="success",
        data=DownloadReportResponseData(
            report_base64="JVBERi0xLjQKJcfs..."  # Ejemplo de string base64
        )
    ).model_dump()

    @extend_schema(
        description="""
        Genera y descarga un informe PDF con los resultados de la predicción para la molécula y enlace indicados.
        
        Generates and downloads a PDF report with the prediction results for the indicated molecule and bond.
        """,
        request=DownloadReportRequest,
        responses={
            200: OpenApiResponse(
                response=DownloadReportResponse,
                description="Respuesta exitosa / Successful response",
                examples=[
                    OpenApiExample(
                        name="Ejemplo exitoso",
                        value=example1,
                        response_only=True
                    )
                ]
            ),
            400: OpenApiResponse(
                response=ErrorResponse,
                description="Error de validación / Validation error",
                examples=[
                    OpenApiExample(
                        name="Error por formato inválido",
                        value={
                            "status": "error",
                            "detail": "Formato no soportado"
                        },
                        response_only=True
                    )
                ]
            )
        },
        examples=[
            OpenApiExample(
                "Entrada de ejemplo / Example input",
                value={"smiles": "CCO", "format": "pdf"},
                request_only=True
            ),
            OpenApiExample(
                "Entrada alternativa / Alternative input",
                value={"smiles": "C1=CC=CC=C1", "format": "pdf"},
                request_only=True
            )
        ]
    )
    def post(self, request):
        try:
            req = DownloadReportRequest(**request.data)  # Validate input
            resp_data = download_report_controller(req)
            resp = APIResponse[DownloadReportResponseData](status="success", data=resp_data)
            return Response(resp.model_dump())
        except ValidationError as e:
            return _handle_validation_error(e)