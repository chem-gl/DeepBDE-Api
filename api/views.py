from rest_framework.views import APIView
from rest_framework.response import Response
from rest_framework import status
from .serializers import (
    PredictSerializer,
    PredictSingleSerializer,
    PredictMultipleSerializer,
    FragmentSmilesSerializer,
    FragmentXYZSerializer,
    PredictCheckSerializer,
    InferAllSerializer,
    DownloadReportSerializer
)

from drf_spectacular.utils import extend_schema, OpenApiExample

# Importar controladores
from api.controllers.predict_controller import (
    predict_bde_controller,
    predict_single_bde_controller,
    predict_multiple_bde_controller
)
from api.controllers.fragment_controller import (
    fragment_smiles_controller,
    fragment_xyz_controller
)
from api.controllers.check_controller import predict_check_controller
from api.controllers.infer_controller import infer_all_controller
from api.controllers.info_controller import (
    info_controller,
    status_controller,
    metrics_controller
)
from api.controllers.report_controller import download_report_controller

# DTOs se importan desde api.model.dto
from api.model.dto import (
    PredictRequest,  
    PredictSingleRequest,  
    PredictMultipleRequest,  
    FragmentSmilesRequest,  
    FragmentXYZRequest,  
    PredictCheckRequest, 
    InferAllRequest,  
    DownloadReportRequest 
)

MOCK_IMAGE = "data:image/png;base64,iVBORw0KGgo..."


class PredictView(APIView):
    """
    Endpoint para obtener la imagen 2D de la molécula y la lista de enlaces con sus índices y átomos involucrados (sin predecir BDEs).
    ---
    Endpoint to get the 2D image of the molecule and the list of bonds with their indices and atoms (no BDE prediction).
    """
    @extend_schema(
        description=(
            "Devuelve información enriquecida de la molécula: SMILES canónico, imagen SVG, metadatos de canvas, "
            "posiciones de átomos y enlaces, id único y SMILES de cada átomo.\n\n"
            "Returns enriched molecule information: canonical SMILES, SVG image, canvas metadata, atom and bond positions, unique id, and atom SMILES."
        ),
        request=PredictSerializer,
        responses={200: None},
        examples=[
            OpenApiExample(
                'Etanol',
                value={
                    "status": "success",
                    "data": {
                        "smiles_canonical": "CCO",
                        "image_svg": "<svg>...</svg>",
                        "canvas": {"width": 300, "height": 300},
                        "atoms": {
                            "0": {"x": 123.45, "y": 67.89},
                            "1": {"x": 156.78, "y": 120.34},
                            "2": {"x": 200.00, "y": 150.00}
                        },
                        "bonds": {
                            "0": {
                                "start": {"x": 123.45, "y": 67.89},
                                "end": {"x": 156.78, "y": 120.34}
                            },
                            "1": {
                                "start": {"x": 156.78, "y": 120.34},
                                "end": {"x": 200.00, "y": 150.00}
                            }
                        },
                        "mol_id": "a1b2c3d4e5f6a7b8",
                        "atom_smiles": ["[C]", "[C]", "[O]"]
                    },
                    "error": None
                },
                response_only=True
            )
        ]
    )
    def post(self, request):
        serializer = PredictSerializer(data=request.data)
        serializer.is_valid(raise_exception=True)
        # Llama al controlador
        dto = PredictRequest(**serializer.validated_data)
        try:
            result = predict_bde_controller(dto)
            return Response({
                "status": "success",
                "data": {
                    "smiles_canonical": result.smiles_canonical,
                    "image_svg": result.image_svg,
                    "canvas": result.canvas,
                    "atoms": result.atoms,
                    "bonds": result.bonds,
                    "mol_id": result.mol_id
                },
                "error": None
            })
        except ValueError as e:
            return Response({
                "status": "error",
                "data": None,
                "error": {
                    "code": "ERR_INVALID_SMILES",
                    "message": str(e)
                }
            }, status=status.HTTP_400_BAD_REQUEST)

class PredictSingleView(APIView):
    """
    Endpoint para predecir la energía de disociación de un enlace específico de la molécula (SMILES y bond_idx).
    ---
    Endpoint to predict the dissociation energy for a specific bond in the molecule (SMILES and bond_idx).
    """
    @extend_schema(
        description="Predice la energía de disociación para un enlace específico de la molécula (SMILES y bond_idx).\n\nPredicts the dissociation energy for a specific bond in the molecule (SMILES and bond_idx).",
        request=PredictSingleSerializer,
        responses={200: PredictSingleSerializer},
        examples=[
            OpenApiExample('Etanol, enlace 1', value={"smiles": "CCO", "bond_idx": 1}, request_only=True),
            OpenApiExample('Benceno, enlace 2', value={"smiles": "c1ccccc1", "bond_idx": 2}, request_only=True),
            OpenApiExample('Metano, enlace 0', value={"smiles": "C", "bond_idx": 0}, request_only=True),
            OpenApiExample('Cafeína, enlace 5', value={"smiles": "Cn1cnc2c1c(=O)n(C)c(=O)n2C", "bond_idx": 5}, request_only=True),
            OpenApiExample('Taxol, enlace 10', value={"smiles": "CC1=C2C(=CC(=O)OC2=CC3=C1C(=O)OC3)OC", "bond_idx": 10}, request_only=True),
            OpenApiExample('Colesterol, enlace 15', value={"smiles": "CC(C)CCCC(C)C1CCC2C3CCC4=CC(=O)CCC4(C)C3CCC12C", "bond_idx": 15}, request_only=True),
            OpenApiExample('Morfina, enlace 8', value={"smiles": "CN1CC[C@]23c4c5ccc(O)c4O[C@H]2[C@@H](O)C=C[C@H]3[C@H]1C5", "bond_idx": 8}, request_only=True),
            OpenApiExample('Punto', value={"smiles": ".", "bond_idx": 1}, request_only=True)
        ]
    )
    def post(self, request):
        serializer = PredictSingleSerializer(data=request.data)
        serializer.is_valid(raise_exception=True)
        dto = PredictSingleRequest(**serializer.validated_data)
        result = predict_single_bde_controller(dto)
        return Response({
            "status": "success",
            "data": {"bde": result.bde},
            "error": None
        })

class PredictMultipleView(APIView):
    """
    Endpoint para predecir las energías de disociación para varios enlaces de la molécula (SMILES y bond_indices).
    ---
    Endpoint to predict dissociation energies for multiple bonds in the molecule (SMILES and bond_indices).
    """
    @extend_schema(
        description="Predice las energías de disociación para varios enlaces de la molécula (SMILES y bond_indices).\n\nPredicts dissociation energies for multiple bonds in the molecule (SMILES and bond_indices).",
        request=PredictMultipleSerializer,
        responses={200: PredictMultipleSerializer},
        examples=[
            OpenApiExample('Etanol, enlaces 1 y 2', value={"smiles": "CCO", "bond_indices": [1,2]}, request_only=True),
            OpenApiExample('Benceno, enlaces 0,1', value={"smiles": "c1ccccc1", "bond_indices": [0,1]}, request_only=True),
            OpenApiExample('Metano, enlace 0', value={"smiles": "C", "bond_indices": [0]}, request_only=True),
            OpenApiExample('Cafeína, enlaces 2,5', value={"smiles": "Cn1cnc2c1c(=O)n(C)c(=O)n2C", "bond_indices": [2,5]}, request_only=True),
            OpenApiExample('Taxol, enlaces 3,10', value={"smiles": "CC1=C2C(=CC(=O)OC2=CC3=C1C(=O)OC3)OC", "bond_indices": [3,10]}, request_only=True),
            OpenApiExample('Colesterol, enlaces 7,15', value={"smiles": "CC(C)CCCC(C)C1CCC2C3CCC4=CC(=O)CCC4(C)C3CCC12C", "bond_indices": [7,15]}, request_only=True),
            OpenApiExample('Morfina, enlaces 4,8', value={"smiles": "CN1CC[C@]23c4c5ccc(O)c4O[C@H]2[C@@H](O)C=C[C@H]3[C@H]1C5", "bond_indices": [4,8]}, request_only=True)
        ]
    )
    def post(self, request):
        serializer = PredictMultipleSerializer(data=request.data)
        serializer.is_valid(raise_exception=True)
        dto = PredictMultipleRequest(**serializer.validated_data)
        result = predict_multiple_bde_controller(dto)
        return Response({
            "status": "success",
            "data": {"bonds": [bond.__dict__ for bond in result.bonds]},
            "error": None
        })

class FragmentSmilesView(APIView):
    """
    Endpoint para generar fragmentos moleculares en formato SMILES a partir de la molécula dada.
    Puede devolver todos los fragmentos por escisión de enlaces simples.
    ---
    Endpoint to generate molecular fragments in SMILES format from the given molecule.
    Can return all fragments by single bond cleavage.
    """
    @extend_schema(
        description="Genera fragmentos moleculares en formato SMILES a partir de la molécula dada. Puede devolver todos los fragmentos por escisión de enlaces simples.\n\nGenerates molecular fragments in SMILES format from the given molecule. Can return all fragments by single bond cleavage.",
        request=FragmentSmilesSerializer,
        responses={200: FragmentSmilesSerializer},
        examples=[
            OpenApiExample('Etanol, todos los enlaces', value={"smiles": "CCO", "all_bonds": True}, request_only=True),
            OpenApiExample('Benceno, solo uno', value={"smiles": "c1ccccc1", "all_bonds": False}, request_only=True),
            OpenApiExample('Metano, todos', value={"smiles": "C", "all_bonds": True}, request_only=True),
            OpenApiExample('Cafeína, todos', value={"smiles": "Cn1cnc2c1c(=O)n(C)c(=O)n2C", "all_bonds": True}, request_only=True),
            OpenApiExample('Taxol, todos', value={"smiles": "CC1=C2C(=CC(=O)OC2=CC3=C1C(=O)OC3)OC", "all_bonds": True}, request_only=True),
            OpenApiExample('Colesterol, todos', value={"smiles": "CC(C)CCCC(C)C1CCC2C3CCC4=CC(=O)CCC4(C)C3CCC12C", "all_bonds": True}, request_only=True),
            OpenApiExample('Morfina, todos', value={"smiles": "CN1CC[C@]23c4c5ccc(O)c4O[C@H]2[C@@H](O)C=C[C@H]3[C@H]1C5", "all_bonds": True}, request_only=True)
        ]
    )
    def post(self, request):
        serializer = FragmentSmilesSerializer(data=request.data)
        serializer.is_valid(raise_exception=True)
        dto = FragmentSmilesRequest(**serializer.validated_data)
        result = fragment_smiles_controller(dto)
        return Response({
            "status": "success",
            "data": {
                "parent": result.parent,
                "fragments": result.fragments,
                "all_fragments": result.all_fragments
            },
            "error": None
        })

class FragmentXYZView(APIView):
    """
    Endpoint para generar fragmentos moleculares en formato XYZ a partir de la molécula dada.
    Puede devolver todos los fragmentos por escisión de enlaces simples.
    ---
    Endpoint to generate molecular fragments in XYZ format from the given molecule.
    Can return all fragments by single bond cleavage.
    """
    @extend_schema(
        description="Genera fragmentos moleculares en formato XYZ a partir de la molécula dada. Puede devolver todos los fragmentos por escisión de enlaces simples.\n\nGenerates molecular fragments in XYZ format from the given molecule. Can return all fragments by single bond cleavage.",
        request=FragmentXYZSerializer,
        responses={200: FragmentXYZSerializer},
        examples=[
            OpenApiExample('Etanol, todos los enlaces', value={"smiles": "CCO", "all_bonds": True}, request_only=True),
            OpenApiExample('Benceno, solo uno', value={"smiles": "c1ccccc1", "all_bonds": False}, request_only=True),
            OpenApiExample('Metano, todos', value={"smiles": "C", "all_bonds": True}, request_only=True),
            OpenApiExample('Cafeína, todos', value={"smiles": "Cn1cnc2c1c(=O)n(C)c(=O)n2C", "all_bonds": True}, request_only=True),
            OpenApiExample('Taxol, todos', value={"smiles": "CC1=C2C(=CC(=O)OC2=CC3=C1C(=O)OC3)OC", "all_bonds": True}, request_only=True),
            OpenApiExample('Colesterol, todos', value={"smiles": "CC(C)CCCC(C)C1CCC2C3CCC4=CC(=O)CCC4(C)C3CCC12C", "all_bonds": True}, request_only=True),
            OpenApiExample('Morfina, todos', value={"smiles": "CN1CC[C@]23c4c5ccc(O)c4O[C@H]2[C@@H](O)C=C[C@H]3[C@H]1C5", "all_bonds": True}, request_only=True)
        ]
    )
    def post(self, request):
        serializer = FragmentXYZSerializer(data=request.data)
        serializer.is_valid(raise_exception=True)
        dto = FragmentXYZRequest(**serializer.validated_data)
        result = fragment_xyz_controller(dto)
        return Response({
            "status": "success",
            "data": result.xyz,
            "error": None
        })

class PredictCheckView(APIView):
    """
    Endpoint para verificar si los productos generados por la escisión de un enlace corresponden a los esperados y predecir la BDE.
    ---
    Endpoint to check if the products generated by bond cleavage match the expected ones and predict the BDE.
    """
    @extend_schema(
        description="Verifica si los productos generados por la escisión de un enlace corresponden a los esperados y predice la BDE.\n\nChecks if the products generated by bond cleavage match the expected ones and predicts the BDE.",
        request=PredictCheckSerializer,
        responses={200: PredictCheckSerializer},
        examples=[
            OpenApiExample('Etanol, productos correctos', value={"smiles": "CCO", "bond_idx": 1, "products": ["[CH3]", "[OH]"]}, request_only=True),
            OpenApiExample('Etanol, productos incorrectos', value={"smiles": "CCO", "bond_idx": 1, "products": ["[C]", "[COH]"]}, request_only=True),
            OpenApiExample('Cafeína', value={"smiles": "Cn1cnc2c1c(=O)n(C)c(=O)n2C", "bond_idx": 5, "products": ["[C]", "[N]"]}, request_only=True),
            OpenApiExample('Taxol', value={"smiles": "CC1=C2C(=CC(=O)OC2=CC3=C1C(=O)OC3)OC", "bond_idx": 10, "products": ["[C]", "[O]"]}, request_only=True)
        ]
    )
    def post(self, request):
        serializer = PredictCheckSerializer(data=request.data)
        serializer.is_valid(raise_exception=True)
        dto = PredictCheckRequest(**serializer.validated_data)
        try:
            result = predict_check_controller(dto)
            return Response({
                "status": "success",
                "data": {"bde": result.bde},
                "error": None
            })
        except Exception as e:
            return Response({
                "status": "error",
                "data": None,
                "error": {
                    "code": "ERR_PRODUCTS_MISMATCH",
                    "message": str(e)
                }
            }, status=status.HTTP_422_UNPROCESSABLE_ENTITY)

class InferAllView(APIView):
    """
    Endpoint para predecir las energías de disociación para todos los enlaces simples de la molécula dada.
    ---
    Endpoint to predict dissociation energies for all single bonds in the given molecule.
    """
    @extend_schema(
        description="Predice las energías de disociación para todos los enlaces simples de la molécula dada.\n\nPredicts dissociation energies for all single bonds in the given molecule.",
        request=InferAllSerializer,
        responses={200: InferAllSerializer},
        examples=[
            OpenApiExample('Etanol', value={"smiles": "CCO"}, request_only=True),
            OpenApiExample('Benceno', value={"smiles": "c1ccccc1"}, request_only=True),
            OpenApiExample('Metano', value={"smiles": "C"}, request_only=True),
            OpenApiExample('Cafeína', value={"smiles": "Cn1cnc2c1c(=O)n(C)c(=O)n2C"}, request_only=True),
            OpenApiExample('Taxol', value={"smiles": "CC1=C2C(=CC(=O)OC2=CC3=C1C(=O)OC3)OC"}, request_only=True),
            OpenApiExample('Colesterol', value={"smiles": "CC(C)CCCC(C)C1CCC2C3CCC4=CC(=O)CCC4(C)C3CCC12C"}, request_only=True),
            OpenApiExample('Morfina', value={"smiles": "CN1CC[C@]23c4c5ccc(O)c4O[C@H]2[C@@H](O)C=C[C@H]3[C@H]1C5"}, request_only=True)
        ]
    )
    def post(self, request):
        serializer = InferAllSerializer(data=request.data)
        serializer.is_valid(raise_exception=True)
        dto = InferAllRequest(**serializer.validated_data)
        result = infer_all_controller(dto)
        return Response({
            "status": "success",
            "data": {"bonds": [bond.__dict__ for bond in result.bonds]},
            "error": None
        })

class InfoView(APIView):
    """
    Endpoint para obtener información sobre la versión del modelo y metadatos relevantes.
    ---
    Endpoint to get information about the model version and relevant metadata.
    """
    @extend_schema(
        description="Devuelve información sobre la versión del modelo y metadatos relevantes.\n\nReturns information about the model version and relevant metadata.",
        responses={200: None},
        examples=[
            OpenApiExample('Punto', value={".": None}, request_only=True)
        ]
    )
    def get(self, request):
        data = info_controller()
        return Response({
            "status": "success",
            "data": data,
            "error": None
        })

class StatusView(APIView):
    """
    Endpoint para obtener el estado actual del sistema y recursos disponibles.
    ---
    Endpoint to get the current system status and available resources.
    """
    @extend_schema(
        description="Devuelve el estado actual del sistema y recursos disponibles.\n\nReturns the current system status and available resources.",
        responses={200: None},
        examples=[
            OpenApiExample('Punto', value={".": None}, request_only=True)
        ]
    )
    def get(self, request):
        data = status_controller()
        return Response({
            "status": "success",
            "data": data,
            "error": None
        })

class MetricsView(APIView):
    """
    Endpoint para obtener métricas Prometheus para monitoreo del sistema.
    ---
    Endpoint to get Prometheus metrics for system monitoring.
    """
    @extend_schema(
        description="Devuelve métricas Prometheus para monitoreo del sistema.\n\nReturns Prometheus metrics for system monitoring.",
        responses={200: None},
        examples=[
            OpenApiExample('Punto', value={".": None}, request_only=True)
        ]
    )
    def get(self, request):
        metrics = metrics_controller()
        return Response(metrics, content_type="text/plain")

class DownloadReportView(APIView):
    """
    Endpoint para generar y descargar un informe PDF con los resultados de la predicción para la molécula y enlace indicados.
    ---
    Endpoint to generate and download a PDF report with the prediction results for the specified molecule and bond.
    """
    @extend_schema(
        description="Genera y descarga un informe PDF con los resultados de la predicción para la molécula y enlace indicados.\n\nGenerates and downloads a PDF report with the prediction results for the specified molecule and bond.",
        request=DownloadReportSerializer,
        responses={200: DownloadReportSerializer},
        examples=[
            OpenApiExample('Etanol', value={"smiles": "CCO", "bond_idx": 1, "format": "pdf"}, request_only=True),
            OpenApiExample('Benceno', value={"smiles": "c1ccccc1", "bond_idx": 2, "format": "pdf"}, request_only=True),
            OpenApiExample('Metano', value={"smiles": "C", "bond_idx": 0, "format": "pdf"}, request_only=True),
            OpenApiExample('Cafeína', value={"smiles": "Cn1cnc2c1c(=O)n(C)c(=O)n2C", "bond_idx": 5, "format": "pdf"}, request_only=True),
            OpenApiExample('Taxol', value={"smiles": "CC1=C2C(=CC(=O)OC2=CC3=C1C(=O)OC3)OC", "bond_idx": 10, "format": "pdf"}, request_only=True),
            OpenApiExample('Colesterol', value={"smiles": "CC(C)CCCC(C)C1CCC2C3CCC4=CC(=O)CCC4(C)C3CCC12C", "bond_idx": 15, "format": "pdf"}, request_only=True),
            OpenApiExample('Morfina', value={"smiles": "CN1CC[C@]23c4c5ccc(O)c4O[C@H]2[C@@H](O)C=C[C@H]3[C@H]1C5", "bond_idx": 8, "format": "pdf"}, request_only=True)
        ]
    )
    def post(self, request):
        serializer = DownloadReportSerializer(data=request.data)
        serializer.is_valid(raise_exception=True)
        dto = DownloadReportRequest(**serializer.validated_data)
        try:
            result = download_report_controller(dto)
            return Response({
                "status": "success",
                "data": result.data,
                "error": None
            })
        except Exception as e:
            return Response({
                "status": "error",
                "data": None,
                "error": {
                    "code": "ERR_REPORT_FAILED",
                    "message": str(e)
                }
            }, status=status.HTTP_400_BAD_REQUEST)