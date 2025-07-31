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

# DTOs se importan desde api.model.dto
from api.model.dto import *

MOCK_IMAGE = "data:image/png;base64,iVBORw0KGgo..."


class PredictView(APIView):
    """
    Endpoint para predecir las energías de disociación de enlaces (BDE) y mostrar la estructura 2D de la molécula dada (SMILES).
    Devuelve una imagen y la lista de enlaces con sus BDEs.
    ---
    Endpoint to predict bond dissociation energies (BDE) and display the 2D structure for the given molecule (SMILES).
    Returns an image and the list of bonds with their BDEs.
    """
    @extend_schema(
        description="Predice las energías de disociación de enlaces (BDE) y muestra la estructura 2D de la molécula dada (SMILES).\n\nPredicts bond dissociation energies (BDE) and displays the 2D structure for the given molecule (SMILES).",
        request=PredictSerializer,
        responses={200: PredictSerializer},
        examples=[
            OpenApiExample('Etanol', value={"smiles": "CCO"}, request_only=True),
            OpenApiExample('Benceno', value={"smiles": "c1ccccc1"}, request_only=True),
            OpenApiExample('Metano', value={"smiles": "C"}, request_only=True),
            OpenApiExample('Cafeína', value={"smiles": "Cn1cnc2c1c(=O)n(C)c(=O)n2C"}, request_only=True),
            OpenApiExample('Taxol', value={"smiles": "CC1=C2C(=CC(=O)OC2=CC3=C1C(=O)OC3)OC"}, request_only=True),
            OpenApiExample('Colesterol', value={"smiles": "CC(C)CCCC(C)C1CCC2C3CCC4=CC(=O)CCC4(C)C3CCC12C"}, request_only=True),
            OpenApiExample('Morfina', value={"smiles": "CN1CC[C@]23c4c5ccc(O)c4O[C@H]2[C@@H](O)C=C[C@H]3[C@H]1C5"}, request_only=True),
            OpenApiExample('Punto', value={"smiles": "."}, request_only=True)
        ]
    )
    def post(self, request):
        serializer = PredictSerializer(data=request.data)
        serializer.is_valid(raise_exception=True)
        return Response({
            "status": "success",
            "data": {
                "image": MOCK_IMAGE,
                "bonds": [{"idx": 1, "atoms": [0,1], "bde": 113.4}]
            },
            "error": None
        })

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
        return Response({
            "status": "success",
            "data": {"bde": 113.7},
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
        return Response({
            "status": "success",
            "data": {"bonds": [
                {"idx": 1, "bde": 112.3},
                {"idx": 2, "bde": 110.5}
            ]},
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
        return Response({
            "status": "success",
            "data": {
                "parent": "CCO",
                "fragments": ["[CH3]", "[OH]"],
                "all_fragments": [
                    {"bond_idx": 0, "fragments": ["CC", "O"]},
                    {"bond_idx": 1, "fragments": ["C", "CO"]}
                ]
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
        xyz = "3\n\nC 0.000 0.000 0.000\nC 1.200 0.000 0.000\nO 2.400 0.000 0.000\n\n1\n\nO 0.000 0.000 0.000"
        return Response({
            "status": "success",
            "data": xyz,
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
        products = serializer.validated_data.get('products', [])
        if products != ["[CH3]", "[OH]"]:
            return Response({
                "status": "error",
                "data": None,
                "error": {
                    "code": "ERR_PRODUCTS_MISMATCH",
                    "message": "Los productos no corresponden al corte indicado."
                }
            }, status=status.HTTP_422_UNPROCESSABLE_ENTITY)
        return Response({
            "status": "success",
            "data": {"bde": 113.7},
            "error": None
        })

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
        return Response({
            "status": "success",
            "data": {"bonds": [
                {"idx": 0, "bde": 110.2},
                {"idx": 1, "bde": 112.3}
            ]},
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
        return Response({
            "status": "success",
            "data": {
                "model_version": "v1.0.2",
                "arxiv": "https://arxiv.org/abs/2306.12345",
                "loaded_at": "2025-07-30T13:00:00Z"
            },
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
        return Response({
            "status": "success",
            "data": {
                "deepbde_submodule": "initialized",
                "gpu_available": True,
                "uptime": "2h34m"
            },
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
        # Mock Prometheus metrics
        return Response("# HELP deepbde_requests_total Total requests\n# TYPE deepbde_requests_total counter\ndeepbde_requests_total 42", content_type="text/plain")

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
        fmt = serializer.validated_data.get('format', '').lower()
        if fmt != 'pdf':
            return Response({
                "status": "error",
                "data": None,
                "error": {
                    "code": "ERR_REPORT_FAILED",
                    "message": "Solo se soporta formato PDF."
                }
            }, status=status.HTTP_400_BAD_REQUEST)
        return Response({
            "status": "success",
            "data": "PDF generado (mock)",
            "error": None
        })
