from rest_framework import serializers

# DTOs para tipado y conversión
from api.model.dto import (
    PredictRequest, PredictSingleRequest, PredictMultipleRequest,
    FragmentSmilesRequest, FragmentXYZRequest, PredictCheckRequest,
    InferAllRequest, DownloadReportRequest
)

class PredictSerializer(serializers.Serializer):
    """Serializer para PredictRequest"""
    smiles = serializers.CharField()

class PredictSingleSerializer(serializers.Serializer):
    """Serializer para PredictSingleRequest"""
    smiles = serializers.CharField()
    bond_idx = serializers.IntegerField()

class PredictMultipleSerializer(serializers.Serializer):
    """Serializer para PredictMultipleRequest"""
    smiles = serializers.CharField()
    bond_indices = serializers.ListField(child=serializers.IntegerField())

class FragmentSmilesSerializer(serializers.Serializer):
    """Serializer para FragmentSmilesRequest"""
    smiles = serializers.CharField()
    all_bonds = serializers.BooleanField(default=False)

class FragmentXYZSerializer(serializers.Serializer):
    """Serializer para FragmentXYZRequest"""
    smiles = serializers.CharField()
    all_bonds = serializers.BooleanField(default=False)

class PredictCheckSerializer(serializers.Serializer):
    """Serializer para PredictCheckRequest"""
    smiles = serializers.CharField()
    bond_idx = serializers.IntegerField()
    products = serializers.ListField(child=serializers.CharField())

class InferAllSerializer(serializers.Serializer):
    """Serializer para InferAllRequest"""
    smiles = serializers.CharField()

class DownloadReportSerializer(serializers.Serializer):
    """Serializer para DownloadReportRequest"""
    smiles = serializers.CharField()
    bond_idx = serializers.IntegerField()
    format = serializers.CharField()
