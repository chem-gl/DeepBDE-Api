from django.urls import path
from .views import (
    MoleculeInfoView, MoleculeSmileCanonicalView, PredictSingleView, PredictMultipleView,
    FragmentView, PredictCheckView, InferAllView, DownloadReportView
)
urlpatterns = [
    path('predict/info/', MoleculeInfoView.as_view(), name='predict'),
    path('predict/info-smile-canonical/', MoleculeSmileCanonicalView.as_view(), name='predict_smile_canonical'),
    path('predict/single/', PredictSingleView.as_view(), name='predict_single'),
    path('predict/multiple/', PredictMultipleView.as_view(), name='predict_multiple'),
    path('fragment/', FragmentView.as_view(), name='fragment'),
    path('predict/check/', PredictCheckView.as_view(), name='predict_check'),
    path('infer/all/', InferAllView.as_view(), name='infer_all'),
    path('download_report/', DownloadReportView.as_view(), name='download_report'),
]
