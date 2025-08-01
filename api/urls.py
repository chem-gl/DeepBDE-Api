from django.urls import path
from .views import (
    PredictView, PredictSingleView, PredictMultipleView,
    FragmentView, PredictCheckView, InferAllView, DownloadReportView
)

urlpatterns = [
    path('predict/', PredictView.as_view(), name='predict'),
    path('predict/single/', PredictSingleView.as_view(), name='predict_single'),
    path('predict/multiple/', PredictMultipleView.as_view(), name='predict_multiple'),
    path('fragment/', FragmentView.as_view(), name='fragment'),
    path('predict/check/', PredictCheckView.as_view(), name='predict_check'),
    path('infer/all/', InferAllView.as_view(), name='infer_all'),
    path('download_report/', DownloadReportView.as_view(), name='download_report'),
]
