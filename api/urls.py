from django.urls import path
from .views import *

urlpatterns = [
    path('predict/', PredictView.as_view()),
    path('predict/single/', PredictSingleView.as_view()),
    path('predict/multiple/', PredictMultipleView.as_view()),
    path('fragment/smiles/', FragmentSmilesView.as_view()),
    path('fragment/xyz/', FragmentXYZView.as_view()),
    path('predict/check/', PredictCheckView.as_view()),
    path('infer/all/', InferAllView.as_view()),
    path('info/', InfoView.as_view()),
    path('status/', StatusView.as_view()),
    path('metrics/', MetricsView.as_view()),
    path('download_report/', DownloadReportView.as_view()),
]
