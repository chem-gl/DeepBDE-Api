from api.model.dto import (
    PredictRequest, PredictResponse, PredictResponseBond,
    PredictSingleRequest, PredictSingleResponse,
    PredictMultipleRequest, PredictMultipleResponse, PredictMultipleBond
)

# Controlador para predicción general, single y multiple

def predict_bde_controller(request: PredictRequest) -> PredictResponse:
    # MOCK
    return PredictResponse(
        image="data:image/png;base64,iVBORw0KGgo...",
        bonds=[PredictResponseBond(idx=1, atoms=[0,1], bde=113.4)]
    )

def predict_single_bde_controller(request: PredictSingleRequest) -> PredictSingleResponse:
    # MOCK
    return PredictSingleResponse(bde=113.7)

def predict_multiple_bde_controller(request: PredictMultipleRequest) -> PredictMultipleResponse:
    # MOCK
    return PredictMultipleResponse(bonds=[
        PredictMultipleBond(idx=1, bde=112.3),
        PredictMultipleBond(idx=2, bde=110.5)
    ])
