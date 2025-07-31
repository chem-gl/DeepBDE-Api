from api.model.dto import (
    PredictRequest, PredictResponse, PredictResponseBond,
    PredictSingleRequest, PredictSingleResponse,
    PredictMultipleRequest, PredictMultipleResponse, PredictMultipleBond
)

# Controlador para predicción general, single y multiple

def predict_bde_controller(request: PredictRequest) -> PredictResponse:
    """
    Devuelve la imagen 2D de la molécula y la lista de enlaces con sus índices y átomos involucrados.
    """
    # Aquí deberías generar la imagen 2D y extraer los enlaces usando RDKit.
    # MOCK: solo estructura de respuesta, sin BDE
    return PredictResponse(
        image="data:image/png;base64,iVBORw0KGgo...",
        bonds=[
            PredictResponseBond(idx=0, atoms=[0,1], bde=None),
            PredictResponseBond(idx=1, atoms=[1,2], bde=None)
        ]
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
