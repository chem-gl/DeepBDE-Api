from api.model.dto import PredictCheckRequest, PredictCheckResponse

def predict_check_controller(request: PredictCheckRequest) -> PredictCheckResponse:
    return PredictCheckResponse(bde=113.7)
