from api.model.dto import InferAllRequest, InferAllResponse, InferAllBond

def infer_all_controller(request: InferAllRequest) -> InferAllResponse:
    # MOCK
    return InferAllResponse(bonds=[
        InferAllBond(idx=0, bde=110.2),
        InferAllBond(idx=1, bde=112.3)
    ])
