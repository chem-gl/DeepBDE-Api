from typing import Dict, Any

def info_controller() -> Dict[str, Any]:
    # MOCK
    return {
        "model_version": "v1.0.2",
        "arxiv": "https://arxiv.org/abs/2306.12345",
        "loaded_at": "2025-07-30T13:00:00Z"
    }

def status_controller() -> Dict[str, Any]:
    # MOCK
    return {
        "deepbde_submodule": "initialized",
        "gpu_available": True,
        "uptime": "2h34m"
    }

def metrics_controller() -> str:
    # MOCK
    return "# HELP deepbde_requests_total Total requests\n# TYPE deepbde_requests_total counter\ndeepbde_requests_total 42"
