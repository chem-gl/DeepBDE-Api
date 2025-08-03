# DeepBDE API

## Disk Cache System

DeepBDE API uses a disk-based cache to speed up repeated queries for molecules, BDE predictions, and fragmentations.

- **Location:** Cache files are stored in `api/controllers/_cache/v1/`.
- **Clearing the cache:** Delete the `_cache` folder or change the cache version in the code (`CACHE_VERSION` in `predict_controller.py`).
- **Automatic regeneration:** The cache is rebuilt as new queries are made.
- **Purpose:** This system improves performance for repeated or similar requests and is safe to delete at any time.

## Main Features

- BDE prediction for single or multiple bonds in molecules from SMILES.
- 2D molecule image generation and bond visualization.
- Molecular fragmentation in SMILES and XYZ formats, with BDE values for each fragment.
- Validation of cleavage products and BDE calculation for proposed products.
- Downloadable TXT or PDF reports with prediction and fragmentation results.
- Endpoints for model info, system status, and Prometheus metrics.
- Interactive documentation examples for simple and complex molecules (ethanol, benzene, caffeine, taxol, cholesterol, morphine, etc.).
- OpenAPI YAML schema always up-to-date and available at `/api/public/openapi.yaml`.

## API Endpoints

- `POST /api/v1/predict/` — Get enriched molecule info: canonical SMILES, SVG image, atom and bond positions, and unique ID.
- `POST /api/v1/predict/single/` — Predict the BDE for a specific bond in a molecule.
- `POST /api/v1/predict/multiple/` — Predict BDEs for multiple bonds in a molecule.
- `POST /api/v1/fragment/` — Generate molecular fragments in SMILES and/or XYZ format, with BDE values for each fragment.
- `POST /api/v1/predict/check/` — Validate cleavage products and compare with predicted fragments, including BDE calculation.
- `POST /api/v1/infer/all/` — Predict BDEs for all single bonds in a molecule.
- `POST /api/v1/download_report/` — Download a TXT (or PDF) report with prediction and fragmentation results.
- `POST /api/v1/predict/info-smile-canonical/` — Get the canonical SMILES and unique ID for a molecule.

## Usage Flow

1. Send a SMILES to `/api/v1/predict/` to get the 2D image and bond indices.
2. Use `/api/v1/predict/single/` to predict the BDE of a specific bond.
3. Use `/api/v1/predict/multiple/` to predict the BDEs of multiple bonds.
4. Generate molecular fragments (SMILES and/or XYZ) with `/api/v1/fragment/` (set `export_smiles` and/or `export_xyz`).
5. Validate cleavage products and compare with predicted fragments using `/api/v1/predict/check/`.
6. Predict all single bond BDEs with `/api/v1/infer/all/`.
7. Download a TXT or PDF report with `/api/v1/download_report/` (set `format` to `txt` or `pdf`).
8. Get canonical SMILES and molecule ID with `/api/v1/predict/info-smile-canonical/`.

## Installation & Deployment

1. Clone the repository:

   ```bash
   # NOTE: Replace with your actual repo URL if different
   git clone https://github.com/chem-gl/DeepBDE-Api.git
   cd DeepBDE-Api
   ```

2. Build the Docker image:

   ```bash
   docker build -t deepbde-api .
   ```

3. Start the container:

   ```bash
   docker run -p 8000:8000 deepbde-api
   ```

4. The API will be available at `http://localhost:8000/`.

## Interactive Documentation

- Swagger UI: `/api/schema/swagger-ui/`
- Redoc: `/api/schema/redoc/`
- OpenAPI YAML schema: `/api/public/openapi.yaml`

## Project Structure

- `api/` — Business logic, views, serializers, DTOs, and controllers.
- `api/model/dto.py` — Dataclass definitions for requests and responses.
- `api/controllers/` — Controllers for each endpoint group.
- `api/controllers/_cache/` — Disk cache for molecule queries and predictions.
- `api/public/openapi.yaml` — Automatically generated OpenAPI schema.
- `predictor/` — Main Django configuration and routes.
- `deepbde/` — DeepBDE model code and utilities.

## Usage Examples

The Swagger documentation includes examples for both simple and complex molecules (ethanol, benzene, caffeine, taxol, cholesterol, morphine, etc.) in all endpoints.

## License

This project is for academic and research use. See the LICENSE file for more details.

---

Developed by the DeepBDE team. For questions or support, contact the repository authors.
