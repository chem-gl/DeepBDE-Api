# DeepBDE API

**DeepBDE API** is a RESTful platform for bond dissociation energy (BDE) prediction and molecular analysis, designed for both human users and frontend applications. The API is fully bilingual (Spanish and English) and documented with OpenAPI/Swagger.

## Main Features

- BDE prediction for single or multiple bonds in molecules from SMILES.
- 2D molecule image generation and bond visualization.
- Molecular fragmentation in SMILES and XYZ formats, with the option to get all possible fragments.
- Validation of cleavage products and BDE calculation for proposed products.
- Downloadable PDF reports with prediction results.
- Endpoints for model info, system status, and Prometheus metrics.
- Interactive documentation examples for simple and complex molecules (ethanol, benzene, caffeine, taxol, cholesterol, morphine, etc.).
- OpenAPI YAML schema always up-to-date and available at `/api/public/openapi.yaml`.

## API Usage Flow

1. Send a SMILES to `/api/v1/predict/` to get the 2D image and bonds with their indices.
2. To predict the BDE of a specific bond, use `/api/v1/predict/single/` with the SMILES and bond index.
3. For multiple bonds, use `/api/v1/predict/multiple/` with the SMILES and a list of indices.
4. To get SMILES fragments, use `/api/v1/fragment/smiles/` (enable `all_bonds` for all fragments).
5. To get fragments in XYZ format, use `/api/v1/fragment/xyz/` (with `all_bonds` if needed).
6. To validate cleavage products and get the BDE, use `/api/v1/predict/check/`.
7. To predict all single bond BDEs, use `/api/v1/infer/all/`.
8. To download a PDF report, use `/api/v1/download_report/`.
9. For model info, use `/api/v1/info/`.
10. For system status, use `/api/v1/status/`.
11. For Prometheus metrics, use `/api/v1/metrics/`.

## Installation & Deployment

1. Clone the repository with all submodules and navigate to the main directory:

   ```bash
   git clone --recurse-submodules https://github.com/chem-gl/DeepBDE-Api.git
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

- `api/` - Business logic, views, serializers, DTOs, and controllers.
- `api/model/dto.py` - Dataclass definitions for requests and responses.
- `api/controllers/` - Controllers for each endpoint group.
- `api/public/openapi.yaml` - Automatically generated OpenAPI schema.
- `predictor/` - Main Django configuration and routes.
- `deepbde/` - DeepBDE model code and utilities.

## Usage Examples

The Swagger documentation includes examples for both simple and complex molecules (ethanol, benzene, caffeine, taxol, cholesterol, morphine, etc.) in all endpoints.

## License

This project is for academic and research use. See the LICENSE file for more details.

---

Developed by the DeepBDE team. For questions or support, contact the repository authors.
