import requests
import base64

# URL del endpoint
url = "http://localhost:8000/api/v1/download_report/"

# Datos de la solicitud
data = {
    "smiles": "CCO",
    "format": "txt"
}

# Encabezados
headers = {
    "accept": "application/json",
    "Content-Type": "application/json",
    "X-CSRFTOKEN": "i5HTbau3vLgHJ97kRWnLsWaeqIXAHSAJGd2flajTzQoJkoUrjpQ6isqc8hTCz5Uw"
}

# Realiza la solicitud POST
response = requests.post(url, json=data, headers=headers)

# Verifica si la solicitud fue exitosa
if response.status_code == 200:
    response_json = response.json()
    if response_json.get("status") == "success":
        # Decodifica el contenido Base64 y guarda el PDF
        pdf_base64 = response_json["data"]["report_base64"]
        try:
            # Limpia caracteres no válidos
            pdf_base64_cleaned = "".join(filter(lambda x: x.isascii(), pdf_base64))
            with open("reporte.pdf", "wb") as pdf_file:
                pdf_file.write(base64.b64decode(pdf_base64_cleaned))
            print("El archivo PDF se ha guardado como 'reporte.pdf'.")
        except Exception as e:
            print(f"Error al decodificar el PDF: {e}")
    else:
        print("Error en la respuesta:", response_json.get("error"))
else:
    print(f"Error en la solicitud: {response.status_code}")