import base64
from api.controllers.predict_controller import report_pdf

def test_report_pdf():
    # SMILES de prueba
    smiles = "CCO"  # Etanol

    # Generar el reporte en formato PDF
    pdf_base64 = report_pdf(smiles)

    # Decodificar el contenido Base64 y guardar el PDF
    pdf_content = base64.b64decode(pdf_base64.encode("latin1"))
    with open("test_report.pdf", "wb") as pdf_file:
        pdf_file.write(pdf_content)

    print("PDF generado y guardado como 'test_report.pdf'.")

if __name__ == "__main__":
    test_report_pdf()
