from api.model.dto import DownloadReportRequest, DownloadReportResponse

def download_report_controller(request: DownloadReportRequest) -> DownloadReportResponse:
    return DownloadReportResponse(data="PDF generado (mock)")


