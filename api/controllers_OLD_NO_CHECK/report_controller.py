from api.model.dto import DownloadReportRequest, DownloadReportResponseData

def download_report_controller(request: DownloadReportRequest) -> DownloadReportResponseData:
    # Mock: returns a base64 string for PDF
    return DownloadReportResponseData(data="JVBERi0xLjQKJcfs...mockbase64pdf...")


