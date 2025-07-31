FROM pytorch/pytorch:2.7.1-cuda11.8-cudnn9-runtime
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
        build-essential \
        git \
        latexmk && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*
WORKDIR /app
COPY requirements.txt /app/
COPY DeepBDE-Fork/requirements.txt /app/DeepBDE-Fork/
RUN pip install --upgrade pip && \
    pip install -r DeepBDE-Fork/requirements.txt && \
    pip cache purge && \
    pip install -r requirements.txt

COPY . /app
COPY entrypoint.sh /entrypoint.sh
RUN chmod +x /entrypoint.sh
EXPOSE 8000
ENTRYPOINT ["/entrypoint.sh"]
