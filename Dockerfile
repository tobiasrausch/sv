# Build: docker build -t trausch/sv:latest .
# Run: docker run -it -p 8888:8888 trausch/sv:latest
# Run & mount the data: docker run -it -p 8888:8888 -v /data/lr:/opt/sv/data/lr trausch/sv:latest
# Cloud: ssh -L 8888:localhost:8888 user@host
# Data download (inside the container): FILE=<id> pixi run download
FROM ubuntu:24.04

LABEL maintainer="Tobias Rausch <rausch@embl.de>"
LABEL org.opencontainers.image.source="https://github.com/tobiasrausch/sv"

RUN apt-get update && apt-get install -y --no-install-recommends curl ca-certificates git bzip2 procps && apt-get clean && rm -rf /var/lib/apt/lists/*

# pixi
RUN curl -fsSL https://pixi.sh/install.sh | bash
ENV PATH=/root/.pixi/bin:${PATH}

WORKDIR /opt/sv

COPY pixi.toml pixi.lock ./
RUN pixi install --locked

# pixi env on PATH
ENV PATH=/opt/sv/.pixi/envs/default/bin:${PATH}

COPY . .

EXPOSE 8888

# launch the browser workbench
CMD ["pixi", "run", "lab"]
