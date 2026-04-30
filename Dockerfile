FROM mambaorg/micromamba:1.5.10

COPY environment.yml /tmp/environment.yml

RUN micromamba install -y -n base -f /tmp/environment.yml && \
    micromamba clean --all --yes

ARG MAMBA_DOCKERFILE_ACTIVATE=1

WORKDIR /workflow
COPY main.nf ./
COPY modules/ ./modules/
COPY configs/ ./configs/

ARG BASECALL_WORKFLOW_COMMIT
ENV BASECALL_WORKFLOW_COMMIT=$BASECALL_WORKFLOW_COMMIT

USER mambauser
