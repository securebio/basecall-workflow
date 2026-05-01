FROM mambaorg/micromamba:1.5.10

COPY environment.yml /tmp/environment.yml

RUN micromamba install -y -n base -f /tmp/environment.yml && \
    micromamba clean --all --yes

ARG MAMBA_DOCKERFILE_ACTIVATE=1

ARG SEQ_IMPORT_SHA=e6ad6266cb1410c595bd93cbb8597cacbb0044e5
RUN --mount=type=secret,id=gh_token \
    pip install "git+https://x-access-token:$(cat /run/secrets/gh_token)@github.com/securebio/nao-mgs-import.git@${SEQ_IMPORT_SHA}#egg=seq_import"
RUN python -c "import seq_import; print(seq_import.__version__)"

WORKDIR /workflow
COPY main.nf ./
COPY modules/ ./modules/
COPY configs/ ./configs/

ARG BASECALL_WORKFLOW_COMMIT
ENV BASECALL_WORKFLOW_COMMIT=$BASECALL_WORKFLOW_COMMIT

USER mambauser
