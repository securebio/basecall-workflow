FROM mambaorg/micromamba:1.5.10

COPY automation/environment.yml /tmp/environment.yml

RUN micromamba install -y -n base -f /tmp/environment.yml && \
    micromamba clean --all --yes

ARG MAMBA_DOCKERFILE_ACTIVATE=1

ARG SEQ_IMPORT_SHA=e6ad6266cb1410c595bd93cbb8597cacbb0044e5
# BuildKit mounts secrets as root-owned files under /run/secrets. This build
# step runs as a non-root user in the micromamba image, so relax permissions for
# this single RUN so pip can read the token.
RUN --mount=type=secret,id=gh_token,mode=0444 \
    pip install --no-cache-dir "git+https://x-access-token:$(cat /run/secrets/gh_token)@github.com/securebio/nao-mgs-import.git@${SEQ_IMPORT_SHA}#egg=seq_import"
RUN python -c "import seq_import"

WORKDIR /workflow
COPY main.nf ./
COPY modules/ ./modules/
COPY configs/ ./configs/

ARG BASECALL_WORKFLOW_COMMIT
ENV BASECALL_WORKFLOW_COMMIT=$BASECALL_WORKFLOW_COMMIT

USER mambauser

COPY automation/__init__.py automation/run_automation.py /workflow/automation/
ENTRYPOINT ["/usr/local/bin/_entrypoint.sh", "python", "-m", "automation.run_automation"]
