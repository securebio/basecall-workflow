# Scoping: GitHub auth for `seq_import` install at docker build time

## Problem statement

PR 4 of the ONT basecalling automation plan needs to `pip install seq_import`
from the private repo `securebio/nao-mgs-import` at `docker build` time, so the
head container's entrypoint can run
`python -m seq_import samplesheet --delivery <delivery>` after `nextflow run`
finishes. We want to pin by git SHA. Doing so requires GitHub auth inside the
docker build, since the source repo is private.

## Constraint: build-time install, not runtime clone

We want the image SHA to pin the `seq_import` version, not have it depend on
runtime config. Build-time install gives us reproducible images; runtime clone
would let a single image float across `seq_import` revisions, which is the
opposite of what we want. The auth options below are therefore all build-time.

(See "Note on runtime-clone" at the bottom — `nao-mgs-partner-reports` does
runtime clone for a different reason; we explicitly considered and rejected
that pattern here.)

## Options to evaluate

### Option 1 — fine-grained PAT

Mechanism. A user generates a fine-grained personal access token scoped to
`securebio/nao-mgs-import` with `Contents: Read`, stored as a GitHub Actions
secret on `basecall-workflow`. The ecr-push workflow forwards it to docker
build via BuildKit `--mount=type=secret`, and the Dockerfile installs with
`pip install git+https://x-access-token:$TOKEN@github.com/securebio/nao-mgs-import.git@<sha>`.

Pros.

- Conceptually simple; minimal upfront infra setup.
- Precedent in this org: `nao-mgs-partner-reports` uses `METADATA_REPO_TOKEN`
  for `nao-mgs-metadata` (though at runtime, not build-time).

Cons.

- Tied to an individual user. If the owner leaves or rotates accounts, the
  token must be regenerated.
- Fine-grained PATs have a mandatory expiry (max 1 year), so this is a
  recurring rotation chore.
- Audit trail attributes pulls to the PAT owner, not to a service identity.
- Org admin approval may be required for each regeneration depending on org
  policy.

Operational notes. If we go this route, document the owner + expiry date in
this repo's README so the rotation doesn't surprise us.

### Option 2 — deploy key

Mechanism. SSH keypair: public half registered as a deploy key on
`securebio/nao-mgs-import`, private half stored as a GitHub Actions secret on
`basecall-workflow`. The ecr-push workflow forwards the private key to docker
build via BuildKit `--mount=type=ssh` (or `--mount=type=secret`), and the
Dockerfile installs with
`pip install git+ssh://git@github.com/securebio/nao-mgs-import.git@<sha>`.

Two sub-options for which key to use:

- **2a — reuse existing `mgs-import-github-deploy`**. This deploy key already
  exists, currently in AWS Secrets Manager, used by the legacy
  `startSeqPipelines` lambda (see
  `../mgs-import/_infrastructure/infrastructure_setup.md`). Cheapest path, but
  creates short-lived dual storage of the same private key (Secrets Manager +
  GHA secrets) until the legacy consumer is decommissioned in PR 6 of the
  plan.
- **2b — new key dedicated to basecall-workflow**. Generate a fresh keypair
  and register it as an additional deploy key on `nao-mgs-import`. Clean
  separation; can rotate independently from the legacy key.

Mechanically, a deploy key is registered once on the source repo; the private
half is just a secret that any number of consumers can hold. Hygiene best
practice is one private key per consumer (so a leak only forces rotation in
that one place), not a GitHub-imposed constraint.

Pros.

- No expiry; survives user turnover.
- Tied to the source repo, not a person.
- Pattern already in use elsewhere in the org (`mgs-import-github-deploy`).

Cons.

- SSH plumbing in docker build (ssh-agent, `known_hosts` pinning for
  `github.com`).
- One-key-per-consumer hygiene means proliferation if many private deps grow
  to many consumers.

### Option 3 — GitHub App

Mechanism. Register a securebio-org GitHub App (e.g. `nao-ci-reader`) with
`Contents: Read` permission, install it on `securebio/nao-mgs-import` (and any
future private deps). The ecr-push workflow uses
`actions/create-github-app-token@v1` to mint a ~1h installation token, then
forwards it to docker build via `--mount=type=secret`. The Dockerfile
installs with
`pip install git+https://x-access-token:$TOKEN@github.com/securebio/nao-mgs-import.git@<sha>`.

Pros.

- Org-level identity (not a user).
- One app grants access to many consumer repos × many source repos: this
  pattern scales if more services need to read more private deps.
- Short-lived tokens (~1h) reduce blast radius of a leak.
- Centralized revocation; uninstall once and all consumers lose access.

Cons.

- Most upfront setup (~30–60 min): an org admin has to register the app,
  generate the private key, and install it on each source repo.
- One extra workflow step (`create-github-app-token`) and a private-key
  secret on each consumer.

## Note on runtime-clone (explicitly out of scope)

`nao-mgs-partner-reports` clones `nao-mgs-metadata` at container startup,
fetching a token from AWS Secrets Manager. We considered mirroring that
pattern here and rejected it: we want the image SHA to pin the `seq_import`
version, and runtime clone would let a single image float across
`seq_import` revisions. Documented here so reviewers don't re-litigate.

## Recommendation

_Leave blank — to be filled in after team discussion._
