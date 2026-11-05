#!/bin/bash
# Build the image locally with `docker build -t peatguard-pipeline:local .`
# or point PEATGUARD_IMAGE at a registry image.
IMAGE="${PEATGUARD_IMAGE:-peatguard-pipeline:local}"
DIR="$(cd "$(dirname "$0")" && pwd)"

DOCKER_ARGS=(
  --platform linux/amd64
  -v "$DIR/output:/workspace/output"
  -v "$DIR/scratch:/workspace/scratch"
  -v "$HOME/.netrc:/root/.netrc:ro"
  -e "PEATGUARD_STORAGE__OUTPUT_DIR=/workspace/output"
  -e "PEATGUARD_STORAGE__SCRATCH_DIR=/workspace/scratch"
  --memory=32g
  --cpus=8
)

case "$1" in
  ingest)     docker run --rm "${DOCKER_ARGS[@]}" "$IMAGE" download --start 2023-01-01 --end 2024-12-31 --level SLC ;;
  insar)      docker run --rm "${DOCKER_ARGS[@]}" "$IMAGE" process --mode insar --workers 1 ;;
  timeseries) docker run --rm "${DOCKER_ARGS[@]}" "$IMAGE" timeseries ;;
  backscatter)docker run --rm "${DOCKER_ARGS[@]}" "$IMAGE" process --mode backscatter ;;
  analyze)    docker run --rm "${DOCKER_ARGS[@]}" "$IMAGE" analyze ;;
  all)        $0 ingest && $0 insar && $0 timeseries && $0 backscatter && $0 analyze ;;
  shell)      docker run --rm -it "${DOCKER_ARGS[@]}" --entrypoint /bin/bash "$IMAGE" ;;
  *)          echo "Usage: $0 {ingest|insar|timeseries|backscatter|analyze|all|shell}" ;;
esac
