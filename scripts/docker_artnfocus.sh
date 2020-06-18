#!/bin/sh

# This script wraps the docker incantation required to map the current directory into the
# container's WORKDIR and run the focus script using any arguments passed to this script.
# Note that filenames passed to this script must be relative to and beneath the directory
# where this script is run.

docker run --rm -ti -v ${PWD}:/artn tepickering/artnfocus:latest $*
