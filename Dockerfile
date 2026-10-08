# Environment only: openEMS and its Python venv.  The CCF code is not in the image: it comes
# from the git clone and is mounted at run time, so editing it never needs a rebuild.
# The tag is the openEMS-Project commit.  Build without a context, so runs/ never reaches Docker:
#   docker build -t ghcr.io/earthyyl/loopgap:9f5cdd4 - < Dockerfile

# Same OS, glibc and Python (3.12) as the WSL build the local results came from
FROM ubuntu:24.04
ENV DEBIAN_FRONTEND=noninteractive
RUN apt-get update && apt-get install -y --no-install-recommends git ca-certificates

# The exact openEMS version the local results came from
RUN git clone https://github.com/thliebig/openEMS-Project.git /src \
 && cd /src && git checkout 9f5cdd4 && git submodule update --init --recursive

# Upstream's own dependency installer; as root it runs apt directly
RUN cd /src && ./scripts/install_deps.sh --auto --python --disable-gui

# Same command as the local build.  It creates /opt/openEMS/venv itself and pip installs
# CSXCAD, openEMS, numpy, matplotlib and h5py into it.
RUN cd /src && ./update_openEMS.sh /opt/openEMS --python --disable-GUI

# What was added to the local venv after the build (cma), plus the build's own picks pinned to
# the local venv's versions, in case PyPI has moved on since 2026-09-27
RUN /opt/openEMS/venv/bin/pip install --no-cache-dir cma==4.5.0 numpy==2.5.3 matplotlib==3.11.2 h5py==3.16.0

# No display on the cluster
ENV PATH=/opt/openEMS/venv/bin:/opt/openEMS/bin:$PATH MPLBACKEND=Agg
