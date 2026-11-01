FROM --platform=linux/amd64 rocker/shiny:4.6.1

# Install system dependencies, Nextflow, Java (required by Nextflow), tmux, git, python, and Docker CLI
RUN apt-get update && apt-get install -y --no-install-recommends \
    curl \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    libgit2-dev \
    libuv1-dev \
    git \
    tmux \
    default-jre-headless \
    python3 \
    python3-pip \
    python3-venv \
    singularity-container \
    squashfuse \
    fuse \
    && rm -rf /var/lib/apt/lists/*

# Install and pin Nextflow version 25.04.7 globally
ENV NXF_VER=25.04.7
RUN curl -fsSL https://github.com/nextflow-io/nextflow/releases/download/v25.04.7/nextflow -o /usr/local/bin/nextflow && \
    chmod +x /usr/local/bin/nextflow

# Use Posit PPM binary mirror for R packages
ENV R_REPOS="https://packagemanager.posit.co/cran/__linux__/jammy/latest"

WORKDIR /app

# Copy renv configuration files if available
COPY renv.lock .Rprofile ./
COPY renv/activate.R renv/

# Restore R environment using renv
RUN R -e "options(repos = c(CRAN = Sys.getenv('R_REPOS'))); install.packages('renv'); renv::restore()"

# Copy application source files
COPY . .

# Create volume mount points and ensure /tmp permissions for tmux sockets
RUN mkdir -p /mnt /work /output && chmod 1777 /tmp

EXPOSE 3838

CMD ["R", "-e", "shiny::runApp('/app', host='0.0.0.0', port=3838)"]
