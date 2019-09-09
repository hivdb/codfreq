FROM python:3.6 as builder
RUN pip wheel pysam

FROM ubuntu:18.04
ENV LANG="C.UTF-8" HTSLIB_CONFIGURE_OPTIONS="--enable-plugins"
RUN apt-get update -q && \
    apt-get install -qy python3 python3-dev curl bowtie2 libperl5.26
RUN curl -sL https://bootstrap.pypa.io/get-pip.py | python3 -
COPY --from=builder /*.whl /
RUN pip install --no-cache-dir /*.whl
WORKDIR /app
COPY scripts/* /app/
COPY bt/* /bt/
