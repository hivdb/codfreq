FROM ubuntu:18.04
ENV LANG="C.UTF-8" HTSLIB_CONFIGURE_OPTIONS="--enable-plugins"
RUN apt-get update -qq && \
    apt-get install -qqy python3 python3-dev curl
RUN curl -sL https://bootstrap.pypa.io/get-pip.py | python3 -
RUN pip install --no-cache-dir pysam
WORKDIR /app
COPY scripts/* /app/
