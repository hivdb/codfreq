FROM python:3.7-alpine as builder
RUN apk add gcc make g++ zlib-dev bzip2-dev xz-dev linux-headers
COPY requirements.txt /tmp/requirements.txt
# RUN echo 'manylinux2014_compatible = True' > /usr/local/lib/python3.7/site-packages/_manylinux.py
RUN mkdir /extras
RUN pip wheel -r /tmp/requirements.txt
RUN apk add ncurses-dev curl-dev coreutils
ARG SAMTOOLS_VERSION=1.10
RUN wget -O /tmp/samtools.tar.bz2 "https://github.com/samtools/samtools/releases/download/1.10/samtools-${SAMTOOLS_VERSION}.tar.bz2" && \
    tar -xjC /tmp -f /tmp/samtools.tar.bz2 && \
    cd /tmp/samtools-${SAMTOOLS_VERSION} && \
    ./configure --prefix=/usr/local/samtools && \
    make && make install

FROM python:3.7-alpine
ENV LANG="C.UTF-8" \
    HTSLIB_CONFIGURE_OPTIONS="--enable-plugins"
RUN apk add --no-cache bash libc6-compat libcurl jq
COPY --from=builder /*.whl /wheels/
RUN echo 'manylinux2014_compatible = True' > /usr/local/lib/python3.7/site-packages/_manylinux.py && \
    pip install --no-cache-dir /wheels/*.whl
COPY --from=builder /usr/local/samtools/ /usr/local/samtools/
ARG MINIMAP2_VERSION=2.17
RUN wget -O /tmp/minimap2.tar.bz2 "https://github.com/lh3/minimap2/releases/download/v${MINIMAP2_VERSION}/minimap2-${MINIMAP2_VERSION}_x64-linux.tar.bz2" && \
    tar -xjC /tmp -f /tmp/minimap2.tar.bz2 && \
    mv /tmp/minimap2-${MINIMAP2_VERSION}_x64-linux /usr/local/minimap2 && \
    ln -s ../minimap2/minimap2 /usr/local/bin/minimap2 && \
    ln -s ../samtools/bin/samtools /usr/local/bin/samtools
WORKDIR /app
COPY bin/align-all-docker /app/bin/
COPY codfreq/ /app/codfreq/
