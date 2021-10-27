FROM python:3.9-alpine as samtools_builder
RUN apk add gcc make g++ zlib-dev bzip2-dev xz-dev linux-headers ncurses-dev curl-dev coreutils
ARG SAMTOOLS_VERSION=1.10
RUN wget -O /tmp/samtools.tar.bz2 "https://github.com/samtools/samtools/releases/download/1.10/samtools-${SAMTOOLS_VERSION}.tar.bz2" && \
    tar -xjC /tmp -f /tmp/samtools.tar.bz2 && \
    cd /tmp/samtools-${SAMTOOLS_VERSION} && \
    ./configure --prefix=/usr/local/samtools && \
    make && make install

FROM python:3.9-alpine as minimap2_installer
RUN apk add gcc make g++ zlib-dev bzip2-dev xz-dev linux-headers ncurses-dev curl-dev coreutils
ARG MINIMAP2_VERSION=2.22
RUN wget -O /tmp/minimap2.tar.bz2 "https://github.com/lh3/minimap2/releases/download/v${MINIMAP2_VERSION}/minimap2-${MINIMAP2_VERSION}_x64-linux.tar.bz2" && \
    tar -xjC /tmp -f /tmp/minimap2.tar.bz2 && \
    mv /tmp/minimap2-${MINIMAP2_VERSION}_x64-linux /usr/local/minimap2

FROM python:3.9-alpine as pysam_builder
RUN apk add gcc make g++ zlib-dev bzip2-dev xz-dev linux-headers ncurses-dev curl-dev coreutils
RUN echo 'manylinux2014_compatible = True' > /usr/local/lib/python3.9/site-packages/_manylinux.py
RUN pip install 'cython>=0.29.23'
RUN pip install 'pysam'

FROM python:3.9-alpine as py_builder
RUN apk add gcc make g++ zlib-dev bzip2-dev xz-dev linux-headers ncurses-dev curl-dev coreutils
COPY --from=pysam_builder /root/.cache/ /root/.cache/
RUN echo 'manylinux2014_compatible = True' > /usr/local/lib/python3.9/site-packages/_manylinux.py
RUN pip install 'cython>=0.29.23'
COPY . /codfreq/
RUN pip install --ignore-installed --target /python-packages /codfreq
RUN mv /python-packages/bin /python-scripts

FROM python:3.9-alpine
ENV LANG="C.UTF-8" \
    HTSLIB_CONFIGURE_OPTIONS="--enable-plugins"
RUN apk add --no-cache bash libc6-compat libcurl jq
COPY --from=py_builder /python-scripts/ /usr/local/bin/
COPY --from=py_builder /python-packages/ /usr/local/lib/python3.9/site-packages/
COPY --from=samtools_builder /usr/local/samtools/ /usr/local/samtools/
COPY --from=minimap2_installer /usr/local/minimap2/ /usr/local/minimap2/
RUN ln -s ../minimap2/minimap2 /usr/local/bin/minimap2 && \
    ln -s ../samtools/bin/samtools /usr/local/bin/samtools
WORKDIR /app
COPY bin/align-all-s3 bin/align-all-local /app/bin/
