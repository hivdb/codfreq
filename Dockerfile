FROM python:3.11-alpine as builder
RUN apk add gcc make g++ zlib-dev bzip2-dev xz-dev linux-headers ncurses-dev curl-dev coreutils

FROM builder as samtools_builder
ARG SAMTOOLS_VERSION=1.17
RUN wget -O /tmp/samtools.tar.bz2 "https://github.com/samtools/samtools/releases/download/${SAMTOOLS_VERSION}/samtools-${SAMTOOLS_VERSION}.tar.bz2" && \
    tar -xjC /tmp -f /tmp/samtools.tar.bz2 && \
    cd /tmp/samtools-${SAMTOOLS_VERSION} && \
    ./configure --prefix=/usr/local/samtools && \
    make && make install

FROM builder as minimap2_installer
ARG MINIMAP2_VERSION=2.26
RUN wget -O /tmp/minimap2.tar.bz2 "https://github.com/lh3/minimap2/releases/download/v${MINIMAP2_VERSION}/minimap2-${MINIMAP2_VERSION}_x64-linux.tar.bz2" && \
    tar -xjC /tmp -f /tmp/minimap2.tar.bz2 && \
    mv /tmp/minimap2-${MINIMAP2_VERSION}_x64-linux /usr/local/minimap2

FROM builder as pysam_builder
RUN echo 'manylinux2014_compatible = True' > /usr/local/lib/python3.11/site-packages/_manylinux.py
RUN pip install 'cython==0.29.35'
RUN pip install 'pysam==0.21.0'

FROM builder as orjson_builder
RUN echo 'manylinux2014_compatible = True' > /usr/local/lib/python3.11/site-packages/_manylinux.py
RUN apk add rust cargo patchelf
RUN pip install 'orjson==3.9.1'

FROM builder as cutadapt_builder
RUN echo 'manylinux2014_compatible = True' > /usr/local/lib/python3.11/site-packages/_manylinux.py
RUN pip install 'cutadapt==4.4'

FROM builder as pydep_builder
COPY --from=pysam_builder /root/.cache/ /root/.cache/
COPY --from=orjson_builder /root/.cache/ /root/.cache/
COPY --from=cutadapt_builder /root/.cache/ /root/.cache/
RUN echo 'manylinux2014_compatible = True' > /usr/local/lib/python3.11/site-packages/_manylinux.py
RUN pip install 'cython==0.29.35'
COPY requirements.txt /codfreq/
RUN pip install -r /codfreq/requirements.txt


FROM builder as py_builder
COPY --from=pydep_builder /root/.cache/ /root/.cache/
RUN echo 'manylinux2014_compatible = True' > /usr/local/lib/python3.11/site-packages/_manylinux.py
RUN pip install 'cython==0.29.35'
COPY . /codfreq/
RUN pip install --ignore-installed --target /python-packages /codfreq
RUN mv /python-packages/bin /python-scripts

FROM builder as fastp_installer
ARG FASTP_VERSION=0.23.4
RUN wget -O /usr/local/bin/fastp http://opengene.org/fastp/fastp.${FASTP_VERSION}
RUN chmod +x /usr/local/bin/fastp

FROM builder as htslib_builder
RUN apk add automake autoconf
ARG HTSLIB_VERSION=1.17
RUN wget -O /tmp/htslib.tar.gz https://github.com/samtools/htslib/releases/download/${HTSLIB_VERSION}/htslib-${HTSLIB_VERSION}.tar.bz2
RUN tar -xjC /tmp -f /tmp/htslib.tar.gz
RUN cd /tmp/htslib-${HTSLIB_VERSION} && \
    autoreconf -i && \
    ./configure --prefix=/usr/local/htslib && \
    make && make install

FROM builder as ivar_builder
RUN apk add automake autoconf
COPY --from=htslib_builder /usr/local/htslib/ /usr/local/htslib/
ARG IVAR_VERSION=1.4.2
RUN wget -O /tmp/ivar.tar.gz https://github.com/andersen-lab/ivar/archive/refs/tags/v${IVAR_VERSION}.tar.gz
RUN tar -xzC /tmp -f /tmp/ivar.tar.gz
RUN export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/htslib/lib && \
    export LIBRARY_PATH=$LIBRARY_PATH:/usr/local/htslib/lib && \
    export CPATH=$CPATH:/usr/local/htslib/include && \
    cd /tmp/ivar-${IVAR_VERSION} && \
    chmod +x autogen.sh && \
    sh ./autogen.sh && \
    ./configure --prefix=/usr/local/ivar && \
    make && make install

FROM python:3.11-alpine
ENV LANG="C.UTF-8" \
    HTSLIB_CONFIGURE_OPTIONS="--enable-plugins"
RUN apk add --no-cache bash libc6-compat libcurl jq zip pigz
COPY --from=py_builder /python-scripts/ /usr/local/bin/
COPY --from=py_builder /python-packages/ /usr/local/lib/python3.11/site-packages/
COPY --from=orjson_builder /usr/lib/libgcc_s.so.1 /usr/lib/libgcc_s.so.1
COPY --from=samtools_builder /usr/local/samtools/ /usr/local/samtools/
COPY --from=minimap2_installer /usr/local/minimap2/ /usr/local/minimap2/
COPY --from=fastp_installer /usr/local/bin/fastp /usr/local/bin/fastp
COPY --from=htslib_builder /usr/local/htslib/lib/libhts.so.3 /usr/lib/libstdc++.so.6 /usr/lib/
COPY --from=ivar_builder /usr/local/ivar/ /usr/local/ivar/
RUN ln -s ../minimap2/minimap2 /usr/local/bin/minimap2 && \
    ln -s ../samtools/bin/samtools /usr/local/bin/samtools && \
    ln -s ../ivar/bin/ivar /usr/local/bin/ivar
WORKDIR /app
COPY bin/align-all-s3 bin/align-all-local /app/bin/
