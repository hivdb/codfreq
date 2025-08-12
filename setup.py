#! /usr/bin/env python

import setuptools
from setuptools.extension import Extension
from Cython.Build import cythonize


extensions = [
    Extension(
        name='codfreq.samfile_helper',
        sources=['codfreq/samfile_helper.py']
    ),
    Extension(
        name='codfreq.paired_reads',
        sources=['codfreq/paired_reads.py'],
        # define_macros=[('CYTHON_TRACE', '1')]
    ),
    Extension(
        name='codfreq.posnas',
        sources=['codfreq/posnas.py'],
        # define_macros=[('CYTHON_TRACE', '1')]
    ),
    Extension(
        name='codfreq.poscodons',
        sources=['codfreq/poscodons.py'],
        # define_macros=[('CYTHON_TRACE', '1')]
    ),
    Extension(
        name='codfreq.sam2codfreq',
        sources=['codfreq/sam2codfreq.py'],
        # define_macros=[('CYTHON_TRACE', '1')]
    ),
    Extension(
        name='codfreq.codonalign_consensus',
        sources=['codfreq/codonalign_consensus.py'],
        # define_macros=[('CYTHON_TRACE', '1')]
    ),
    Extension(
        name='codfreq.sam2consensus',
        sources=['codfreq/sam2consensus.py'],
        # define_macros=[('CYTHON_TRACE', '1')]
    )
]


if __name__ == '__main__':
    setuptools.setup(
        ext_modules=cythonize(  # type: ignore[no-untyped-call]
            extensions,
            compiler_directives={
                'language_level': '3',
                'profile': False,
                'linetrace': False
            }
        )
    )
