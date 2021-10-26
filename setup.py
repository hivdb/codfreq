#! /usr/bin/env python
# -*- coding: UTF-8 -*-

import os
# import re
# import ast
import setuptools
from setuptools.extension import Extension
from Cython.Build import cythonize

# _version_re = re.compile(r'VERSION\s+=\s+(.*)')
#
# with open('codfreq/version.py', 'rb') as f:
#     version = str(ast.literal_eval(_version_re.search(
#         f.read().decode('utf-8')).group(1)))


extensions = [
    # Extension(
    #     name='codfreq.codfreq_types',
    #     sources=['codfreq/codfreq_types.py']
    # ),
    Extension(
        name='codfreq.posnas',
        sources=['codfreq/posnas.py']
    ),
    Extension(
        name='codfreq.poscodons',
        sources=['codfreq/poscodons.py']
    ),
    Extension(
        name='codfreq.sam2codfreq',
        sources=['codfreq/sam2codfreq.py']
    ),
    Extension(
        name='codfreq.codonalign_consensus',
        sources=['codfreq/codonalign_consensus.py']
    ),
    Extension(
        name='codfreq.sam2consensus',
        sources=['codfreq/sam2consensus.py']
    )
]


def strip_comments(line):
    return line.split('#', 1)[0].strip()


def pep508(line):
    if line.startswith('-i '):
        return ''
    if not line:
        return line
    if '://' in line:
        url, package_name = line.split('#egg=', 1)
        return '{} @ {}'.format(package_name.strip(), url.strip())
    return line


def req(filename):
    with open(os.path.join(os.getcwd(), filename)) as fp:
        requires = set(
            strip_comments(pep508(ln))
            for ln in fp.readlines()
        )
        requires -= set([''])
    return list(requires)


setup_params = dict(
    name="codfreq",
    version="0.1",
    url="https://github.com/hivdb/codfreq",
    author='Philip Tzou',
    author_email="philiptz@stanford.edu",
    packages=[
        'codfreq', 'codfreq/cmdwrappers'
    ],
    package_data={
        'codfreq': ['py.typed']
    },
    install_requires=req('requirements.txt'),
    ext_modules=cythonize(
        extensions,
        compiler_directives={
            'language_level': '3',
            'profile': False
        }
    ),
    # tests_require=reqs('test-requirements.txt'),
    # include_package_data=True,
    entry_points={'console_scripts': [
        'fastq2codfreq = codfreq.align:align_cmd',
        'make-response = codfreq.make_response:make_response'
    ]},
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Developers',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: MIT License',
        'Operating System :: OS Independent',
        'Programming Language :: Python',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
        'Topic :: Scientific/Engineering :: Bio-Informatics'],
    # test_suite="nose.collector",
    zip_safe=True)

if __name__ == '__main__':
    setuptools.setup(**setup_params)
