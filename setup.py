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
    Extension(
        name='codfreq.posnas',
        sources=['codfreq/posnas.py']
    ),
    Extension(
        name='codfreq.poscodons_new',
        sources=['codfreq/poscodons_new.py']
    ),
    Extension(
        name='codfreq.pairwise_consensus',
        sources=['codfreq/pairwise_consensus.py']
    )
]


def strip_comments(line):
    if line.startswith('-i '):
        return ''
    return line.split('#', 1)[0].strip()


def req(filename):
    with open(os.path.join(os.getcwd(), filename)) as fp:
        requires = set([strip_comments(ln) for ln in fp.readlines()])
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
    # install_requires=req('requirements.txt'),
    ext_modules=cythonize(
        extensions,
        compiler_directives={'language_level': "3"}
    ),
    # tests_require=reqs('test-requirements.txt'),
    # include_package_data=True,
    # entry_points={'console_scripts': [
    #     'align = codfreq.align:main',
    # ]},
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
