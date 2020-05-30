import os
from .base import docker_execute


def sam2bam(sam, bam, outdir):
    docker_execute([
        'sam2bam',
        '/shared/output/{}'.format(sam),
        '/shared/output/{}'.format(bam)
    ], {
        outdir: {
            'bind': '/shared/output',
            'mode': 'rw'
        }
    })


def bam2fasta(refseq, bam, fasta, outdir):
    refdir, refname = os.path.split(refseq)
    docker_execute([
        'bam2fasta',
        '/shared/ref/{}'.format(refname),
        '/shared/output/{}'.format(bam),
        '/shared/output/{}'.format(fasta)
    ], {
        refdir: {
            'bind': '/shared/ref',
            'mode': 'rw'
        },
        outdir: {
            'bind': '/shared/output',
            'mode': 'rw'
        }
    })
