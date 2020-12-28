def load(fp):
    sequences = []
    header = None
    curseq = ''
    for line in fp:
        if line.startswith('>'):
            if header and curseq:
                sequences.append({
                    'header': header,
                    'sequence': curseq.upper()
                })
            header = line[1:].strip()
            curseq = ''
        elif line.startswith('#'):
            continue
        else:
            curseq += line.strip()
    if header and curseq:
        sequences.append({
            'header': header,
            'sequence': curseq.upper()
        })
    return sequences
