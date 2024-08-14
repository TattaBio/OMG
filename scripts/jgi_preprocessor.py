#!/usr/bin/env python
"""
Preprocessing script for JGI IMG files, to the format expected by omg_generator.py.

This script is written by Antônio Camargo (https://github.com/apcamargo).

Usage: python jgi_preprocessor.py TAXON_OID
Where TAXON_OID is the taxon oid of the genome/sample to be preprocessed.
The script will create a directory named IMG/TAXON_OID and store the following files:
- TAXON_OID.fna: fasta file with intergenic regions
- TAXON_OID.faa: fasta file with CDS
- TAXON_OID.tsv: TSV file containing the coordinates of the CDS and intergenic regions.
This script requires all raw IMG files to be stored in the IMG_RAW directory, where each TAXON_OID has its own directory.
"""

import bz2
import gzip
import lzma
import textwrap
from collections import defaultdict
from contextlib import contextmanager
from copy import deepcopy
from enum import Enum, auto
from pathlib import Path
import sys

TAXON_OID = sys.argv[1]
MIN_SCAFFOLD_LENGTH = 2_000
MIN_INTERGENIC_LENGTH = 0
REMOVE_EDGE_CDS = True

# Paths
DIRNAME = "IMG"  # where all processed files will be stored.
RAW_DIRNAME = "IMG_RAW"  # where all raw files are stored.
OUTPUT_DIRECTORY = Path(DIRNAME) / TAXON_OID
OUTPUT_TSV = OUTPUT_DIRECTORY / f"{TAXON_OID}.tsv"
OUTPUT_FNA = OUTPUT_DIRECTORY / f"{TAXON_OID}.fna"
OUTPUT_FAA = OUTPUT_DIRECTORY / f"{TAXON_OID}.faa"
INPUT_FNA = f"{RAW_DIRNAME}/{TAXON_OID}.fna"
INPUT_FAA = f"{RAW_DIRNAME}/{TAXON_OID}.faa"
INPUT_GFF = f"{RAW_DIRNAME}/{TAXON_OID}.gff"


class Compression(Enum):
    bzip2 = auto()
    gzip = auto()
    xz = auto()
    zstd = auto()
    uncompressed = auto()


def is_compressed(filepath: Path):
    with open(filepath, "rb") as fin:
        signature = fin.peek(8)[:8]
        if tuple(signature[:2]) == (0x1F, 0x8B):
            return Compression.gzip
        elif tuple(signature[:3]) == (0x42, 0x5A, 0x68):
            return Compression.bzip2
        elif tuple(signature[:7]) == (0xFD, 0x37, 0x7A, 0x58, 0x5A, 0x00, 0x00):
            return Compression.xz
        elif tuple(signature[:4]) == (0x28, 0xB5, 0x2F, 0xFD):
            return Compression.zstd
        else:
            return Compression.uncompressed


@contextmanager
def open_file(filepath):
    filepath_compression = is_compressed(filepath)
    if filepath_compression == Compression.gzip:
        fin = gzip.open(filepath, "rt")
    elif filepath_compression == Compression.bzip2:
        fin = bz2.open(filepath, "rt")
    elif filepath_compression == Compression.xz:
        fin = lzma.open(filepath, "rt")
    else:
        fin = open(filepath, "r")
    try:
        yield fin
    finally:
        fin.close()


class Sequence:
    def __init__(self, header: str, seq: str, compress: bool = False):
        self._compress = compress
        self._header = header
        if self._compress:
            self._seq = gzip.compress(seq.encode("ascii"), 1)
        else:
            self._seq = seq.encode("ascii")

    @property
    def header(self):
        return self._header

    @property
    def accession(self):
        return self._header.split()[0]

    @property
    def seq(self):
        if self._compress:
            return gzip.decompress(self._seq).decode()
        else:
            return self._seq.decode()

    @property
    def seq_ascii(self):
        return self.seq.upper().encode("ascii")

    def count(self, substring: str):
        return self.seq.count(substring)

    def rc(self):
        tab = self.seq.maketrans("ACTGNactgn", "TGACNtgacn")
        return Sequence(self.header, self.seq.translate(tab)[::-1], self._compress)

    def has_dtr(self, min_length: int = 21):
        substring = self.seq.casefold()[:min_length]
        pos = self.seq.casefold().rfind(substring)
        if pos < len(self) / 2:
            return False
        substring = self.seq.casefold()[pos:]
        return self.seq.casefold()[: len(substring)] == substring

    def has_itr(self, min_len: int = 21):
        rev = self.rc().seq
        return self.seq.casefold()[:min_len] == rev.casefold()[:min_len]

    def __str__(self):
        return f">{self.header}\n{textwrap.fill(self.seq, 60)}\n"

    def __repr__(self):
        if len(self) > 40:
            start = self.seq[:34]
            end = self.seq[-3:]
            seq = f"{start}...{end}"
        else:
            seq = self.seq
        return f"Sequence({self.accession}, {seq})"

    def __len__(self):
        return len(self.seq)

    def __getitem__(self, k: int):
        return Sequence(self.header, self.seq[k], self._compress)

    def __eq__(self, other: object):
        if other.__class__ is self.__class__:
            return self.seq.casefold() == other.seq.casefold()
        elif other.__class__ is str:
            return self.seq.casefold() == other.casefold()
        return NotImplemented

    def __hash__(self):
        return hash(self.seq.casefold())

    def __add__(self, other: object):
        if other.__class__ is not self.__class__:
            return NotImplemented
        compress = other._compress or self._compress
        return Sequence(
            f"{self.accession}+{other.accession}", f"{self.seq}{other.seq}", compress
        )


def read_fasta(filepath, uppercase=False, strip_n=False, compress=False):
    with open_file(filepath) as fin:
        last = None
        while True:
            if not last:
                for l in fin:
                    if l[0] == ">":
                        last = l[:-1]
                        break
            if not last:
                break
            name, seqs, last = last[1:], [], None
            for l in fin:
                if l[0] == ">":
                    last = l[:-1]
                    break
                seqs.append(l[:-1])
            seqs = "".join(seqs)
            if uppercase:
                seqs = seqs.upper()
            if strip_n:
                seqs = seqs.strip("nN")
            if len(seqs):
                yield Sequence(name, seqs, compress)
            if not last:
                break


def merge_intervals(intervals):
    intervals = deepcopy(intervals)
    intervals.sort()
    stack = []
    stack.append(intervals[0])
    for i in intervals[1:]:
        if stack[-1][0] <= i[0] <= stack[-1][-1]:
            stack[-1][-1] = max(stack[-1][-1], i[-1])
        else:
            stack.append(i)
    return stack


def interval_difference(a, b):
    values_b = []
    for interval in b:
        values_b.extend(range(interval[0], interval[1] + 1))
    values_b = set(values_b)
    result = []
    for interval in a:
        values_a = set(range(interval[0], interval[1] + 1))
        diff = sorted(values_a.difference(values_b))
        intervals = []
        i = 0
        while i < len(diff):
            start = diff[i]
            while i < len(diff) - 1 and diff[i + 1] == diff[i] + 1:
                i += 1
            end = diff[i]
            intervals.append([start, end])
            i += 1
        result.extend(intervals)
    return result


def main():
    # Create output directory
    if OUTPUT_DIRECTORY.exists():
        for filepath in OUTPUT_DIRECTORY.glob("*"):
            filepath.unlink()
    else:
        OUTPUT_DIRECTORY.mkdir(parents=True)

    selected_sequences = set()
    for record in read_fasta(INPUT_FNA):
        if len(record) - record.count("N") >= MIN_SCAFFOLD_LENGTH:
            selected_sequences.add(record.accession)

    cds_coordinates = defaultdict(list)
    with open_file(INPUT_GFF) as fin:
        for line in fin:
            if line.startswith("#"):
                continue
            elif line.split("\t")[2] != "CDS":
                continue
            gene = line.split("\t")[8].split(";")[0].split("=")[1]
            scaffold = line.split("\t")[0]
            if scaffold in selected_sequences:
                start = int(line.split("\t")[3])
                end = int(line.split("\t")[4])
                strand = line.split("\t")[6]
                cds_coordinates[scaffold].append((gene, strand, [start, end]))

    # Remove scaffolds with less than 4 CDS
    cds_coordinates = {k: v for k, v in cds_coordinates.items() if len(v) >= 4}

    # Remove edge CDS
    if REMOVE_EDGE_CDS:
        for k, v in cds_coordinates.items():
            min_index = min(range(len(v)), key=lambda i: v[i][2][0])
            v.pop(min_index)
            max_index = max(range(len(v)), key=lambda i: v[i][2][1])
            v.pop(max_index)

    # Get intergenic regions
    mix_max_coordinates = {}
    for k, v in cds_coordinates.items():
        min_coord = min([j for i in v for j in i[2]])
        max_coord = max([j for i in v for j in i[2]])
        mix_max_coordinates[k] = [min_coord, max_coord]

    intergenic_coordinates = {}
    for k, v in cds_coordinates.items():
        cds_intervals = [i[2] for i in v]
        cds_intervals = merge_intervals(cds_intervals)
        intergenic_coordinates[k] = interval_difference(
            [mix_max_coordinates[k]], cds_intervals
        )
    for k, v in intergenic_coordinates.items():
        v = [
            (f"IG_{e:06d}", "+", i)
            for e, i in enumerate(v, 1)
            if i[1] - i[0] + 1 >= MIN_INTERGENIC_LENGTH
        ]
        intergenic_coordinates[k] = v

    with open(OUTPUT_FNA, "w") as fout:
        for record in read_fasta(INPUT_FNA, uppercase=True):
            if record.accession in intergenic_coordinates:
                for name, strand, (start, end) in intergenic_coordinates[
                    record.accession
                ]:
                    interval_record = record[start - 1 : end]
                    interval_record._header = f"{TAXON_OID}|{record.accession}|IG|{name}|{strand}|{start}:{end}"
                    fout.write(str(interval_record))

    genes_dict = {}
    for k, v in cds_coordinates.items():
        for gene, strand, (start, end) in v:
            genes_dict[gene] = [k, strand, start, end]
    with open(OUTPUT_FAA, "w") as fout:
        for record in read_fasta(INPUT_FAA, uppercase=True):
            if record.accession in genes_dict:
                scaffold, strand, start, end = genes_dict[record.accession]
                record._header = f"{TAXON_OID}|{scaffold}|CDS|{record.accession}|{strand}|{start}:{end}"
                fout.write(str(record))

    all_coordinates = {}
    for k in cds_coordinates.keys():
        all_coordinates[k] = []
        for name, strand, (start, end) in cds_coordinates[k]:
            header = f"{TAXON_OID}|{k}|CDS|{name}|{strand}|{start}:{end}"
            all_coordinates[k].append(
                [TAXON_OID, k, name, "CDS", strand, start, end, header]
            )
        for name, strand, (start, end) in intergenic_coordinates[k]:
            header = f"{TAXON_OID}|{k}|IG|{name}|{strand}|{start}:{end}"
            all_coordinates[k].append(
                [TAXON_OID, k, name, "IG", strand, start, end, header]
            )
    # Sort values by start
    for k, v in all_coordinates.items():
        all_coordinates[k] = sorted(v, key=lambda x: x[-3])

    with open(OUTPUT_TSV, "w") as fout:
        fout.write(
            "taxon_oid\tscaffold\taccession\tfeature_type\tstrand\tstart\tend\theader\n"
        )
        for v in all_coordinates.values():
            for row in v:
                row = "\t".join([str(i) for i in row])
                fout.write(f"{row}\n")


if __name__ == "__main__":
    main()
