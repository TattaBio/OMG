#!/usr/bin/env python
"""
Preprocessing script for JGI IMG files, to the format expected by omg_generator.py.

Adapted from jgi_preprocessor.py.

Usage: python mgnify_preprocessor.py erz_list.txt
Where erz_list.txt is a file containing the list of ERZs to be preprocessed.
The script will create a directory called MGnify in the current directory and store the processed files:
- ERZ.fna.gz: FASTA file containing the intergenic regions.
- ERZ.faa.gz: FASTA file containing the CDS regions.
- ERZ.tsv.gz: TSV file containing the coordinates of the CDS and intergenic regions.
This script requires all ERZ files to be stored in the MGnify_RAW directory, where each ERZ has its own directory.
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
from tqdm import tqdm
import sys
import os

ERZ_LISTFILE = sys.argv[1]
MIN_SCAFFOLD_LENGTH = 2_000
MIN_INTERGENIC_LENGTH = 0
REMOVE_EDGE_CDS = True

# Paths
DIRNAME = "MGnify"  # where all processed files will be stored.
RAW_DIRNAME = "MGnify_RAW"  # where all raw files are stored.


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


def read_fasta_full_header(filepath, uppercase=False, strip_n=False, compress=False):
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
    # Read ERZ list
    with open(ERZ_LISTFILE, "r") as f:
        erz_list = f.read().splitlines()
    for input_erz in tqdm(erz_list, total=len(erz_list)):
        # Paths
        output_directory = Path(DIRNAME) / input_erz
        output_tsv = output_directory / f"{input_erz}.tsv.gz"
        output_fna = output_directory / f"{input_erz}.fna.gz"
        output_faa = output_directory / f"{input_erz}.faa.gz"

        input_fna = f"{RAW_DIRNAME}/{input_erz}/{input_erz}_FASTA.fasta.gz"
        input_faa = f"{RAW_DIRNAME}/{input_erz}/{input_erz}_FASTA_CDS.faa.gz"
        if not os.path.isfile(input_fna) or not os.path.isfile(input_faa):
            print(f"Missing input files for {input_erz}")
            continue
        # Create output directory
        if output_directory.exists():
            for filepath in output_directory.glob("*"):
                filepath.unlink()
        else:
            output_directory.mkdir(parents=True)

        selected_sequences = set()
        for record in read_fasta(input_fna):
            if len(record) - record.count("N") >= MIN_SCAFFOLD_LENGTH:
                selected_sequences.add(record.accession)

        cds_coordinates = defaultdict(list)
        for record in read_fasta_full_header(input_faa):
            gene = record.accession
            elems = record.header.split(" # ")

            if len(elems) == 5:  # prodigal results
                scaffold = "_".join(gene.split("_")[:-1])
                if elems[3] == "1":
                    strand = "+"
                elif elems[3] == "-1":
                    strand = "-"
                start = int(elems[1])
                end = int(elems[2])
                cds_coordinates[scaffold].append((gene, strand, [start, end]))
            else:  # FRAGENESCAN results
                id_parts = record.accession.split("_")
                scaffold = "_".join(id_parts[:-3])
                strand = id_parts[-1]
                start = int(id_parts[-3])
                end = int(id_parts[-2])
                cds_coordinates[scaffold].append((gene, strand, [start, end]))

        # Remove scaffolds with less than 4 CDS
        cds_coordinates = {k: v for k, v in cds_coordinates.items() if len(v) >= 4}

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

        # Remove edge CDS
        if REMOVE_EDGE_CDS:
            for k, v in cds_coordinates.items():
                min_index = min(range(len(v)), key=lambda i: v[i][2][0])
                v.pop(min_index)
                max_index = max(range(len(v)), key=lambda i: v[i][2][1])
                v.pop(max_index)

        with gzip.open(output_fna, "wt") as fout:
            for record in read_fasta(input_fna, uppercase=True):
                if record.accession in intergenic_coordinates:
                    for name, strand, (start, end) in intergenic_coordinates[
                        record.accession
                    ]:
                        interval_record = record[start - 1 : end]
                        interval_record._header = f"{input_erz}|{record.accession}|IG|{name}|{strand}|{start}:{end}"
                        # assert that the length of the sequence is equal to the length of the interval
                        assert len(interval_record) == end - start + 1
                        # assert interval sequence is longer than 0
                        assert len(interval_record) >= MIN_INTERGENIC_LENGTH

                        fout.write(str(interval_record))

        genes_dict = {}
        for k, v in cds_coordinates.items():
            for gene, strand, (start, end) in v:
                genes_dict[gene] = [k, strand, start, end]
        with gzip.open(output_faa, "wt") as fout:
            for record in read_fasta(input_faa, uppercase=True):
                if record.accession in genes_dict:
                    scaffold, strand, start, end = genes_dict[record.accession]
                    record._header = f"{input_erz}|{scaffold}|CDS|{record.accession}|{strand}|{start}:{end}"

                    # assert scaffold is a substring of record.accession
                    assert scaffold in record.accession
                    fout.write(str(record))

        all_coordinates = {}
        for k in cds_coordinates.keys():
            all_coordinates[k] = []
            for name, strand, (start, end) in cds_coordinates[k]:
                header = f"{input_erz}|{k}|CDS|{name}|{strand}|{start}:{end}"
                all_coordinates[k].append(
                    [input_erz, k, name, "CDS", strand, start, end, header]
                )
            for name, strand, (start, end) in intergenic_coordinates[k]:
                header = f"{input_erz}|{k}|IG|{name}|{strand}|{start}:{end}"
                all_coordinates[k].append(
                    [input_erz, k, name, "IG", strand, start, end, header]
                )
        # Sort values by start
        for k, v in all_coordinates.items():
            all_coordinates[k] = sorted(v, key=lambda x: x[-3])

        with gzip.open(output_tsv, "wt") as fout:
            fout.write(
                "input_erz\tscaffold\taccession\tfeature_type\tstrand\tstart\tend\theader\n"
            )
            for v in all_coordinates.values():
                for row in v:
                    row = "\t".join([str(i) for i in row])
                    fout.write(f"{row}\n")


if __name__ == "__main__":
    main()
