from typing import Optional
import random
import logging
import os
import uuid
import datasets
from tqdm import tqdm
import gzip
from Bio import SeqIO
from datetime import datetime

BLACKLIST_PATH = '/mnt/raid0/datasets/raw/pathogenic_scaffolds_uniq.txt'


_MAX_IGS_SEQ_LENGTH = 4000  # Length in base pairs.
_MAX_CDS_SEQ_LENGTH = 45000  # Length in base pairs.
_MAX_ELEMS = 1000  # Max CDS + IGS elems before splitting the scaffold.
# Max percentage of Ns/Xs in a sequence to be considered valid.
_MAX_PCT_INVALID_SEQS = 0.2

current_time = datetime.now().strftime("%b%d_%H-%M-%S")
logfilename = "corpus_"+current_time+".log"
logging.basicConfig(filename=logfilename, filemode='w', level=logging.DEBUG)


def is_valid_sequence(seq: str, is_CDS: bool, max_pct_invalid: float = _MAX_PCT_INVALID_SEQS):
    """Checks for Xs in the sequence."""
    # Valid of less than max percent of characters are X
    if is_CDS:
        return seq.count('X') < max_pct_invalid * len(seq)
    else:
        return seq.count('N') < max_pct_invalid * len(seq)


class Scaffold():
    """Scaffold storing all CDS and IGS elements."""

    def __init__(self):
        self.num_elements = 0
        # CDS members.
        self.cds_list = []
        self.cds_ids = []
        self.cds_position_ids = []
        self.cds_orientations = []
        # IGS members.
        self.igs_list = []
        self.igs_ids = []
        self.igs_position_ids = []
        # Initialize by adding an IGS sep token as the first token.
        self.add_igs_separator()

    def as_example_dict(self):
        """Returns the scaffold as an example dict."""
        return {
            "CDS_position_ids": self.cds_position_ids,
            "IGS_position_ids": self.igs_position_ids,
            "CDS_ids": self.cds_ids,
            "IGS_ids": self.igs_ids,
            "CDS_seqs": self.cds_list,
            "IGS_seqs": self.igs_list,
            "CDS_orientations": self.cds_orientations,
        }

    def add_igs_separator(self):
        """Add a separator igs token."""
        self.igs_list.append("<sep>")
        self.igs_ids.append("SEP")
        self.igs_position_ids.append(self.num_elements)
        # Increment the number of elements.
        self.num_elements += 1

    def add_cds(self, seq: str, id: str, sample_name: str):
        """Adds a CDS elem."""
        self.cds_list.append(seq)
        self.cds_ids.append(id)
        self.cds_position_ids.append(self.num_elements)
        # Increment the number of elements.
        self.num_elements += 1
        # Encode positive direction as 1, negative direction as 0.
        orientation = 1 if "|+|" in id else 0
        self.cds_orientations.append(orientation)

    def add_igs(self, seq: str, id: str):
        """Adds a IGS elem."""
        self.igs_list.append(seq)
        self.igs_ids.append(id)
        self.igs_position_ids.append(self.num_elements)
        # Increment the number of elements.
        self.num_elements += 1


class OMGDataset(datasets.GeneratorBasedBuilder):
    """OMG dataset."""
    VERSION = datasets.Version("1.0.0")

    def __init__(self, *args, min_elems: Optional[int] = None, min_cds: Optional[int] = None, **kwargs):
        super().__init__(*args, **kwargs)
        self._min_elems = min_elems
        self._min_cds = min_cds

        with open(BLACKLIST_PATH, 'r') as file:
            lines = file.read().splitlines()
            self.blacklist = set(lines)

    def _info(self):
        # Seqs and ids are of type large_string for large datasets.
        # Orientations are stored as bool (1 for positive dir, 0 for reverse).
        return datasets.DatasetInfo(
            description="The OMG Dataset",
            features=datasets.Features(
                {
                    "CDS_position_ids": datasets.features.Sequence(datasets.Value("int32")),
                    "IGS_position_ids": datasets.features.Sequence(datasets.Value("int32")),
                    "CDS_ids": datasets.features.Sequence(datasets.Value("string")),
                    "IGS_ids": datasets.features.Sequence(datasets.Value("string")),
                    "CDS_seqs": datasets.features.Sequence(datasets.Value("large_string")),
                    "IGS_seqs": datasets.features.Sequence(datasets.Value("large_string")),
                    "CDS_orientations": datasets.features.Sequence(datasets.Value("bool")), }
            ),
            supervised_keys=None,
            homepage="TODO",
            citation="TODO",
        )

    def _split_generators(self, data_dir):
        manual_dir = data_dir.manual_dir
        all_samples = [os.path.join(manual_dir, d) for d in os.listdir(
            manual_dir) if os.path.isdir(os.path.join(manual_dir, d))]
        # Random shuffle the sample directories to avoid grouping large samples.
        random.Random(42).shuffle(all_samples)
        splits = [
            datasets.SplitGenerator(
                name=datasets.Split.TRAIN,
                # These kwargs will be passed to _generate_examples
                gen_kwargs={
                    "sample_dirs": all_samples,
                },
            )]
        return splits

    def _load_seq_dict(self, filename):
        """Parse the fasta file and return a mapping from id to sequence and remove stop codon char in case present ("*")."""
        with gzip.open(filename, 'rt') as handle:
            fasta_dict = {record.id: str(record.seq).replace(
                '*', '') for record in SeqIO.parse(handle, "fasta")}
        return fasta_dict

    def _generate_examples(self, sample_dirs):
        """Generate examples give list of sample directories."""
        for sample_dir in tqdm(sample_dirs, total=len(sample_dirs)):
            sample_name = os.path.basename(sample_dir)

            # Check if sample_dir contains all the necessary files.
            order_file = os.path.join(sample_dir, sample_name+".tsv.gz")
            cds_file = os.path.join(sample_dir, sample_name+".faa.gz")
            igs_file = os.path.join(sample_dir, sample_name+".fna.gz")
            if not os.path.isfile(order_file) or not os.path.isfile(cds_file) or not os.path.isfile(igs_file):
                logging.warning(
                    f"directory is not complete: {sample_dir} requires .tsv.gz, .faa.gz, and .fna.gz files")
                continue

            cds_dict = self._load_seq_dict(cds_file)
            igs_dict = self._load_seq_dict(igs_file)

            with gzip.open(order_file, "rt") as f:
                lines = f.readlines()
            lines = lines[1:]  # remove header
            # Original element ids. List of cds/igs keys in the correct genomic order.
            orig_ids = [line.strip().split("\t")[-1] for line in lines]
            num_elems = len(orig_ids)
            # Need at least _min_elems elements to be considered, -1 because of the separator token that is not part of the data.
            if num_elems < self._min_elems-1:
                logging.info(
                    f"{sample_dir} only has {num_elems} elements, skipping")
                continue

            # Filter elements from blacklist
            orig_ids = [elem for elem in orig_ids if elem.split(
                "|")[1] not in self.blacklist]

            current_scaff = Scaffold()

            for i, elem in enumerate(orig_ids):
                # Get the current and next element.
                next_elem = orig_ids[i+1] if i + 1 < len(orig_ids) else None
                current_scaff_name = "|".join(elem.split("|")[:2])
                next_scaff_name = "|".join(
                    next_elem.split("|")[:2]) if next_elem else None

                # Check if the current element is CDS or IGS.
                is_cds = '|CDS|' in elem
                assert is_cds or '|IG|' in elem

                max_elem_len = _MAX_CDS_SEQ_LENGTH if is_cds else _MAX_IGS_SEQ_LENGTH
                start, end = elem.split("|")[-1].split(":")
                elem_len = int(end)-int(start)
                # Get the sequence of the current element, or None.
                elem_seq = cds_dict.get(
                    elem, None) if is_cds else igs_dict.get(elem, None)
                valid_elem = (
                    # Element Was found.
                    elem_seq is not None
                    # Element does not contain many Ns/Xs
                    and is_valid_sequence(elem_seq, is_cds)
                    # Element is not too long.
                    and elem_len <= max_elem_len
                )
                # Add the element to the scaffold if it is valid.
                if valid_elem:
                    if is_cds:
                        current_scaff.add_cds(
                            seq=elem_seq, id=elem, sample_name=sample_name)
                    else:
                        current_scaff.add_igs(seq=elem_seq, id=elem)

                # Add logging for invalid elements.
                if elem_seq is None:
                    logging.info(
                        f"Could not find {elem} in {cds_file} ignoring scaffold {current_scaff_name}")
                if elem_len > max_elem_len:
                    logging.info(
                        f"Scaffold {current_scaff_name} has an element {elem} that is too long, splitting the scaffold.")

                # Whether we should start a new scaffold.
                start_new_scaffold = (
                    # If the next element is a different scaffold.
                    current_scaff_name != next_scaff_name
                    # If we have reached the max number of elements.
                    or current_scaff.num_elements >= _MAX_ELEMS
                    # If the current element has too many Ns/Xs, too long, or not found.
                    # We don't add the element to the current scaffold, and start a new one.
                    or not valid_elem
                )
                # Whether the current scaffold is valid for yielding
                valid_scaffold = (
                    # Has enough elements.
                    current_scaff.num_elements >= self._min_elems
                    # Has enough CDS.
                    and len(current_scaff.cds_ids) >= self._min_cds
                    # The current element was found. Otherwise we discard the scaffold.
                    and elem_seq is not None
                )

                if start_new_scaffold:
                    if valid_scaffold:  # Yield the current scaffold.
                        unique_id = str(uuid.uuid4())
                        key = f"{sample_name}-{unique_id}"
                        yield key, current_scaff.as_example_dict()
                    # Discard and reset the current scaffold whether or not we yielded.
                    current_scaff = Scaffold()
