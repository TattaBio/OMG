"""Example script for tokenizing the OMG dataset."""
import argparse
import datasets
import numpy as np
from transformers import AutoTokenizer
from itertools import chain


def tokenize_function(examples, tokenizer, max_seq_length: int):
    """Tokenization for mixed modality data."""
    # Remove batch dim.
    example = {k: v[0] for k, v in examples.items()}

    # Get CDS features.
    cds_seqs = example['CDS_seqs']
    cds_positions = example['CDS_position_ids']
    cds_orientations = example['CDS_orientations']

    # Get IGS features.
    igs_seqs = example['IGS_seqs']
    igs_positions = example['IGS_position_ids']

    # NOTE: This assumes the tokenizer handles all nucleotide sequences as lower-case.
    # This is to avoid collisions with amino-acids A, T, C, G (for example Glycine).
    igs_seqs = [seq.lower() for seq in igs_seqs]

    # Add orientation tokens.
    orientation_map = {True: '<+>', False: '<->'}
    cds_seqs = [f"{orientation_map[orient]}{seq}" for orient, seq in zip(
        cds_orientations, cds_seqs)]

    # Tokenize the CDS and IGS sequences.
    tokenized_cds = tokenizer(cds_seqs)
    tokenized_igs = tokenizer(igs_seqs) if igs_seqs else {
        k: [] for k in tokenized_cds.keys()}

    def _interleave(cds_elems, igs_elems):
        """Interleave the cds and igs elements based on the position ids."""
        num_elems = len(cds_elems) + len(igs_elems)
        elems = [None] * num_elems
        for i, elem in zip(cds_positions, cds_elems):
            elems[i] = elem
        for i, elem in zip(igs_positions, igs_elems):
            elems[i] = elem
        return elems

    input_ids = _interleave(
        tokenized_cds['input_ids'], tokenized_igs['input_ids'])
    attention_mask = _interleave(
        tokenized_cds['attention_mask'], tokenized_igs['attention_mask'])

    # Flatten the sequence elements
    input_ids = np.array(list(chain(*input_ids)), dtype=np.int16)
    attention_mask = np.array(list(chain(*attention_mask)), dtype=bool)

    # Pad to multiple of max_seq_length
    pad_length = max_seq_length - len(input_ids) % max_seq_length
    input_ids = np.pad(input_ids, (0, pad_length),
                       constant_values=tokenizer.pad_token_id)
    attention_mask = np.pad(attention_mask, (0, pad_length), constant_values=0)

    # Reshape to (num_examples, max_seq_length)
    input_ids = input_ids.reshape(-1, max_seq_length)
    attention_mask = attention_mask.reshape(-1, max_seq_length)

    tokenized_data = {
        'input_ids': input_ids,
        'attention_mask': attention_mask,
    }
    return tokenized_data


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--dataset_name",
        type=str,
        default='tattabio/OMG',
        help="The name of the dataset to tokenize",
    )
    parser.add_argument(
        "--tokenizer_name",
        type=str,
        default='tattabio/gLM',
        help="The name of the tokenizer to use.",
    )
    parser.add_argument(
        "--save_dir",
        type=str,
        default="tokenized_OMG",
    )
    parser.add_argument(
        "--max_seq_length",
        type=int,
        default=2048,
        help="The max sequence length after tokenization.",
    )
    parser.add_argument(
        "--num_proc",
        type=int,
        default=20,
        help="The number of processes to use for tokenization.",
    )
    args = parser.parse_args()
    return args


def main():
    args = parse_args()

    ds = datasets.load_dataset(args.dataset_name)
    tokenizer = AutoTokenizer.from_pretrained(args.tokenizer_name)

    tokenize_function_kwargs = {
        "tokenizer": tokenizer,
        "max_seq_length": args.max_seq_length,
    }

    # Tokenize the dataset.
    # NOTE: Currently, each scaffold is tokenized separately.
    # Padding can be reduced by packing scaffolds with a <sep> token.
    tokenized_ds = ds.map(
        tokenize_function,
        batched=True,
        batch_size=1,
        num_proc=args.num_proc,
        fn_kwargs=tokenize_function_kwargs,
        remove_columns=ds['train'].column_names,
    )

    tokenized_ds = tokenized_ds.with_format(
        type="torch", columns=tokenized_ds['train'].column_names)

    tokenized_ds.save_to_disk(args.save_dir, num_proc=args.num_proc)
    print(f"Saved to: {args.save_dir}")


if __name__ == "__main__":
    main()
