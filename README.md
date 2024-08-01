
<h1 align="center">OMG:  An Open MetaGenomic Dataset</h1>

<p align="center">
    <a href="https://www.biorxiv.org/content/XXX">
        <img alt="bioRxiv URL" src="https://img.shields.io/badge/bioRxiv-XXX.svg">
    </a>
    <a href="https://huggingface.co/datasets/tattabio/OMG">
        <img alt="Huggingface URL" src="https://huggingface.co/datasets/huggingface/badges/resolve/main/dataset-on-hf-md.svg">
    </a>
</p>

<h4 align="center">
    <p>
        <a href="#usage">Usage</a> |
        <a href="#citing">Citing</a>
    <p>
</h4>

<h3 align="center">
    <a href="https://huggingface.co/spaces/dgeb"><img style="float: middle; padding: 10px 10px 10px 10px;" width="100" height="100" src="./docs/images/tatta_logo.png" /></a>
</h3>

The OpenMetaGenome (OMG) dataset is ...

## Usage

Download OMG from the Huggingface Hub: 
```python
import datasets

ds = datasets.load('tattabio/OMG')
```

## Format

Each row of the dataset represents a genomic scaffold, as a list of amino acid coding sequences (CDS) and nucleotide intergenic sequences (IGS). 

| Feature | Description | Example |
|---|---|---|
| `CDS_seqs` | A list of strings representing the amino acid sequences of the CDS. | `['MALTKVEKRNR...', 'MLGIDNIERVK...', 'MATIKVKQVR...', 'MNLSNIKPAS...']` |
| `IGS_seqs` | A list of strings representing the nucleotide sequences of the IGS. | `['<sep>', 'AATTTAAGGAA', 'TTTTAAAAGTATCGAAAT', 'TTTTTAAAGAAAA', 'AAATAATCTT', 'AAAAAATA']` |
| `CDS_position_ids` | A list of integers representing the positions of the coding sequences (CDS) in the sequence. | `[2, 4, 6, 8]` |
| `IGS_position_ids` | A list of integers representing the positions of the intergenic sequences (IGS) in the sequence. | `[0, 1, 3, 5, 7, 9]` |
| `CDS_ids` | A list of strings representing the identifiers of the CDS. The format of the identifier is `accession|contig_id|feature_type|gene_id|strand|start:end`.  | `['7000000126|C1821366|CDS|C1821366__gene_115413|+|84:437', '7000000126|C1821366|CDS|C1821366__gene_115414|+|456:977', '7000000126|C1821366|CDS|C1821366__gene_115415|+|991:1167', '7000000126|C1821366|CDS|C1821366__gene_115416|+|1243:1689']` |
| `IGS_ids` | A list of strings representing the identifiers of the IGS. The format of the identifier is `SEP` for the separator or `accession|contig_id|feature_type|feature_id|strand|start:end`. | `['SEP', '7000000126|C1821366|IG|IG_000001|+|73:83', '7000000126|C1821366|IG|IG_000002|+|438:455', '7000000126|C1821366|IG|IG_000003|+|978:990', '7000000126|C1821366|IG|IG_000004|+|1168:1242', '7000000126|C1821366|IG|IG_000005|+|1690:1697']` |
| `CDS_orientations` | A list of booleans indicating the orientation of the CDS. `True` represents the forward strand, and `False` represents the reverse strand. | `[True, True, True, True]` |

## Acknowledgements



## Citing

OMG was introduced in "[The OMG dataset: An Open MetaGenomic corpus for mixed-modality genomic language modeling
]()", feel free to cite:

```
@article{WestRoberts2024,
  title = {Diverse Genomic Embedding Benchmark for functional evaluation across the tree of life},
  url = {http://dx.doi.org/10.1101/2024.07.10.602933},
  DOI = {10.1101/2024.07.10.602933},
  publisher = {Cold Spring Harbor Laboratory},
  author = {West-Roberts,  Jacob and Kravitz,  Joshua and Jha,  Nishant and Cornman,  Andre and Hwang,  Yunha},
  year = {2024},
  month = jul 
}
```
