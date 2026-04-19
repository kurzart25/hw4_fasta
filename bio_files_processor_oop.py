"""
bio_files_processor_oop.py

OOP-based utilities for:
- FASTA: converting multiline format to single-line format
- BLAST (text output): extracting the top-hit description for each query
- GenBank: selecting neighboring CDS features and exporting protein sequences to FASTA

Design principles:
- using pathlib instead of os.path
- configuration via dataclasses
- stateful classes with reusable methods
- strict input validation
- GenBank parsing implemented with Biopython
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple, Union


PathLike = Union[str, Path]


@dataclass(frozen=True, slots=True)
class FastaOnelineConfig:
    encoding: str = "utf-8"
    line_ending: str = "\n"
    suffix: str = "_oneline"
    default_ext: str = ".fasta"


class FastaOnelineConverter:
    """Stateful converter: multi-line FASTA -> one-line FASTA."""

    def __init__(self, config: Optional[FastaOnelineConfig] = None) -> None:
        self.config = config or FastaOnelineConfig()

    def convert(self, input_fasta: PathLike, output_fasta: Optional[PathLike] = None) -> Path:
        inp = Path(input_fasta)
        if not inp.exists():
            raise FileNotFoundError(f"Input FASTA not found: {inp}")
        if inp.is_dir():
            raise IsADirectoryError(f"Input FASTA is a directory: {inp}")

        out = self._resolve_output_path(inp, output_fasta)

        header: Optional[str] = None
        seq_parts: List[str] = []

        with (
            inp.open("r", encoding=self.config.encoding) as fin,
            out.open("w", encoding=self.config.encoding) as fout,
        ):
            for raw in fin:
                line = raw.rstrip("\r\n")
                if not line:
                    continue

                if line.startswith(">"):
                    if header is not None:
                        fout.write(header + self.config.line_ending)
                        fout.write("".join(seq_parts) + self.config.line_ending)
                    header = line
                    seq_parts = []
                else:
                    seq_parts.append(line.strip())

            if header is not None:
                fout.write(header + self.config.line_ending)
                fout.write("".join(seq_parts) + self.config.line_ending)

        return out

    def _resolve_output_path(self, inp: Path, output_fasta: Optional[PathLike]) -> Path:
        if output_fasta is not None:
            return Path(output_fasta)
        ext = inp.suffix if inp.suffix else self.config.default_ext
        return inp.with_name(f"{inp.stem}{self.config.suffix}{ext}")


@dataclass(frozen=True, slots=True)
class BlastParseConfig:
    encoding: str = "utf-8"
    section_marker: str = "Sequences producing significant alignments:"
    header_keyword: str = "Description"
    end_patterns: Tuple[str, ...] = ("Max Score", "Total Score", "Query Cover", "E value", "Per.")
    top_hit_only: bool = True


@dataclass(frozen=True, slots=True)
class BlastTopHit:
    query_index: int
    description: str


class BlastTopHitExtractor:
    """
    Extracts top hit descriptions per query from BLAST text output.

    Note: not thread-safe due to stateful _last_hits.
    """

    def __init__(self, config: Optional[BlastParseConfig] = None) -> None:
        self.config = config or BlastParseConfig()
        self._last_hits: List[BlastTopHit] = []

    @property
    def last_hits(self) -> List[BlastTopHit]:
        return list(self._last_hits)

    def parse(self, input_file: PathLike) -> List[BlastTopHit]:
        inp = Path(input_file)
        if not inp.exists():
            raise FileNotFoundError(f"BLAST output not found: {inp}")
        if inp.is_dir():
            raise IsADirectoryError(f"BLAST output is a directory: {inp}")

        hits: List[BlastTopHit] = []
        in_alignments_section = False
        found_first_hit = False
        desc_col_end: Optional[int] = None
        query_index = -1

        with inp.open("r", encoding=self.config.encoding) as fin:
            for raw in fin:
                line = raw.rstrip("\r\n")

                if self.config.section_marker in line:
                    in_alignments_section = True
                    found_first_hit = False
                    desc_col_end = None
                    query_index += 1
                    continue

                if not in_alignments_section:
                    continue

                if self.config.header_keyword in line and desc_col_end is None:
                    desc_col_end = self._infer_desc_col_end(line)
                    continue

                if line.startswith("---"):
                    continue

                if not line.strip():
                    if found_first_hit:
                        in_alignments_section = False
                    continue

                if desc_col_end is None:
                    continue

                if not found_first_hit:
                    desc = line[:desc_col_end].strip()
                    if desc:
                        hits.append(BlastTopHit(query_index=query_index, description=desc))
                        found_first_hit = True
                        if self.config.top_hit_only:
                            in_alignments_section = False

        if not hits:
            raise ValueError("No hit descriptions found in BLAST output.")

        self._last_hits = hits
        return list(hits)

    def write(self, output_file: PathLike, hits: Optional[Sequence[BlastTopHit]] = None) -> Path:
        out = Path(output_file)
        data = list(hits) if hits is not None else self._last_hits
        if not data:
            raise ValueError("No hits to write (provide hits or call parse() first).")

        with out.open("w", encoding=self.config.encoding) as fout:
            for h in data:
                fout.write(h.description + "\n")
        return out

    def parse_and_write(self, input_file: PathLike, output_file: PathLike) -> Path:
        hits = self.parse(input_file)
        return self.write(output_file, hits)

    def _infer_desc_col_end(self, header_line: str) -> Optional[int]:
        for p in self.config.end_patterns:
            idx = header_line.find(p)
            if idx != -1:
                return idx

        desc_start = header_line.find(self.config.header_keyword)
        if desc_start == -1:
            return None

        search_start = desc_start + len(self.config.header_keyword)
        for i in range(search_start, len(header_line) - 3):
            if header_line[i: i + 2] == "  " and header_line[i + 2: i + 4].strip():
                return i
        return None


@dataclass(frozen=True, slots=True)
class GenBankNeighborsConfig:
    encoding: str = "utf-8"
    output_fasta: str = "selected_genes.fasta"
    n_before: int = 1
    n_after: int = 1
    feature_type: str = "CDS"
    translation_qualifier: str = "translation"
    gene_qualifier: str = "gene"
    locus_tag_qualifier: str = "locus_tag"
    product_qualifier: str = "product"


@dataclass(frozen=True, slots=True)
class CdsProtein:
    gene: Optional[str]
    locus_tag: Optional[str]
    product: Optional[str]
    translation: str

    def fasta_header(self) -> str:
        name = self.gene or self.locus_tag or "unknown_gene"
        prod = self.product or "unknown product"
        return f">{name} {prod}"


class GenBankNeighborsExtractor:
    """
    Extract neighboring CDS translations around target genes/locus_tags from GenBank.

    Requires Biopython. Neighbors are strictly the CDS entries adjacent to the
    target gene (the target itself is not included in the output).
    """

    def __init__(self, config: Optional[GenBankNeighborsConfig] = None) -> None:
        self.config = config or GenBankNeighborsConfig()
        self._last_all: List[CdsProtein] = []
        self._last_selected_indices: List[int] = []

        if self.config.n_before < 1 or self.config.n_after < 1:
            raise ValueError("n_before and n_after must be >= 1")

    @property
    def last_all_genes(self) -> List[CdsProtein]:
        return list(self._last_all)

    @property
    def last_selected_indices(self) -> List[int]:
        return list(self._last_selected_indices)

    def extract(
        self,
        input_gbk: PathLike,
        genes: Union[str, Sequence[str]],
        output_fasta: Optional[PathLike] = None,
        n_before: Optional[int] = None,
        n_after: Optional[int] = None,
    ) -> Path:
        inp = Path(input_gbk)
        if not inp.exists():
            raise FileNotFoundError(f"GenBank file not found: {inp}")
        if inp.is_dir():
            raise IsADirectoryError(f"GenBank path is a directory: {inp}")

        genes_set: Set[str] = {genes} if isinstance(genes, str) else set(genes)

        nb = n_before if n_before is not None else self.config.n_before
        na = n_after if n_after is not None else self.config.n_after
        if nb < 1 or na < 1:
            raise ValueError("n_before and n_after must be >= 1")

        all_genes = self._parse_genbank_cds(inp)
        if not all_genes:
            raise ValueError("No CDS translations found in GenBank file.")

        targets = self._find_target_indices(all_genes, genes_set)
        if not targets:
            raise ValueError(f"No target genes/locus_tags found: {sorted(genes_set)}")

        neighbor_indices = self._neighbor_indices(targets, len(all_genes), nb, na)
        out = Path(output_fasta) if output_fasta is not None else Path(self.config.output_fasta)

        self._write_proteins(out, [all_genes[i] for i in sorted(neighbor_indices)])

        self._last_all = all_genes
        self._last_selected_indices = sorted(neighbor_indices)

        return out

    def _parse_genbank_cds(self, gbk_path: Path) -> List[CdsProtein]:
        try:
            from Bio import SeqIO
        except ImportError as e:
            raise ImportError(
                "Biopython is required for GenBank parsing. "
                "Install it with: pip install biopython"
            ) from e

        proteins: List[CdsProtein] = []
        for record in SeqIO.parse(gbk_path, "genbank"):
            for feat in record.features:
                if feat.type != self.config.feature_type:
                    continue
                q: Dict[str, List[str]] = feat.qualifiers

                translation = self._first_qual(q, self.config.translation_qualifier)
                if not translation:
                    continue

                proteins.append(
                    CdsProtein(
                        gene=self._first_qual(q, self.config.gene_qualifier),
                        locus_tag=self._first_qual(q, self.config.locus_tag_qualifier),
                        product=self._first_qual(q, self.config.product_qualifier),
                        translation=translation,
                    )
                )
        return proteins

    @staticmethod
    def _first_qual(qualifiers: Dict[str, List[str]], key: str) -> Optional[str]:
        val = qualifiers.get(key)
        if not val:
            return None
        return val[0] if isinstance(val, list) else str(val)

    @staticmethod
    def _find_target_indices(all_genes: Sequence[CdsProtein], genes_set: Set[str]) -> List[int]:
        return [
            i for i, g in enumerate(all_genes)
            if g.gene in genes_set or g.locus_tag in genes_set
        ]

    @staticmethod
    def _neighbor_indices(
        targets: Sequence[int], n_total: int, n_before: int, n_after: int
    ) -> Set[int]:
        """Return indices of neighbors (target gene itself is excluded)."""
        selected: Set[int] = set()
        for idx in targets:
            for i in range(max(0, idx - n_before), idx):
                selected.add(i)
            for i in range(idx + 1, min(n_total, idx + 1 + n_after)):
                selected.add(i)
        return selected

    def _write_proteins(self, output_fasta: Path, proteins: Sequence[CdsProtein]) -> None:
        with output_fasta.open("w", encoding=self.config.encoding) as fout:
            for p in proteins:
                fout.write(p.fasta_header() + "\n")
                fout.write(p.translation + "\n")


