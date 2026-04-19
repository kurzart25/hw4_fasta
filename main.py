from __future__ import annotations

from abc import ABC, abstractmethod
from collections import Counter
from dataclasses import dataclass
from typing import ClassVar, Dict, Union

Index = Union[int, slice]


@dataclass(frozen=True, slots=True)
class BiologicalSequence(ABC):
    _seq: str

    def __post_init__(self) -> None:
        if not isinstance(self._seq, str) or not self._seq:
            raise ValueError("Sequence must be a non-empty string.")
        object.__setattr__(self, "_seq", self._seq.upper())
        self.check_alphabet()

    def __len__(self) -> int:
        return len(self._seq)

    def __getitem__(self, item: Index) -> Union[str, "BiologicalSequence"]:
        if isinstance(item, slice):
            sliced = self._seq[item]
            if not sliced:
                raise ValueError("Slice resulted in an empty sequence.")
            return type(self)(sliced)
        return self._seq[item]

    def __str__(self) -> str:
        return f"{self.__class__.__name__}(len={len(self)}, seq='{self._seq}')"

    def __repr__(self) -> str:
        return f"{self.__class__.__name__}({self._seq!r})"

    @property
    def seq(self) -> str:
        return self._seq

    @abstractmethod
    def check_alphabet(self) -> None:
        """Raise ValueError if sequence contains invalid symbols."""
        ...


@dataclass(frozen=True, slots=True)
class NucleicAcidSequence(BiologicalSequence, ABC):
    alphabet: ClassVar[set[str]] = set()
    complement_map: ClassVar[Dict[str, str]] = {}

    def check_alphabet(self) -> None:
        if not self.alphabet:
            raise NotImplementedError(
                "NucleicAcidSequence is abstract; define alphabet in subclass."
            )
        bad = sorted({c for c in self._seq if c not in self.alphabet})
        if bad:
            raise ValueError(
                f"Invalid nucleotide symbols: {bad}. Allowed: {sorted(self.alphabet)}"
            )

    def complement(self) -> "NucleicAcidSequence":
        if not self.complement_map:
            raise NotImplementedError(
                "NucleicAcidSequence is abstract; define complement_map in subclass."
            )
        return type(self)("".join(self.complement_map[c] for c in self._seq))

    def reverse(self) -> "NucleicAcidSequence":
        return type(self)(self._seq[::-1])

    def reverse_complement(self) -> "NucleicAcidSequence":
        return self.complement().reverse()


@dataclass(frozen=True, slots=True)
class DNASequence(NucleicAcidSequence):
    alphabet: ClassVar[set[str]] = set("ACGTN")
    complement_map: ClassVar[Dict[str, str]] = {
        "A": "T",
        "T": "A",
        "G": "C",
        "C": "G",
        "N": "N",
    }

    def transcribe(self) -> "RNASequence":
        return RNASequence(self._seq.replace("T", "U"))


@dataclass(frozen=True, slots=True)
class RNASequence(NucleicAcidSequence):
    alphabet: ClassVar[set[str]] = set("ACGUN")
    complement_map: ClassVar[Dict[str, str]] = {
        "A": "U",
        "U": "A",
        "G": "C",
        "C": "G",
        "N": "N",
    }


@dataclass(frozen=True, slots=True)
class AminoAcidSequence(BiologicalSequence):
    alphabet: ClassVar[set[str]] = set("ACDEFGHIKLMNPQRSTVWY")

    def check_alphabet(self) -> None:
        bad = sorted({c for c in self._seq if c not in self.alphabet})
        if bad:
            raise ValueError(
                f"Invalid amino acid symbols: {bad}. Allowed: {sorted(self.alphabet)}"
            )

    def aa_composition(self) -> Dict[str, int]:
        return dict(Counter(self._seq))