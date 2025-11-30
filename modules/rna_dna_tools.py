from typing import Any, List, Union

VALID = set("ACGTUacgtu")
DNA_COMP = str.maketrans("ATGCatgc", "TACGtacg")
RNA_COMP = str.maketrans("AUGCaugc", "UACGuacg")
DNA2RNA = str.maketrans("Tt", "Uu")


def _kind(seq: str) -> str:
    """Heuristic: returns 'RNA' if sequence contains U/u, otherwise 'DNA'."""
    return "RNA" if any(c in "Uu" for c in seq) else "DNA"


def _ok(seq: str) -> bool:
    """
    Quick validity check:
    - non-empty string
    - only A,C,G,T,U (any case)
    - no mixing T and U together
    """
    if not isinstance(seq, str) or not seq:
        return False
    s = set(seq)
    if not s <= VALID:
        return False
    return not (s & set("Tt") and s & set("Uu"))


def _check(seq: str) -> None:
    """
    Strict validator. Raises ValueError on:
    - empty / non-string
    - characters outside ACGTU
    - mixing T and U in the same sequence
    """
    if not isinstance(seq, str) or not seq:
        raise ValueError("Sequence must be a non-empty string.")
    bad_pos = [i for i, c in enumerate(seq) if c not in VALID]
    if bad_pos:
        bad = sorted({seq[i] for i in bad_pos})
        raise ValueError(
            f"Invalid chars {bad} at positions {bad_pos}; allowed A,C,G,T,U."
        )
    s = set(seq)
    if s & set("Tt") and s & set("Uu"):
        raise ValueError("T and U cannot be mixed in one sequence.")


def _reverse(s: str) -> str:
    """Reverse string."""
    return s[::-1]


def _complement(s: str) -> str:
    """DNA or RNA complement (auto-detected by presence of U/u)."""
    table = RNA_COMP if _kind(s) == "RNA" else DNA_COMP
    return s.translate(table)


def _revcomp(s: str) -> str:
    """Reverse complement for DNA or RNA."""
    return _reverse(_complement(s))


def _transcribe(s: str) -> str:
    """DNA -> RNA (replace T/t with U/u)."""
    return s.translate(DNA2RNA)


OPS = {
    "transcribe": _transcribe,
    "reverse": _reverse,
    "complement": _complement,
    "reverse_complement": _revcomp,
}


def _apply(seq: str, proc: str) -> Union[str, bool]:
    """
    Apply a single operation to one sequence.
    Returns:
      - bool for proc == 'is_nucleic_acid'
      - str for all other recognized operations
    """
    if proc == "is_nucleic_acid":
        return _ok(seq)
    _check(seq)
    try:
        fn = OPS[proc]
    except KeyError:
        avail = ", ".join(["is_nucleic_acid"] + sorted(OPS))
        raise ValueError(f"Unknown procedure {proc!r}. Available: {avail}.")
    return fn(seq)


def run_dna_rna_tools(*args: str, **kwargs: Any) -> Any:
    """
    Dispatcher for sequence utilities.

    Call styles:
      run_dna_rna_tools("ATGC", "reverse")
      run_dna_rna_tools("ATGC", "AUGC", proc="reverse_complement")

    Behavior:
      - If 'proc' is given as kwarg → all *args are treated as sequences.
      - Else → last positional argument is the procedure name,
               all previous positional args are sequences.
      - Returns a single value for one input sequence, or a list for multiple.

    Raises:
      ValueError/TypeError on invalid usage; ValueError on unknown procedure.
    """
    proc = kwargs.pop("proc", None)
    if proc is not None:
        seqs = list(args)
    else:
        if len(args) < 2:
            raise ValueError("Provide sequences and a procedure name.")
        *seqs, proc = args
    if kwargs:
        unknown = ", ".join(kwargs.keys())
        raise TypeError(f"Unknown keyword arguments: {unknown}")
    if not isinstance(proc, str):
        raise TypeError("Procedure name must be a string.")
    if not all(isinstance(s, str) for s in seqs):
        raise TypeError("All sequences must be strings.")

    out: List[Any] = [_apply(s, proc) for s in seqs]
    return out[0] if len(out) == 1 else out
