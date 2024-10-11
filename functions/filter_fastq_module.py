from run_dna_rna_tools_module import gc_content


def is_good_gc_content(seq: str, gc_bound: tuple | float) -> bool:
    """Checks to ensure that the gc content of the reid fits
    thresholds."""
    gc_lower, gc_upper = gc_bound\
        if isinstance(gc_bound, tuple) else (0, gc_bound)
    gc_bounds_seq = gc_content(seq)
    return gc_lower <= gc_bounds_seq <= gc_upper


def is_good_length(seq: str, length_bounds: tuple | float) -> bool:
    """Checks to make sure the length of the reid
     is within the threshold values."""
    length_lower, length_upper = length_bounds\
        if isinstance(length_bounds, tuple) else (0, length_bounds)
    return length_lower <= len(seq) <= length_upper


def is_good_quality(quality: str, quality_threshold: float) -> bool | None:
    """Checks to ensure that the reid quality
     on the phred33 scale is within the threshold values."""
    if not quality:
        return
    quality_transform = [ord(quality) - 33 for quality in quality.upper()]
    return sum(quality_transform) / len(quality_transform) >= quality_threshold
