from run_dna_rna_tools_module import gc_content


def is_good_gc_bounds(seq: str, gc_bound: tuple | float) -> bool:
    gc_lower = 0
    gc_upper = gc_bound
    if isinstance(gc_bound, tuple):
        gc_lower, gc_upper = gc_bound
    gc_bounds_seq = gc_content(seq)
    return gc_lower <= gc_bounds_seq <= gc_upper

def is_good_length_bounds(seq: str, length_bounds: tuple | float) -> bool:
    """Возвращает булево значение """
    length_lower = 0
    length_upper = length_bounds
    if isinstance(length_bounds, tuple):
        length_lower, length_upper = length_bounds
    return length_lower <= len(seq) <= length_upper


def is_good_quality_threshold(quality: str, quality_threshold: float) -> bool:
    if not quality:
        raise ValueError("Quality is empty string")
    quality_transform = [ord(quality)-33 for quality in quality.upper()]
    return sum(quality_transform) / len(quality_transform) >= quality_threshold

