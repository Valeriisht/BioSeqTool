import run_dna_rna_tools_module as tools

def is_good_gc_bounds(seq:str, gc_bound:tuple=(0, 100)) -> bool:
    gc_bounds_seq = tools.gc_content(seq)
    return gc_bounds_seq in gc_bound

def is_good_length_bounds(seq:str, length_bounds:float=0.2**32) -> bool:
    """Возвращает булевое значение """
    return len(seq) < length_bounds

def is_good_quality_threshold(quality_seq:str, quality_threshold:float=0) -> bool:
    quality_transform = [ord(quality) for quality in quality_seq.upper()]
    return sum(quality_transform)/len(quality_transform) > quality_threshold