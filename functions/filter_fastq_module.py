import run_dna_rna_tools_module as tools
def is_good_gc_bounds(seq:str, gc_bound:tuple=(0, 100)) -> bool:
    gc_bounds_seq = tools.gc_content(seq)
    return gc_bounds_seq in gc_bound

def is_length_bounds(seq:str, length_bounds:float=0.2**32) -> bool:
    return len(seq) < length_bounds