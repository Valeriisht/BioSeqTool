import pytest
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqUtils import gc_fraction as GC
from io import StringIO
import tempfile

import sys
import os

sys.path.append(os.path.dirname(os.getcwd()))

from NASeqTool import filter_fastq


# ======================
# ======================


# enviromental settings
@pytest.fixture
def tmp_fastq_file():
    """Create a temporary FASTQ file with various test sequences"""
    records = [
        # Perfect record - passes all filters
        SeqRecord(
            Seq("ATGC"),  # GC=50%, length=4
            id="perfect_record",
            description="Passes all filters",
            letter_annotations={"phred_quality": [40, 40, 40, 40]},  # avg=40
        ),
        # Fails GC and quality
        SeqRecord(
            Seq("AAAAAA"),  # GC=0%, length=6
            id="low_gc_poor_quality",
            description="Fails GC and quality filters",
            letter_annotations={"phred_quality": [10, 10, 10, 10, 10, 10]},  # avg=10
        ),
        # Edge case GC
        SeqRecord(
            Seq("GGGA"),  # GC=75%, length=4
            id="high_gc_record",
            description="Edge case GC content",
            letter_annotations={"phred_quality": [30, 30, 30, 30]},  # avg=30
        ),
        # Fails length
        SeqRecord(
            Seq("ATGCATGC"),  # GC=50%, length=8
            id="long_sequence",
            description="Fails length filter",
            letter_annotations={
                "phred_quality": [20, 20, 20, 20, 20, 20, 20, 20]
            },  # avg=20
        ),
        # Additional test cases
        SeqRecord(
            Seq("GATGACA"),  # GC=42.8%, length=7
            id="medium_gc",
            description="Medium GC content",
            letter_annotations={
                "phred_quality": [25, 25, 25, 25, 25, 25, 25]
            },  # avg=25
        ),
    ]

    # Return a file-like object that can be used as a temporary storage area.
    with tempfile.NamedTemporaryFile(
        mode="w+t", buffering=-1, prefix=".fastq", delete=False
    ) as tmp:
        SeqIO.write(records, tmp, "fastq")

        # The flushing system basically makes it  possible to record a message
        # at the end of a request and access it next request and only next request.

        tmp.flush()
        yield tmp.name  # sleep
        os.unlink(tmp.name)


@pytest.fixture
def output_fastq_file():
    """Create a temporary output file path"""
    with tempfile.NamedTemporaryFile(
        mode="w+t", buffering=-1, prefix=".fastq", delete=False
    ) as tmp:
        tmp.close()
        yield tmp.name  # sleep
        os.unlink(tmp.name)


# 1. Launch with default params
def test_fastq_file(tmp_fastq_file, output_fastq_file):
    """Test with default filters"""

    filter_fastq(tmp_fastq_file, output_fastq_file)
    records = list(SeqIO.parse(output_fastq_file, "fastq"))
    assert len(records) == 5
    assert set([rec.id for rec in records]) == {
        "perfect_record",
        "low_gc_poor_quality",
        "high_gc_record",
        "long_sequence",
        "medium_gc",
    }


# 2. Launch with gc-content params
# ================ gc-content


@pytest.mark.parametrize(
    "gc_bounds, expect_ids",
    [
        (
            (0, 1),
            {
                "perfect_record",
                "low_gc_poor_quality",
                "high_gc_record",
                "long_sequence",
                "medium_gc",
            },
        ),
        ((0.7, 1), {"high_gc_record"}),
        ((0.40, 0.60), {"perfect_record", "long_sequence", "medium_gc"}),
        ((0, 0), {"low_gc_poor_quality"}),  # Ни одна запись не имеет ровно 0% GC
    ],
)
def test_gc_bounds(tmp_fastq_file, output_fastq_file, gc_bounds, expect_ids):
    filter_fastq(tmp_fastq_file, output_fastq_file, gc_bounds=gc_bounds)

    records = list(SeqIO.parse(output_fastq_file, "fastq"))
    assert set(rec.id for rec in records) == expect_ids


# 3. Launch with length params
# ================ length


@pytest.mark.parametrize(
    "length_bounds, expect_ids",
    [
        (
            (0, 100),
            {
                "perfect_record",
                "low_gc_poor_quality",
                "high_gc_record",
                "long_sequence",
                "medium_gc",
            },
        ),
        ((4, 4), {"perfect_record", "high_gc_record"}),
        ((6, 7), {"low_gc_poor_quality", "medium_gc"}),
        ((8, 100), {"long_sequence"}),
    ],
)
def test_length_bounds(tmp_fastq_file, output_fastq_file, length_bounds, expect_ids):
    filter_fastq(tmp_fastq_file, output_fastq_file, length_bounds=length_bounds)

    records = list(SeqIO.parse(output_fastq_file, "fastq"))
    assert set(rec.id for rec in records) == expect_ids


# 4. Launch with quality params
# ================ quality


@pytest.mark.parametrize(
    "quality_threshold, expect_ids",
    [
        (
            0,
            {
                "perfect_record",
                "low_gc_poor_quality",
                "high_gc_record",
                "long_sequence",
                "medium_gc",
            },
        ),
        (20, {"perfect_record", "high_gc_record", "long_sequence", "medium_gc"}),
        (30, {"perfect_record", "high_gc_record"}),
        (40, {"perfect_record"}),
    ],
)
def test_quality_threshold(
    tmp_fastq_file, output_fastq_file, quality_threshold, expect_ids
):
    filter_fastq(tmp_fastq_file, output_fastq_file, quality_threshold=quality_threshold)

    records = list(SeqIO.parse(output_fastq_file, "fastq"))
    assert set(rec.id for rec in records) == expect_ids


# 5. Test Error
def test_error_handling(tmp_fastq_file, output_fastq_file):
    """Test handling of non-existent input file"""

    with pytest.raises(SystemError, match="File does not exist"):
        filter_fastq("nonexistent.fastq", output_fastq_file)

    with pytest.raises(ValueError, match="GC bound must be between 0 and 100"):
        filter_fastq(tmp_fastq_file, output_fastq_file, gc_bounds=(-10, 20))

    with pytest.raises(ValueError, match="Length bounds are invalid"):
        filter_fastq(tmp_fastq_file, output_fastq_file, length_bounds=(100, 10))


# 6. Creation file
def test_output_file_creation(tmp_fastq_file, output_fastq_file):
    """Test checks the correct path"""
    filter_fastq(tmp_fastq_file, output_fastq_file)
    assert os.path.exists(output_fastq_file)
