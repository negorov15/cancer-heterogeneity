import pytest
import pandas as pd
from main import get_transcript, get_exons, get_strand
from calculate_distance import calculate_distance

# Test data
VARIANTS = [100, 104, 109, 200, 204, 209, 300, 304, 309, 400, 605]
CHROMOSOME = "chr1"
ANNOTATION = pd.DataFrame({
    "chr_name": ["chr1"] * 8,
    "genomic start": [100, 200, 300, 600, 602, 603, 700, 800],
    "genomic end": [109, 209, 309, 609, 611, 612, 709, 809],
    "genomic strand": ["+", "-", "+", "+", "+", "+", "+", "+"],
    "transcript_id": ["TRANSCRIPT_1", "TRANSCRIPT_2", "TRANSCRIPT_3", "TRANSCRIPT_4", "TRANSCRIPT_5", "TRANSCRIPT_6", "TRANSCRIPT_4", "TRANSCRIPT_4"]
})


def test_get_transcript():
    assert get_transcript(VARIANTS[0], ANNOTATION) == ["TRANSCRIPT_1"]
    assert get_transcript(VARIANTS[1], ANNOTATION) == ["TRANSCRIPT_1"]
    assert get_transcript(VARIANTS[2], ANNOTATION) == ["TRANSCRIPT_1"]
    assert get_transcript(VARIANTS[3], ANNOTATION) == ["TRANSCRIPT_2"]
    assert get_transcript(VARIANTS[4], ANNOTATION) == ["TRANSCRIPT_2"]
    assert get_transcript(VARIANTS[5], ANNOTATION) == ["TRANSCRIPT_2"]
    assert get_transcript(VARIANTS[6], ANNOTATION) == ["TRANSCRIPT_3"]
    assert get_transcript(VARIANTS[7], ANNOTATION) == ["TRANSCRIPT_3"]
    assert get_transcript(VARIANTS[8], ANNOTATION) == ["TRANSCRIPT_3"]
    assert get_transcript(VARIANTS[9], ANNOTATION) == None
    assert get_transcript(VARIANTS[10], ANNOTATION) == ["TRANSCRIPT_4", "TRANSCRIPT_5", "TRANSCRIPT_6"]
    with pytest.raises(TypeError):
        get_transcript("200", ANNOTATION)


def test_get_exons():
    assert get_exons("TRANSCRIPT_1", ANNOTATION) == [(100, 109)]
    assert get_exons("TRANSCRIPT_2", ANNOTATION) == [(200, 209)]
    assert get_exons("TRANSCRIPT_3", ANNOTATION) == [(300, 309)]
    assert get_exons("TRANSCRIPT_4", ANNOTATION) == [(600, 609), (700, 709), (800, 809)]
    assert get_exons("TRANSCRIPT_10", ANNOTATION) == None
    with pytest.raises(TypeError):
            get_exons(200, ANNOTATION)


def test_get_strand():
    assert get_strand("TRANSCRIPT_1", ANNOTATION) == "+"
    assert get_strand("TRANSCRIPT_2", ANNOTATION) == "-"
    assert get_strand("TRANSCRIPT_3", ANNOTATION) == "+"
    assert get_strand("TRANSCRIPT_4", ANNOTATION) == "+"
    assert get_strand("TRANSCRIPT_10", ANNOTATION) == None
    with pytest.raises(TypeError):
                get_strand(200, ANNOTATION)


EXONS = [(100, 109), (200, 209), (300, 309)]
EXONS_UNSORTED = [(300, 309), (200, 209), (100, 109)]

def test_calculate_distance():
    # Last exon distance (E3): sorted
    assert calculate_distance(300, EXONS, "+") == 9
    assert calculate_distance(309, EXONS, "+") == 0
    assert calculate_distance(300, EXONS, "-") == 20
    assert calculate_distance(309, EXONS, "-") == 29
    assert calculate_distance(304, EXONS, "+") == 5
    assert calculate_distance(304, EXONS, "-") == 24
    # Last exon distance (E3): unsorted
    assert calculate_distance(300, EXONS_UNSORTED, "+") == 9
    assert calculate_distance(309, EXONS_UNSORTED, "+") == 0
    assert calculate_distance(300, EXONS_UNSORTED, "-") == 20
    assert calculate_distance(309, EXONS_UNSORTED, "-") == 29
    assert calculate_distance(304, EXONS_UNSORTED, "+") == 5
    assert calculate_distance(304, EXONS_UNSORTED, "-") == 24
    # Middle exon distance (E2): sorted
    assert calculate_distance(200, EXONS, "+") == 19
    assert calculate_distance(209, EXONS, "+") == 10
    assert calculate_distance(200, EXONS, "-") == 10
    assert calculate_distance(209, EXONS, "-") == 19   
    assert calculate_distance(204, EXONS, "+") == 15
    assert calculate_distance(204, EXONS, "-") == 14
    # Middle exon distance (E2): unsorted
    assert calculate_distance(200, EXONS_UNSORTED, "+") == 19
    assert calculate_distance(209, EXONS_UNSORTED, "+") == 10
    assert calculate_distance(200, EXONS_UNSORTED, "-") == 10
    assert calculate_distance(209, EXONS_UNSORTED, "-") == 19   
    assert calculate_distance(204, EXONS_UNSORTED, "+") == 15
    assert calculate_distance(204, EXONS_UNSORTED, "-") == 14
    # First exon distance (E1): sorted
    assert calculate_distance(100, EXONS, "+") == 29
    assert calculate_distance(109, EXONS, "+") == 20
    assert calculate_distance(100, EXONS, "-") == 0
    assert calculate_distance(109, EXONS, "-") == 9
    assert calculate_distance(104, EXONS, "+") == 25
    assert calculate_distance(104, EXONS, "-") == 4
    # First exon distance (E1): unsorted
    assert calculate_distance(100, EXONS_UNSORTED, "+") == 29
    assert calculate_distance(109, EXONS_UNSORTED, "+") == 20
    assert calculate_distance(100, EXONS_UNSORTED, "-") == 0
    assert calculate_distance(109, EXONS_UNSORTED, "-") == 9
    assert calculate_distance(104, EXONS_UNSORTED, "+") == 25
    assert calculate_distance(104, EXONS_UNSORTED, "-") == 4
    # Out of bounds distance
    assert calculate_distance(400, EXONS, "+") == None
    assert calculate_distance(400, EXONS, "-") == None
    # Invalid strand
    with pytest.raises(ValueError):
        calculate_distance(200, EXONS, "*")