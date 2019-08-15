from io import StringIO
import tempfile
import pysam
from sketchy.utils import *


def create_tmp_fastq(contents: str) -> pysam.FastxFile:
    with tempfile.NamedTemporaryFile(mode="r+") as tmp:
        tmp.write(contents)
        tmp.truncate()
        return pysam.FastxFile(tmp.name)


def test_filterFastq_emptyFileWritesEmpty():
    contents = ""
    fastq = create_tmp_fastq(contents)
    output = StringIO()
    records = {"foo"}

    filter_fastq(fastq, output, records)
    output.seek(0)
    actual = output.read()
    expected = ""

    assert actual == expected


def test_filterFastq_invalidFormatWritesEmpty():
    contents = "I am not a correctly formatted file!"
    fastq = create_tmp_fastq(contents)
    output = StringIO()
    records = {"foo"}

    filter_fastq(fastq, output, records)
    output.seek(0)
    actual = output.read()
    expected = ""

    assert actual == expected


def test_filterFastq_noRequestedRecordsWritesEmpty():
    contents = "@read1\nATGC\n+\n$$$$\n"
    fastq = create_tmp_fastq(contents)
    output = StringIO()
    records = {"read2"}

    filter_fastq(fastq, output, records)
    output.seek(0)
    actual = output.read()
    expected = ""

    assert actual == expected


def test_filterFastq_oneRequestedRecordsWritesOneRecord():
    contents = "@read1\nATGC\n+\n$$$$\n"
    fastq = create_tmp_fastq(contents)
    output = StringIO()
    records = {"read1"}

    filter_fastq(fastq, output, records)
    output.seek(0)
    actual = output.read()
    expected = contents

    assert actual == expected


def test_filterFastq_oneRequestedRecordOneNotRequestedWritesOnlyRequested():
    contents = "@read1\nATGC\n+\n$$$$\n@read2\nATCC\n+\n$!$$\n"
    fastq = create_tmp_fastq(contents)
    output = StringIO()
    records = {"read2"}

    filter_fastq(fastq, output, records)
    output.seek(0)
    actual = output.read()
    expected = "@read2\nATCC\n+\n$!$$\n"

    assert actual == expected


def test_filterFastq_twoRequestedRecordOneNotRequestedWritesOnlyRequested():
    contents = "@read1\nATGC\n+\n$$$$\n@read2\nATCC\n+\n$!$$\n@read3\nGGGG\n+\n$$$$\n"
    fastq = create_tmp_fastq(contents)
    output = StringIO()
    records = {"read1", "read3"}

    filter_fastq(fastq, output, records)
    output.seek(0)
    actual = output.read()
    expected = "@read1\nATGC\n+\n$$$$\n@read3\nGGGG\n+\n$$$$\n"

    assert actual == expected


def test_filterFastq_fastaIsAlsoParsed():
    contents = ">read1\nATGC\n>read2\nATCC\n>read3\nGGGG\n\n"
    fastq = create_tmp_fastq(contents)
    output = StringIO()
    records = {"read1", "read3"}

    filter_fastq(fastq, output, records)
    output.seek(0)
    actual = output.read()
    expected = ">read1\nATGC\n>read3\nGGGG\n"

    assert actual == expected
