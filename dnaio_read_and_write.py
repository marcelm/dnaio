#!/usr/bin/env python3

import dnaio
import sys

with dnaio.open(sys.argv[1]) as in_records:
    with dnaio.open(sys.argv[2], mode="w", fileformat="fastq") as out_records:
        for record in in_records:
            out_records.write(record)
