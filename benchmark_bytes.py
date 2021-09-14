import timeit
import statistics


def benchmark(name, line, setup, loops=100_000, runs=10):
    print(f"{name}")
    results = [timeit.timeit(line, setup, number=loops) for _ in range(runs)]
    # Calculate microseconds
    results = [(result / loops) * 1_000_000 for result in results]
    print(f"average: {round(statistics.mean(results), 3)}, "
          f"range: {round(min(results), 3)}-{round(max(results),3)} "
          f"stdev: {round(statistics.stdev(results),3)}")


SETUP = """
from dnaio import Sequence

seq = Sequence("chr1_245318753_245318470_0_0_0_0_0:0:0_0:0:0_0/1",
               "CTGAAATTTCATGAGGCCAGAGTAGGTCAACAGGGGCTAATGAGACAAGATCTAACTAAACA"
               "GAGTAAAGAGTTTCTCGGTCACAAGGGAAACTTTTACTATTTGAAATCAGCCTGAGCCAAGA"
               "TTGATGAGGAAAAAAAAACAAAAACCAA",
               "I?>DC:>@?IDC9??G?>EH9E@66=9<?@E?DC:@<@BBFG>=FIC@F9>7CG?IC?I;CD"
               "9>>>A@C7>>8>>D9GCB<;?DD>C;9?>5G>?H?=6@>:G6B<?==A7?@???8IF<75C="
               "@A:BEA@A;C89D:=1?=<A>D=>B66C")
"""

if __name__ == "__main__":
    benchmark("Calculate fastq bytes", "seq.fastq_bytes()", SETUP)
    benchmark("Calculate fastq bytes 2", "seq.fastq_bytes2()", SETUP)
    benchmark("Calculate fastq bytes 3", "seq.fastq_bytes3()", SETUP)
    benchmark("Calculate fastq bytes 4", "seq.fastq_bytes4()", SETUP)
    benchmark("Calculate fastq bytes 5", "seq.fastq_bytes5()", SETUP)
    benchmark("Encode qualities ascii", "seq.qualities.encode('ascii')", SETUP)
    benchmark("Encode qualities ascii cython", "seq.qualities_bytes()", SETUP)
    benchmark("Encode qualities latin-1", "seq.qualities.encode('latin-1')",
              SETUP)

