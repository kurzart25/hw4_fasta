import sys, os
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), 'modules')))

from filter_fastq import filter_fastq, filter_fastq_stream

input_file = "example_data/example_fastq.fastq"

# In-memory 
print(">>> Testing filter_fastq (in-memory)")
out1 = filter_fastq(
    input_fastq=input_file,
    output_fastq="filtered_memory.fq",
    gc_bounds=(40, 60),
    length_bounds=(50, 500),
    quality_threshold=20.0,
    on_duplicate="skip",
    dup_dir="duplicates"
)
print("Output (in-memory):", out1)

#Streaming  
print("\n>>> Testing filter_fastq_stream (stream)")
out2 = filter_fastq_stream(
    input_fastq=input_file,
    output_fastq="filtered_stream.fq",
    gc_bounds=(40, 60),
    length_bounds=(50, 500),
    quality_threshold=20.0,
)
print("Output (stream):", out2)

def count_reads(path):
    n = 0
    with open(path, "r", encoding="utf-8") as fh:
        for line in fh:
            if line.startswith("@"):
                n += 1
    return n

print("\nKept reads (in-memory):", count_reads(out1))
print("Kept reads (stream):   ", count_reads(out2))
