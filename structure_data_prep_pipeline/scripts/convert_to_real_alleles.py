input_file = snakemake.input[0]
output_file = snakemake.output[0]

with open(input_file, "r") as f, open(output_file, "w") as w:
    for line in f:
        fields = line.strip().split()
        ref, alt = fields[1], fields[2]
        converted = [
            ref if val == "A" else alt if val == "B" else val
            for val in fields
        ]
        w.write("\t".join(converted) + "\n")
