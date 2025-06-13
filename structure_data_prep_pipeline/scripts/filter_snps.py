input_file = snakemake.input[0]
output_file = snakemake.output[0]

with open(input_file, "r") as f, open(output_file, "w") as w:
    count = 0
    for line in f:
        count += 1
        fields = line.strip().split()
        if count == 1:
            w.write("\t".join(fields) + "\n")
        elif len(fields[1]) == 1 and len(fields[2]) == 1:
            w.write("\t".join(fields) + "\n")
