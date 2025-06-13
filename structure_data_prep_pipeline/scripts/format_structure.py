import pandas as pd

file = snakemake.input[0]
output = snakemake.output[0]

df = pd.read_csv(file, sep='\t', lineterminator='\n', header=None)

positions_processed = 0
for col_index, col_data in df.iteritems():
    if col_index == 0:
        continue
    unique = sorted(list(set(col_data.values)))
    if len(unique) == 3:
        df[col_index] = df[col_index].replace([unique[0], unique[1], unique[2]], [-9, 1, 0])
    else:
        df[col_index] = df[col_index].replace([unique[0], unique[1]], [1, 0])
    positions_processed += 1

df.to_csv(output, sep='\t', index=False, header=None)
