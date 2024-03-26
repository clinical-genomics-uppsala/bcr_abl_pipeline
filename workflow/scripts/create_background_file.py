import gzip
import statistics
import sys

gvcf_files = snakemake.input.gvcfs
background_file = open(snakemake.output.background_file, "w")
min_dp = snakemake.params.min_dp
max_af = snakemake.params.max_af

background_dict = {}

for file_name in gvcf_files:
    print(file_name)
    with gzip.open(file_name, "rt") as infile:
        file_content = infile.read().split("\n")
        header = True
        for line in file_content:
            if header:
                if line[:6] == "#CHROM":
                    header = False
                continue
            columns = line.strip().split("\t")
            if len(columns) <= 1:
                continue
            if columns[6] == "PASS":
                chrom = columns[0][3:]
                pos = columns[1]
                key = chrom + "_" + pos
                format = columns[8].split(":")
                data = columns[9].split(":")
                AD_id = 0
                for f in format:
                    if f == "AD":
                        break
                    AD_id += 1
                AD_info = data[AD_id].split(",")
                if columns[4] == ".":
                    alt_AF = float(data[AD_id + 2])
                else:
                    # ref_AD = int(AD_info[0])
                    # import pdb; pdb.set_trace()
                    DP = int(data[AD_id + 1])
                    alt_AD = DP - sum([int(x) for x in AD_info])
                    for AD in AD_info[1:]:
                        alt_AD += int(AD)
                    alt_AF = 0.0
                    if DP > min_dp:
                        alt_AF = float(alt_AD / float(DP))
                    else:
                        continue
                if alt_AF > 1 - float(max_af):
                    alt_AF = float(1 - alt_AF)
                if alt_AF > float(max_af):
                    continue
                if key in background_dict:
                    background_dict[key].append(alt_AF)
                else:
                    background_dict[key] = [alt_AF]

background_file.write("Median\tSD\tNrSamples\n")
for key in background_dict:
    background_dict[key].sort()
    nr_obs = len(background_dict[key])
    if nr_obs >= 4:
        median_background = statistics.median(background_dict[key])
        """This is the sample variance s² with Bessel’s correction, also known as variance with N-1 degrees of freedom.
        Provided that the data points are representative (e.g. independent and identically distributed),
        the result should be an unbiased estimate of the true population variance."""
        stdev_background = statistics.stdev(background_dict[key])
        background_file.write(
            "chr"
            + key.split("_")[0]
            + "\t"
            + key.split("_")[1]
            + "\t"
            + str(median_background)
            + "\t"
            + str(stdev_background)
            + "\t"
            + str(nr_obs)
            + "\n"
        )

background_file.close()
