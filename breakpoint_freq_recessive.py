"""
Python script to determine the frequency of breakpoints in the TAMU region.
This is fake data!

OUTPUTS:
- An excel file with 2 sheets
    - Deletions sheet characterizing breakpoint frequencies for deletions
    - Duplications sheet characterizing breakpoint frequencies for duplications
- Two .png files for:
    - a bar graph for deletion frequencies
    - a bar graph for duplication frequencies
"""

#!usr/bin/python

import csv
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import xlsxwriter

"""
Read the data and assemble subsets for future use.
"""

# Read csv data into a dataframe
all_data = pd.read_csv("TAMU_data.csv")
all_data = pd.DataFrame(all_data)

# Drop unlabeled column as index is present
all_data = all_data.drop(columns=["Unnamed: 0"])

# Create 2 data subsets: gene_probes, control_probes
control_probe = all_data.loc[:, "non_TAMU_probe_0":"non_TAMU_probe_49"]
gene_probe = all_data.loc[:, "TAMU_probe_0":"TAMU_probe_49"]

# Create 2 maintenance subsets: probe_header, ethnicity
probe_header = list(gene_probe.columns.values)
ethnicity_data = all_data.iloc[:, 0]

# Calculate number of samples for each ethnicity
total_A = len(all_data[all_data.ethnicity == "A"])
total_B = len(all_data[all_data.ethnicity == "B"])
total_C = len(all_data[all_data.ethnicity == "C"])

"""
Read depth is linearly proportional to copy number.

Assumption: CN = ~1 is a deletion.
Assumption: CN = ~3 is a duplication.

Assumption: There is variability in quality of sample, variability in lab procedures, and variability in the efficiency of each probe relative to each other.
Assumption: Based on above assumption, calculate copy number for each sample per TAMU_probe by comparison to its respective non_TAMU_probe.

Calculate copy number for each sample per probe.
"""

# Control probe: CN = 2x, where x = copy number of 1
# Gene probe: CN = nx, n needs to be calculated

# Convert to np array for ease of calculation
array_gene_probe = gene_probe.to_numpy()
array_control_probe = control_probe.to_numpy()

# CNsg = Read(sg) / Read(sc) * 2, where s = sample, g = gene_probe, c = control_probe
# Method: Conservative: +- 0.05


# Set thresholds for data to tabulate constitutive 'good' behavior
threshold_min_del = 0  # Min CN
threshold_max_del = 1.05
threshold_min_dup = 2.95
threshold_max_dup = 22.6  # Max CN

# Create Copy Number Array
cn_data = array_gene_probe / array_control_probe * 2

"""
Create a Boolean dataset based on threshold.
"""

# Create a boolean deletion array
bool_del = (cn_data > threshold_min_del) & (cn_data < threshold_max_del)

# Create a boolean duplication array
bool_dup = (cn_data > threshold_min_dup) & (cn_data < threshold_max_dup)

# convert boolean array back to dataframe
bool_del_df = pd.DataFrame(bool_del, columns=probe_header)
bool_dup_df = pd.DataFrame(bool_dup, columns=probe_header)

# Transpose boolean array
bool_del_dfT = bool_del_df.T
bool_dup_dfT = bool_dup_df.T

"""
Create a function with 2 parameters: dataframe and a minimum contguous value.
This value has been predetermined and assumed to be 4.
"""


def contiguous(df, min_contiguous=4):
    out = []
    for sample in df.columns:  # Iterate over columns (sample)
        s = df[sample].eq(True)
        hold = (df[sample] != df[sample].shift()).cumsum()[s]
        mask = hold.groupby(hold).transform("count").ge(min_contiguous)
        filter = hold[mask].reset_index()
        output = filter.groupby(sample)["index"].agg(["first", "last"])
        output.insert(0, "samp", sample)
        out.append(output)
    return pd.concat(out, ignore_index=True)


"""
Organize ethnicity data to fetch in the future.
"""

# Rename index on ethnicity data for future merging
ethnicity_data = ethnicity_data.rename_axis("samp").reset_index()


"""
Characterize Deletions.
"""

# Run function on boolean deletion data
del_probe_locs = contiguous(bool_del_dfT)

# Create merged df with ethnicity data
merged_del = pd.merge(del_probe_locs, ethnicity_data, how="left", on="samp")

# Organize dataframe by ethnicity and probes
deletion_frequency_df = (
    merged_del.groupby(["ethnicity", "first", "last"])
    .size()
    .reset_index()
    .rename(columns={0: "count"})
)

# Set conditions for categorizing data by ethnicity
conditions = [
    (deletion_frequency_df["ethnicity"] == "A"),
    (deletion_frequency_df["ethnicity"] == "B"),
    (deletion_frequency_df["ethnicity"] == "C"),
]

# Store number of samples per ethnicity (values = [4988, 2543, 2469])
values = [total_A, total_B, total_C]

# Create new column for total population per ethnicity
deletion_frequency_df["total_pop"] = np.select(conditions, values)

# Create a new column that joins the value of 5'breakpoint and 3'breakpoint
deletion_frequency_df["breakpoint"] = deletion_frequency_df[["first", "last"]].agg(
    "-".join, axis=1
)

# Get percentage of population per ethnicity per breakpoint for deletions
deletion_frequency_df = deletion_frequency_df.assign(
    pop_percent=lambda a: a["count"] / a["total_pop"] * 100
)

deletion_frequency_df = deletion_frequency_df.rename(
    columns={
        "first": "5' Breakpoint",
        "last": "3' Breakpoint",
        "ethnicity": "Ethnicity",
        "pop_percent": "Percent Population of Ethnicity",
        "breakpoint": "Breakpoint Region",
    }
)

final_deletion_output = deletion_frequency_df[
    ["Ethnicity", "5' Breakpoint", "3' Breakpoint", "Percent Population of Ethnicity"]
]

# Print deletion data
print("Deletion Frequency is as follows:\n", final_deletion_output)


"""
Characterize Duplications.
"""

# Run function on boolean duplication data
dup_probe_locs = contiguous(bool_dup_dfT)

# Create merged df with ethnicity data
merged_dup = pd.merge(dup_probe_locs, ethnicity_data, how="left", on="samp")

# Organize dataframe by ethnicity and probes
duplication_frequency_df = (
    merged_dup.groupby(["ethnicity", "first", "last"])
    .size()
    .reset_index()
    .rename(columns={0: "count"})
)

# Set conditions for categorizing data by ethnicity
conditions = [
    (duplication_frequency_df["ethnicity"] == "A"),
    (duplication_frequency_df["ethnicity"] == "B"),
    (duplication_frequency_df["ethnicity"] == "C"),
]

# Store number of samples per ethnicity (values = [4988, 2543, 2469])
values = [total_A, total_B, total_C]

# Create new column for total population per ethnicity
duplication_frequency_df["total_pop"] = np.select(conditions, values)

# Create a new column that joins the value of 5'breakpoint and 3'breakpoint
duplication_frequency_df["breakpoint"] = duplication_frequency_df[
    ["first", "last"]
].agg("-".join, axis=1)

# Get percentage of population per ethnicity per breakpoint for duplications
duplication_frequency_df = duplication_frequency_df.assign(
    pop_percent=lambda duplication_frequency_df: duplication_frequency_df["count"]
    / duplication_frequency_df["total_pop"]
    * 100
)

duplication_frequency_df = duplication_frequency_df.rename(
    columns={
        "first": "5' Breakpoint",
        "last": "3' Breakpoint",
        "ethnicity": "Ethnicity",
        "pop_percent": "Percent Population of Ethnicity",
        "breakpoint": "Breakpoint Region",
    }
)

final_duplication_output = duplication_frequency_df[
    ["Ethnicity", "5' Breakpoint", "3' Breakpoint", "Percent Population of Ethnicity"]
]

# Print duplication data
print("Duplication Frequency is as follows:\n", final_duplication_output)


"""
Plot Deletions.
"""

plot_del_df = deletion_frequency_df[
    ["Ethnicity", "Breakpoint Region", "Percent Population of Ethnicity"]
]
cat_plot_del = sns.catplot(
    x="Ethnicity",
    y="Percent Population of Ethnicity",
    hue="Breakpoint Region",
    data=plot_del_df,
    kind="bar",
).set(title="Deletion Frequency")


"""
Plot Duplications.
"""

plot_dup_df = duplication_frequency_df[
    ["Ethnicity", "Breakpoint Region", "Percent Population of Ethnicity"]
]
cat_plot_dup = sns.catplot(
    x="Ethnicity",
    y="Percent Population of Ethnicity",
    hue="Breakpoint Region",
    data=plot_dup_df,
    kind="bar",
).set(title="Duplication Frequency")


"""
Generate outs.
"""

# Write output to an excel file
writer = pd.ExcelWriter("Frequency.xlsx", engine="xlsxwriter")
final_deletion_output.to_excel(
    writer, sheet_name="Deletions", startrow=0, header=True, index=True
)
final_duplication_output.to_excel(
    writer, sheet_name="Duplications", startrow=0, header=True, index=True
)
writer.save()

# Write plots to .png files
cat_plot_del.figure.savefig("del_freq.png")
cat_plot_dup.figure.savefig("dup_freq.png")
