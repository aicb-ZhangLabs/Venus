import pandas as pd


def main():
    """Converts genomic coordinates into genes for final output"""
    with open("../data/guide-Aligned.out.sam", "r") as file:
        quality_reads = []
        for line in file:
            if not line.startswith("@"):
                quality_reads.append(line.split("\t")[0])
                # print(line.split("\t")[0])
        quality_reads = list(set(quality_reads))

    integration_sites = pd.read_csv("../data/Chimeric.out.junction", sep="\t")
    integrA = integration_sites[(integration_sites["read_name"].isin(quality_reads)) &
                            (integration_sites["chr_donorA"] != "NC_001802.1")][["chr_donorA",
                                                                                 "brkpt_donorA",
                                                                                 "read_name"]]

    integrA.rename(columns={'chr_donorA': 'chr', 'brkpt_donorA': 'coord'}, inplace=True)
    integrB = integration_sites[(integration_sites["read_name"].isin(quality_reads)) &
                            (integration_sites["chr_acceptorB"] != "NC_001802.1")][["chr_acceptorB",
                                                                                    "brkpt_acceptorB",
                                                                                    "read_name"]]
    integrB.rename(columns={'chr_acceptorB': 'chr', 'brkpt_acceptorB': 'coord'}, inplace=True)
    integration = pd.concat([integrA, integrB])

    integration.replace(["NC_000001.11", "NC_000002.12", "NC_000003.12", "NC_000004.12",
                         "NC_000005.10", "NC_000006.12", "NC_000007.14", "NC_000008.11",
                         "NC_000009.12", "NC_000010.11", "NC_000011.10", "NC_000012.12",
                         "NC_000013.11", "NC_000014.9", "NC_000015.10", "NC_000016.10",
                         "NC_000017.11", "NC_000018.10", "NC_000019.10", "NC_000020.11",
                         "NC_000021.9", "NC_000022.11", "NC_000023.11", "NC_000024.10"],
                        ["chr1", "chr2", "chr3", "chr4",
                         "chr5", "chr6", "chr7", "chr8",
                         "chr9", "chr10", "chr11", "chr12",
                         "chr13", "chr14", "chr15", "chr16",
                         "chr17", "chr18", "chr19", "chr20",
                         "chr21", "chr22", "chrX", "chrY"], inplace=True)
    print(integration)

    return


# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    main()