import pandas as pd


def main():
    def output_infection():
        with open("../data/mega-virus_aligned.out.sam", "r") as file:
            count = 0
            virus_species = []
            for line in file:
                if not line.startswith("@"):
                    virus_species.append(line.split("\t")[2])
                    count += 1
            total = len(virus_species)
            virus_species = list(set(virus_species))
            print(len(virus_species))

        species_name = []
        species_count = []
        species_ratio = []
        stopper = 0
        reference = pd.read_csv('../data/new_virus.species.txt', sep='\t', names=["id", "name"])
        reference = reference.set_index("id")
        for species in virus_species:
            count = 0
            with open("../data/mega-virus_aligned.out.sam", "r") as file:
                # species_count.append(subprocess.check_output("grep -v @ " + "../data/mega-virus_aligned.out.sam"
                #                                              + " | grep '" + species + "' | wc -l", shell=True))
                for line in file:
                    if (not line.startswith("@")) and species in line:
                        count += 1
            species_count.append(count)
            species_name.append(reference.loc[species, "name"])
            species_ratio.append(count / total * 100.0)

            # if stopper >= 5:
            #     exit
            stopper += 1

            print(species)
            print(species_name)
            print(species_count)
            print(species_ratio)
        output = pd.DataFrame({"Name": species_name,
                               "Count": species_count,
                               "Percentage": species_ratio})
        output.sort_values(by="Percentage", ascending=False, inplace=True)
        output.to_csv(path_or_buf="../data/output.csv", index=False)
        # print(virus_species[:5])
        # print(species_name[:5])
        # print(species_count[:5])
        # print(len(species_count))
        # print(len(species_name))
        return output

    print(output_infection())


if __name__ == '__main__':
    main()
