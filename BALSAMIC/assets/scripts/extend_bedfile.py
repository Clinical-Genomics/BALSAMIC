import click


@click.command()
@click.argument("input_bedfile", type=click.Path(exists=True))
@click.argument("output_bedfile", type=click.Path())
@click.option("--extend-to-min-region-size", default=100, help="Will extend regions shorter than the specified size to this minimum size.")
def extend_bedfile(input_bedfile: str, output_bedfile: str, extend_to_min_region_size: int):
    """
    Process a BED file to ensure regions are at least a minimum size.

    Args:
        input_bedfile (str): Input BED file path.
        output_bedfile (str): Output BED file path.
        min_region_size (int): Minimum region size to enforce.
    """
    with open(input_bedfile, "r") as infile, open(output_bedfile, "w") as outfile:
        for line in infile:
            fields = line.strip().split("\t")

            chrom: str = fields[0]
            start = int(fields[1])
            end = int(fields[2])

            region_length: int = end - start
            if region_length < extend_to_min_region_size:
                center = (start + end) // 2
                half_size = extend_to_min_region_size // 2
                start = max(0, center - half_size)
                end = center + half_size
                if extend_to_min_region_size % 2 != 0:
                    end += 1

            outfile.write(f"{chrom}\t{start}\t{end}\n")


if __name__ == "__main__":
    extend_bedfile()
