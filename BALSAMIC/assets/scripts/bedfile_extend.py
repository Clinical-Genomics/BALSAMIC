import click


@click.command()
@click.argument('input_bedfile', type=click.Path(exists=True))
@click.argument('output_bedfile', type=click.Path())
@click.option('--min_region_size', default=20, help='Minimum region size to enforce.')
def process_bedfile(input_bedfile: str, output_bedfile: str, min_region_size: int):
    """
    Process a BED file to ensure regions are at least a minimum size.

    Args:
        input_bedfile (str): Input BED file path.
        output_bedfile (str): Output BED file path.
        min_region_size (int): Minimum region size to enforce.
    """
    with open(input_bedfile, 'r') as infile, open(output_bedfile, 'w') as outfile:
        for line in infile:
            fields = line.strip().split('\t')
            if len(fields) != 3:
                continue  # Skip lines that don't have exactly 3 columns

            chrom = fields[0]
            start = int(fields[1])
            end = int(fields[2])

            region_length = end - start
            if region_length < min_region_size:
                center = (start + end) // 2
                half_size = min_region_size // 2
                start = max(0, center - half_size)
                end = center + half_size
                if min_region_size % 2 != 0:
                    end += 1

            outfile.write(f"{chrom}\t{start}\t{end}\n")


if __name__ == '__main__':
    process_bedfile()