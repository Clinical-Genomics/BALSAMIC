
import click
import pysam

@click.command()
@click.argument('input_bam', type=click.Path(exists=True))
@click.argument('output_bam', type=click.Path())
@click.option('--max_quality', default=70, type=int, help='Maximum quality value to cap to.')
def cap_base_qualities(input_bam: str, output_bam: str, max_quality: int):
    """
    Cap the base qualities in a BAM file.

    Args:
        input_bam (str): Input BAM file path.
        output_bam (str): Output BAM file path.
        max_quality (int): Maximum quality value to cap to.
    """
    # Open input BAM file for reading
    with pysam.AlignmentFile(input_bam, "rb") as in_bam:
        # Open output BAM file for writing
        with pysam.AlignmentFile(output_bam, "wb", header=in_bam.header) as out_bam:
            # Iterate over reads in the input BAM file
            for read in in_bam:
                # Cap base qualities to the maximum value
                capped_qualities = [min(q, max_quality) for q in read.query_qualities]
                # Update the base qualities in the read
                read.query_qualities = capped_qualities
                # Write the modified read to the output BAM file
                out_bam.write(read)

if __name__ == '__main__':
    cap_base_qualities()
