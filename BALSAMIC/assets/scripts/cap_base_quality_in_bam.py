import click
import pysam
import numpy as np


@click.command()
@click.argument("input_bam", type=click.Path(exists=True))
@click.argument("output_bam", type=click.Path())
@click.option(
    "--max-quality",
    default=70,
    type=int,
    help="Maximum quality value to cap to.",
)
def cap_base_qualities(input_bam: str, output_bam: str, max_quality: int):
    """
    Cap the base qualities in a BAM file.

    Args:
        input_bam (str): Input BAM file path.
        output_bam (str): Output BAM file path.
        max_quality (int): Maximum quality value to cap to.
    """
    # Open input BAM file for reading
    samfile = pysam.AlignmentFile(input_bam, "rb")
    out_bam = pysam.AlignmentFile(output_bam, "wb", header=samfile.header)
    for read in samfile.fetch():
        qualities = np.array(read.query_qualities)
        capped_qualities = np.minimum(qualities, max_quality)
        # Update the base qualities in the read
        read.query_qualities = capped_qualities.tolist()
        # Write the modified read to the output BAM file
        out_bam.write(read)


if __name__ == "__main__":
    cap_base_qualities()
