import os
from io import StringIO
from typing import List, Optional

from bam2sdf import bam2sdf
from sdf2gdf import sdf2gdf
from common import sm_tag, bam_getter

@bam_getter
def bam2gdf(
        target_reg: str,
        control_reg: str,
        genome_build: str,
        target_gene: str,
        control_gene: str,
        create_gdf_file: str,
        output_dir: str,
        bam_file: List[str],
        bam_dir: Optional[str] = None,
        bam_list: Optional[str] = None,
        **kwargs
    ) -> None:
    """Convert BAM files to a GDF file.

    This command calculates read depth from BAM files and then outputs a
    GDF (GATK-DepthOfCoverage Format) file, which is one of the input 
    files for the Stargazer program. Even though ``gatk DepthOfCoverage`` 
    could still be used to make GDF files, we recommend that you use this 
    command because the former is too heavy (i.e. requires too much memory) 
    for such a simple task (i.e. counting reads). The latter uses 
    ``samtools depth`` under the hood, which is way faster and requires 
    way less memory. Another nice about using ``bam2gdf`` instead of 
    ``samtools depth`` is that everything is already parametrized for 
    compatibility with Stargazer. 

    .. note::
        You do NOT need to install ``samtools`` to run this command.

    Args:
        genome_build (str):
            Genome build ('hg19' or 'grc38').
        target_gene (str):
            Name of target gene (e.g. 'cyp2d6').
        control_gene (str):
            Name or region of control gene (e.g. ‘vdr’, 
            ‘chr12:48232319-48301814’)
        output_file (str):
            Write output to this file.
        bam_file (list[str]):
            Input BAM files.
        bam_dir (str, optional):
            Use all BAM files in this directory as input.
        bam_list (str, optional):
            List of input BAM files, one file per line.
    """
    # Parse keyward arguments from the decorator.
    input_files = kwargs["input_files"]
    
    sdf = bam2sdf(genome_build, target_gene, target_reg, control_gene, control_reg, input_files)
    sm = [sm_tag(x) for x in input_files]
    result = sdf2gdf(None, sm, f=StringIO(sdf))
    output_path = f"{output_dir}/{create_gdf_file}"
    with open(output_path, "w") as f:
        f.write(result)
