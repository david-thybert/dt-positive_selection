#!/Users/dthybert/bin//nextflow

params.ortho = "$projectDir/data/orthologues.txt"
params.nuc = "$projectDir/data/nuc/"
params.pep = "$projectDir/data/pep/"
params.out = "$projectDir/out/"
params.prank_command = "$projectDir/ext/prank-msa/prank"



log.info"""\


    =======================================================
        POSSEL : A POSITIVE SELECTiON ANALYSIS PIPELINE
    =======================================================

    Parameters: 
    ortho: $params.ortho
    nuc: $params.nuc
    pep: $params.pep
    out: $params.out
    prank location: $params.prank_command 

"""

process FormatFasta{

    publishDir params.out, mode: 'copy'

    input:
        path ortho
        val nucleotides
        val peptides
        val out_base

    output:
        path "*.nuc.fasta" , emit: ortho_nuc
    script:
    """
        python $projectDir/scripts/format_fasta.py --ortho $ortho --nuc $nucleotides --pep $peptides --out ./
    """
}


process AlignPeptide{

    publishDir params.out, mode: 'copy'

    input:
        path ortho_nuc

    output:
        path "*.nuc.fasta.best.fas", emit: align_pep

    script:
    """
        $params.prank_command  -d=$ortho_nuc -o=$ortho_nuc -codon -F
    """
}


process PositiveSelection{

    input:
        path m_align

    output:
        path
    script:
    """
        hyphy aBSREL --alignment
    """

}

// mafft --quiet --auto $ortho_pep > $ortho_pep-aln
workflow{
    
    ortho_dir = FormatFasta(params.ortho, params.nuc, params.pep, params.out)
    align_pep_ch = AlignPeptide(ortho_dir.ortho_nuc.flatten())
    PositiveSelection(align_pep_ch)
    align_pep_ch.align_pep.view()


}