#!/Users/dthybert/bin//nextflow

params.ortho = "$projectDir/data/orthologues.txt"
params.nuc = "$projectDir/data/nuc/"
params.pep = "$projectDir/data/pep/"
params.out = "$projectDir/out/"
params.prank_command = "$projectDir/ext/prank-msa/prank"
params.raxml_command = "$projectDir/ext/raxml/raxml-ng"
params.fgrd_species = "$projectDir/data/foreground.txt"
params.ancestral = "False"


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
    raxml location: $params.raxml_command

"""

process FormatFasta{
/* This process format fasta files to be compatible with positive selction analysis
   I look at each 1:1 orthougroups in orthology matrix provided in input and for each 
   orthogroup it selects the correpsonding sequecing in the genome or considered species
   the process check a consistancy between CDS nucleotid sequence and protein nucleotid 
   sequence

   input:
   path ortho : path to the orthology matrix
   val nucleotides: directory where to find the nucletide based CDS
   val peptides: directory where to find the peptide based sequeces

   output:
   path ortho_nuc: fasta file formated for positive selection
   */
    publishDir params.out, mode: 'copy'

    input:
        path ortho
        val nucleotides
        val peptides

    output:
        path "*.nuc.fasta" , emit: ortho_nuc
    script:
    """
        python $projectDir/scripts/format_fasta.py --ortho $ortho --nuc $nucleotides --pep $peptides --out ./
    """
}


process AlignSequence{
/*This process run the multiple sequence alignment of the orthogroup 
  sequence. The alignment is made using a codon aware parameter from pranks
  to improve the aligment of the protein coding sequence.
  
  input:
  path ortho_nuc: Path to the fasta file storing all the sequence form a similar orthogroup

  output:
  path align_seq: Path to multiple sequence alignment of the orthogroup.
*/

    publishDir params.out, mode: 'copy'

    input:
        path ortho_nuc

    output:
        path "*.nuc.fasta.best.fas", emit: align_seq

    script:
    """
        $params.prank_command  -d=$ortho_nuc -o=$ortho_nuc -codon -F
    """
}

process ExtractFourFoldDegeSites{
/*This process extract the four fold degenrated site form consrevd amino acid across 
  the species from the multiple alignemnt and produce an alignment file based exclusively 
  on the third codon site. This file can be used for modeling neutral evolution 

  input:
  path align_seq : path to the orthogroup multiple sequence alignment at the fasta format
  output:
  path four_four_deg: path to the multiple sequence alignment of four fold degenrated sites.
*/
    publishDir params.out, mode: 'copy'

    input:
        path align_seq

    output:
        path "*.ff", emit: four_four_deg
    
    script:
    """
        python $projectDir/scripts/four_fold_degen_site.py --nucal $align_seq --out ${align_seq}.ff 
    """
}

process RaxmlPhylogeny{
/*Ths process is running the phylogenitc tree reconstruction using Raxml
  and the four fold degenrated site alignment
  input:
  path four_fold_sites: path to the four fold degenerate site multiple alignment
                        file in fasta format
  output:
  path tree: the path tot he tree file in nwk format.
*/
    publishDir params.out, mode: 'copy'

    input:
        path four_fold_sites
    
    output:
        path  "*ff.raxml.bestTree", emit: tree
    
    script:
    """
        $params.raxml_command --msa $four_fold_sites --model GTR+G
    """

}

process Check_align_and_trees{
/* This process look for all alignemnt file with a correpsonding tree.
   and return all pair alignemnt tree in file 
   input:
   path align : list of alignments 
   path trees : list of trees
   output:
   path pair_align_tree: file storing all pair alignment/tree with a pair per line  
*/
    publishDir params.out, mode: 'copy'

    input:
        path align
        path trees
    
    output:
         path  "pair_align_tree.txt", emit: pair_align_tree
   
    script:
        """
            python $projectDir/scripts/check_align_and_trees.py  --ff_deg "$align" --trees "$trees" --base_dir ${params.out} --out pair_align_tree.txt
        """
}


process TagForgroundInTree{
/*This process tag all species belonging to the forground set
  in each phylogenetic trees.

  input:
  path tree: the phylogenetic tree
  path foreground: the file describing the forground species(one species per line)
  output:
  path tagged_tree: tree tagged with the forground
*/
    publishDir params.out, mode: 'copy'

    input:
        path tree
        //path foreground
    output:
        path "*.bestTree.tagged", emit: tagged_tree

    script:
    """
        python $projectDir/scripts/tag_tree_nodes.py --tree $tree --spe_tag $params.fgrd_species --out ${tree}.tagged --tag {fgrd} --ancestral $params.ancestral
    """

}
process PositiveSelectionABSREL{
/* this process rune the positive selction analysis using 
   the aBSREL model.

   input:
   tuple val(align), val(tree): align is the path to the alignment file and tree the path tot he correspding tree file
   output:
   path pos_sel: path to the json file storing positive selection results.
*/

    publishDir params.out, mode: 'copy'

    input:
        tuple val(align), val(tree)

    output:
        path  "*.ABSREL.json", emit: pos_sel

    script://{fgrd}
    """
       hyphy aBSREL --alignment $params.out/$align  --tree $params.out/$tree --branches fgrd --output ${align}.ABSREL.json
    """
    
}

nextflow.enable.dsl=2

workflow{

    ortho_dir = FormatFasta(params.ortho, params.nuc, params.pep)
    align_nuc_ch = AlignSequence(ortho_dir.ortho_nuc.flatten())
    
    // build tree based on four fold degenrate sites
    four_four_deg_ch = ExtractFourFoldDegeSites(align_nuc_ch.flatten())
    tree_ch = RaxmlPhylogeny(four_four_deg_ch.flatten())
    // tag trees with foreground tag
    //foreground_ch = channel.fromPath(params.fgrd_species)
    tree_tagged_ch = TagForgroundInTree(tree_ch.flatten())
    
    // wait that all  alignment and trees are finished 
    // and keep only alignemnt with a correpssondant tree
    lst_align_nuc_ch = align_nuc_ch.collect(sort: true)
    lst_tree_ch = tree_tagged_ch.collect(sort:true)
    pair_align_tree_file_ch = Check_align_and_trees(lst_align_nuc_ch, lst_tree_ch)
    
    //
    pair_align_tree_ch = pair_align_tree_file_ch.splitCsv(sep:"\t")
    pos_sel_res_ch  = PositiveSelectionABSREL(pair_align_tree_ch)
}