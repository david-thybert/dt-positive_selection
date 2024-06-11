#!/Users/dthybert/bin//nextflow

params.ortho = "$projectDir/data/orthologues.txt"
params.nuc = "$projectDir/data/nuc/"
params.pep = "$projectDir/data/pep/"
params.out = "$projectDir/out/"
params.prank_command = "$projectDir/ext/prank-msa/prank"
params.mafft_comand = "$projectDir/ext/mafft-mac/mafft.bat"
params.zorro_command = "$projectDir/ext/zorro-master/bin/zorro"
params.fasttree_command = "$projectDir/ext/zorro-master/bin/FastTree" 
params.zorro_thr = "5.0"
params.pal2nal = "$projectDir/ext/pal2nal.v14/pal2nal.pl"
params.codeml_command = "$projectDir/ext/paml-master/bin/codeml"
params.raxml_command = "$projectDir/ext/raxml/raxml-ng"
params.fgrd_species = "$projectDir/data/foreground.txt"
params.ancestral = "False"
params.possel_method = "ABSREL" // ABSREL | PAML_BRST
if (params.possel_method == "PAML_BRST" )
{params.tag_fgrd = "#1"}
else
{params.tag_fgrd = "{fgrd}"}



log.info"""\

    ==================================================================
    |                                                                |
    |    *****      *****      ******    ******   ******   *         |
    |    *    *   *       *   *         *         *        *         |
    |    *    *  *         *  *         *         *        *         |
    |    *****   *         *   ******    ******   ******   *         |
    |    *       *         *         *         *  *        *         |
    |    *        *       *          *         *  *        *         |
    |    *          *****      ******     *****   ******   ******    |
    |                                                                |
    |                                                                |
    |             A POSitive SELection analysis pipeline             |
    |                                                                |
    ==================================================================
                   
"""
/*
    Parameters: 
    ortho: $params.ortho
    nuc: $params.nuc
    pep: $params.pep
    out: $params.out
    prank location: $params.prank_command 
    raxml location: $params.raxml_command
*/

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
   path ortho_nuc: nucleotide fasta file formated per orthogroup
   path ortho_pep: peptide fasta file formated per orthogroup
   */
    publishDir params.out, mode: 'copy'

    input:
        path ortho
        val nucleotides
        val peptides

    output:
        path "*.nuc.fasta" , emit: ortho_nuc
        path "*.pep.fasta" , emit: ortho_pep
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
        path ortho_pep

    output:
        path "*.ali", emit: aligned_pep

    script:
    """
     $params.mafft_comand --amino --localpair $ortho_pep > ${ortho_pep}.ali
    """
    //$params.prank_command  -d=$ortho_nuc -o=$ortho_nuc -codon -F
}


process FilterNonConfidentColumns{
/*This process filter non confident column using zorro

  input:
  path align_seq: path to the alignment file
  output:
  path align_filt: path to the filterted alignment
  path reg_filt: path to th efile describing the filtered regions
*/
    publishDir params.out, mode: 'copy'

    input:
        path align_seq

    output:
        path "*.filt.fa", emit: align_filt
        path "*.conf_reg", emit : reg_filt
    
    script:
    """
        python $projectDir/scripts/zorro_wrapper.py --mult $align_seq --zorro $params.zorro_command --tree_cmd $params.fasttree_command --out ${align_seq}.filt.fa --threshold $params.zorro_thr
    """

}

process PepAli_2_DNAAli{
/*This process convert a peptide multiple sequence alignemnt into a DNA multiple sequence alignemnt
  input:
  val pair_pepali_nuc : mapping betwen id, path to pep align, path to nuc file [id, pep_al, nuc_file]
  output:
  path align_nuc : path to the aligned nucleotide sequences
*/
    publishDir params.out, mode: 'copy'

    input:
        val pair_pepali_nuc

    output:
        path "*.nuc.ali.fasta", emit: align_nuc
    
    script:
    """
        perl $params.pal2nal -output fasta ${pair_pepali_nuc[1]} ${pair_pepali_nuc[2]} > ./${pair_pepali_nuc[0]}.nuc.ali.fasta
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
        python $projectDir/scripts/tag_tree_nodes.py --tree $tree --spe_tag $params.fgrd_species --out ${tree}.tagged --tag "$params.tag_fgrd"  --ancestral $params.ancestral
    """

}
process TestSaturation{
/*this process test whether a gene multiple alignemnt reached a substitution saturation
   input:
   path align_nuc_filt: path of the multiple DNA sequence aligment 
   output:
   path saturation_info: path the information reltive to saturation
*/
    publishDir params.out, mode: 'copy'
    input:
        path align_nuc_filt
    output:
        path "*.sat", emit: saturation_info
    script:
    """
        python $projectDir/scripts/test_saturation_substitution.py --mult $align_nuc_filt --out ${align_nuc_filt}.sat
    """

}

process PositiveSelectionABSREL{
/* this process rune the positive selction analysis using 
   the aBSREL model.

   input:
   val pair_nuc_tree_ch: mapping betwen id, path to dna align, path tree [id, DNA_al, tree]
   output:
   path pos_sel: path to the json file storing positive selection results.
*/

    publishDir params.out, mode: 'copy'

    input:
         val pair_nuc_tree_ch

    output:
        path  "*.ABSREL.json", emit: pos_sel

    script://{fgrd}
    """
       hyphy aBSREL --alignment ${pair_nuc_tree_ch[1]} --tree ${pair_nuc_tree_ch[2]} --branches fgrd --output ./${pair_nuc_tree_ch[0]}.ABSREL.json
    """
    
}

process CombinePosselInfoABSREL{
/* this process combine all the json file produced by
   Hyphy aBSREL model.

   input:
   path pos_sel_json: the list of hyphy json file from ABSREL model  to be combined 
   path sat_subst : the list of file with substitution saturation annotation
   output:
   path pos_sel: path to the csv combining and multitest correcting the output of possel
*/
    publishDir params.out, mode: 'copy'

    input:
        path pos_sel_json
        path sat_subst
    output:
        path "*.possel", emit: pos_sel

    script:
    """
         python $projectDir/scripts/combine_pos_sell_info_absrel.py --files_possel "$pos_sel_json" --files_sat_subst "$sat_subst" --prefix_out "aBSREL"
    """
}

process CreateCTLFile_null{
/*
 */
    publishDir params.out, mode: 'copy'

    input:
        val params_ctl
    
    output:
        path "*.${type_ctl}.ctl", emit: ctl_files
    
    script:
    """
        python $projectDir/scripts/create_ctl_paml.py --alignment ${params_ctl[1]}  --tree ${params_ctl[2]} --out_PAML ${params_ctl[0]}.null.ctd --out_file ${params_ctl[0]}.null.ctl --type null
    """
}

process CreateCTLFile_alt{
/*
 */
    publishDir params.out, mode: 'copy'

    input:
        val params_ctl
    
    output:
        path "*.alt.ctl", emit: ctl_files
    
    script:
    """
        python $projectDir/scripts/create_ctl_paml.py --alignment ${params_ctl[1]}  --tree ${params_ctl[2]} --out_PAML ${params_ctl[0]}.alt.ctd --out_file ${params_ctl[0]}.alt.ctl --type alt
    """
}
process RunCodeML_null{
/*
*/
    publishDir params.out, mode: 'copy'
    
    input:
        path ctl_file
    
    output:
        //path "*.${type_ctl}.ctd", emit: ctd_files
        path "*.null.ctd", emit: ctd_files
    script:
    """
        python $projectDir/scripts/run_codeml.py --codeml $params.codeml_command --ctl $ctl_file
    """
}

process RunCodeML_alt{
/*
*/
    publishDir params.out, mode: 'copy'
    
    input:
        path ctl_file
    
    output:
       // path "*.${type_ctl}.ctd", emit: ctd_files
       path "*.alt.ctd", emit: ctd_files
    
    script:
    """
        python $projectDir/scripts/run_codeml.py --codeml $params.codeml_command --ctl $ctl_file
    """
}

process Fasta2Phylip{
/*
*/
    publishDir params.out, mode: 'copy'

    input:
        path fasta_mult
    output:
        path "*.phy", emit: phy_mult
    script:
    """
        python $projectDir/scripts/fasta2phy.py --mult_fasta $fasta_mult --out_phy ${fasta_mult}.phy
    """
}

process CalculateCodemlPval{
/*
*/
    publishDir params.out, mode: 'copy'

    input:
        val pair_alt_null_ctd
    output:
        path "*.pml.brst", emit: pml_brst
    script:
    """
        python $projectDir/scripts/calculate_codeml_pval.py --alt_codeml ${pair_alt_null_ctd[1]}  --null_codeml ${pair_alt_null_ctd[2]} --out ${pair_alt_null_ctd[0]}.pml.brst
    """   
}

process CombinePosselInfoPML_BRST{
/* this process combine all the tsv file produced when implementing the branch site model 
with codeml

   input:
   path pos_sel_tsv: the list of hyphy tsv file from  codeml branch site 
   path sat_subst : the list of file with substitution saturation annotation
   output:
   path pos_sel: path to the csv combining and multitest correcting the output of possel
*/
    publishDir params.out, mode: 'copy'

    input:
        path pos_sel_tsv
        path sat_subst
    output:
        path "*.possel", emit: pos_sel

    script:
    """
         python $projectDir/scripts/combine_pos_sell_info_brst.py --files_possel "$pos_sel_tsv" --files_sat_subst "$sat_subst" --prefix_out "PML_BRST"
    """
}



//combine_pos_sell_info_brst
nextflow.enable.dsl=2

workflow{

    ortho_dir = FormatFasta(params.ortho, params.nuc, params.pep)
    
    // Align protein sequences
    align_pep_ch = AlignSequence(ortho_dir.ortho_pep.flatten())
    
    /////////// Filtering non informative steps //////////////
    // remove non confident alignment 
    align_pep_filt = FilterNonConfidentColumns(align_pep_ch)

    // combine two chanel to be used later
    ortho_dir_nuc_id = ortho_dir.ortho_nuc.flatten().map { [it.toString().split("/")[-1].split(".nuc")[0], it]}
    align_pep_filt_id = align_pep_filt.align_filt.map { [it.toString().split("/")[-1].split(".pep")[0], it]}
    pair_pepal_nuc_ch = align_pep_filt_id.combine(ortho_dir_nuc_id, by: 0)

    //convert prot alignemnt in DNA alignemnt
    nuc_ali_filt = PepAli_2_DNAAli(pair_pepal_nuc_ch)
    
    // build tree based on four fold degenrate sites
    four_four_deg_ch = ExtractFourFoldDegeSites(nuc_ali_filt.flatten())
    tree_ch = RaxmlPhylogeny(four_four_deg_ch.flatten())

    // test for saturation
    sat_info_ch = TestSaturation(four_four_deg_ch.flatten())

    if ( params.possel_method == "ABSREL" ) {// ABSREL model
        
        tag_fgrd = "{fgrd}"
        // tag trees with foreground tag
        //foreground_ch = channel.fromPath(params.fgrd_species)
        tree_tagged_ch = TagForgroundInTree(tree_ch.flatten())
        // combining 
        tree_tagged_id_ch = tree_tagged_ch.flatten().map { [it.toString().split("/")[-1].split(".nuc")[0], it]}
        nuc_ali_filt_id_ch = nuc_ali_filt.flatten().map { [it.toString().split("/")[-1].split(".nuc")[0], it]}
        pair_nuc_tree_ch = nuc_ali_filt_id_ch.combine(tree_tagged_id_ch, by: 0)
        //pair_nuc_tree_ch.view()

        // positive selection
        pos_sel_res_ch  = PositiveSelectionABSREL(pair_nuc_tree_ch)
    
        // combine result and multiple test correciton
        pos_sel_comb_ch = CombinePosselInfoABSREL(pos_sel_res_ch.collect(), sat_info_ch.collect())
    }
    else if(params.possel_method == "PAML_BRST" ){

        tag_fgrd = "#1"
        //foreground_ch = channel.fromPath(params.fgrd_species)
        tree_tagged_ch = TagForgroundInTree(tree_ch.flatten())
        
        // fasta2phylyp
        nuc_ali_filt_phy_ch = Fasta2Phylip(nuc_ali_filt.flatten())

        //combine alignment and tree
        tree_tagged_id_ch = tree_tagged_ch.flatten().map { [it.toString().split("/")[-1].split(".nuc")[0], it]}
        nuc_ali_filt_id_ch = nuc_ali_filt_phy_ch.flatten().map { [it.toString().split("/")[-1].split(".nuc")[0], it]}
        pair_nuc_tree_ch = nuc_ali_filt_id_ch.combine(tree_tagged_id_ch, by: 0)
       
        //create CTL null hyp and run codeml
        null_ctl_ch = CreateCTLFile_null(pair_nuc_tree_ch)
        null_ctd_ch = RunCodeML_null(null_ctl_ch)
        
        //create CTL alt hyp and run codeml
        alt_ctl_ch = CreateCTLFile_alt(pair_nuc_tree_ch)
        alt_ctd_ch = RunCodeML_alt(alt_ctl_ch)
        
        // map two ctd file to calulate based ont he gene id.
        null_ctd_id_ch = null_ctd_ch.flatten().map { [it.toString().split("/")[-1].split(".null")[0], it]}
        alt_ctd_id_ch = alt_ctd_ch.flatten().map { [it.toString().split("/")[-1].split(".alt")[0], it]}
        pair_null_alt_ctd_ch = alt_ctd_id_ch.combine(null_ctd_id_ch, by: 0)

        // calculate pval and store file
        pml_brst_ch = CalculateCodemlPval(pair_null_alt_ctd_ch)

        // combine all outfile into a possel file
        pos_sel_comb_ch = CombinePosselInfoPML_BRST(pml_brst_ch.collect(), sat_info_ch.collect())
    }


}