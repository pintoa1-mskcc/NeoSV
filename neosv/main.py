import os
import sys
from .args import create_arg_parser
from .input import bedpe_load, ensembl_load, get_window_range
from .sv_utils import sv_pattern_infer_bedpe, remove_duplicate
from .annotation_utils import sv_to_sveffect
from .fusion_utils import sv_to_svfusion
from .sequence_utils import set_nt_seq, set_aa_seq, wt_and_mut_altered_aas,generate_neoepitopes #wt_and_mut_altered_NTs
from .netMHC import netmhc_pep_prep_fasta,netmhc_pep_prep
from .output import write_annot


def main():
    # generate the args
    args = create_arg_parser()
    # load a pyensembl Genome class
    ensembl = ensembl_load(args.release, args.gtffile, args.cdnafile, args.cachedir)

    if args.svfile.endswith('.bedpe'):
        # load vcf file and transform SVs into BEDPE class
        sv_beds = bedpe_load(args.svfile)
        # transform SV from BEDPE class to StructuralVariant class
        svs = [sv_pattern_infer_bedpe(sv_bed) for sv_bed in sv_beds]
    else:
        sys.exit('The input file must end with .bedpe. Other format is not allowed.')
    # remove duplicated SVs
    svs = remove_duplicate(svs)

    # annotate SVs and write to file_anno
    sv_effects = [sv_to_sveffect(sv, ensembl, args.complete) for sv in svs]
    sv_effects_flat = [sv_effect_unit for sv_effect in sv_effects for sv_effect_unit in sv_effect]
    file_anno = os.path.join(args.outdir, args.prefix + '.anno.txt')
    write_annot(file_anno, sv_effects_flat)

    if not args.anno:
        # specify the lengths of neoepitopes you want to predict
        pad_length = int(args.pad_length)
        # transform StructuralVariant class to SVFusion class
        sv_fusions = [sv_to_svfusion(sv, ensembl) for sv in svs]
        # remove those empty SVFusions
        sv_fusions = [sv_fusion for sv_fusion in sv_fusions if not sv_fusion.is_empty()]

        # predict neoepitopes for each SVFusion
        for sv_fusion in sv_fusions:
            sv_fusion.nt_sequence = set_nt_seq(sv_fusion)
            sv_fusion.aa_sequence_noutr = set_aa_seq(sv_fusion.nt_sequence_cds)
            
            sv_fusion.aa_sequence = set_aa_seq(sv_fusion.nt_sequence)
            
            sv_fusion = wt_and_mut_altered_aas(sv_fusion,pad_length)

        # predict binding affinity using netMHC
        file_netmhc_in_wt = os.path.join(args.outdir, args.prefix + '.WT.net.in.txt')
        file_netmhc_in = os.path.join(args.outdir, args.prefix + '.net.in.txt')
        if sv_fusions:
            netmhc_pep_prep_fasta(file_netmhc_in_wt,file_netmhc_in, sv_fusions)
        else:
            open(file_netmhc_in, 'w').close()
            open(file_netmhc_in_wt, 'w').close()

if __name__ == '__main__':
    main()
