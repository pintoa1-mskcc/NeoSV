from Bio.Seq import Seq


def set_aa_seq(nt_sequence):
    """
    :function: translate the nucelotide sequence in svfusion to amino acid sequence
    :param: svfusion: a SVFusion class
    :return: the corresponding amino acid sequence
    :NOTE: if there is a stop codon, the process will stop before translating all nucleotides
    """
    dna_seq = Seq(trim_to_3x(nt_sequence))
    mrna_seq = dna_seq.transcribe()
    aa_seq = mrna_seq.translate(to_stop=True)
    return str(aa_seq)

def set_aa_seq_no_utr(nt_sequence):
    """
    :function: translate the nucelotide sequence in svfusion to amino acid sequence
    :param: svfusion: a SVFusion class
    :return: the corresponding amino acid sequence
    :NOTE: if there is a stop codon, the process will stop before translating all nucleotides
    """
    dna_seq = Seq(trim_to_3x(nt_sequence))
    mrna_seq = dna_seq.transcribe()
    aa_seq = mrna_seq.translate(to_stop=True)
    return str(aa_seq)


def set_nt_seq(svfusion):
    """
    :function: given a svfusion, paste the nucleotides in cc_1, cc_2 and 3utr
    :param svfusion: a SVFusion class
    :return: the nucleotide sequence from start codon to the end of 3utr
    """
    return svfusion.nt_sequence_cds + svfusion.nt_sequence_3utr


def trim_to_3x(nt_sequence):
    """
    :function: trim a sequence to be divisible by 3
    :return: trimmed sequence
    """
    remainder = len(nt_sequence) % 3
    if remainder:
        return nt_sequence[: -remainder]
    else:
        return nt_sequence
    
def wt_and_mut_altered_aas(svfusion,pad_len):
    mt = svfusion.aa_sequence_noutr
    
    if svfusion.cc_1.transcript.protein_sequence == svfusion.aa_sequence:
        svfusion.aa_sequence = None
        return svfusion
    if svfusion.cc_1.part == "5" or svfusion.cc_1.transcript.protein_sequence[0:10] == mt[0:10]:
            
        wt1 = svfusion.cc_1.transcript.protein_sequence
        wt2 = svfusion.cc_2.transcript.protein_sequence
        
    else:
        wt2 = svfusion.cc_1.transcript.protein_sequence
        wt1 = svfusion.cc_2.transcript.protein_sequence
        
    len_from_start1 = len_from_end2 = 0

    ## from start
    for i in range(0, len(wt1)):
        len_from_start1 = i
        if wt1[i:i + 1] != mt[i:i + 1]:
            break

    ## from end
    wt2_rev = wt2[::-1]
    mt_rev = mt[::-1]
    

    
    for i in range(0, len(wt2)):
        len_from_end2 = i
        
        if wt2_rev[i:i + 1] != mt_rev[i:i + 1]:
            break
        
    if svfusion.cc_1.part == "5"  or svfusion.cc_1.transcript.protein_sequence[0:10] == mt[0:10] :
            
        svfusion.wt_altered_aa1 = wt1[max(0, len_from_start1 - pad_len):min(len(wt1), len_from_start1 )]
        
        if svfusion.frame_effect != "In-Frame":
            svfusion.wt_altered_aa2 = wt2[max(0, len_from_end2 - pad_len):min(len(wt2), len_from_end2 + pad_len-1)]
        else:
            svfusion.wt_altered_aa2 = wt2[min(0, len_from_end2 - pad_len):min(len(wt2), 0 + pad_len)]
        
        svfusion.mt_altered_aa = filter_mutatedseq(svfusion,pad_len)
    else:
        if svfusion.frame_effect != "In-frame":
            
            svfusion.wt_altered_aa1 = wt2[max(0, len_from_end2 - pad_len):max(len(wt2), len_from_end2 + pad_len-1)]
            
        else:
            svfusion.wt_altered_aa1 = wt2[min(0, len_from_end2 - pad_len):min(len(wt2), min(0, len_from_end2 - pad_len + 1) + pad_len)]
        
        svfusion.wt_altered_aa2 = wt1[max(0, len_from_start1 - pad_len):max(len(wt1), len_from_start1 + pad_len)]
        svfusion.mt_altered_aa = filter_mutatedseq(svfusion,pad_len)
        
    if svfusion.wt_altered_aa1 == svfusion.wt_altered_aa2: 
        svfusion.wt_altered_aa2 = None
    return(svfusion)

def filter_mutatedseq(svfusion,pad_len):
    mt = svfusion.aa_sequence
    
    if svfusion.cc_1.part == "5" or svfusion.cc_1.transcript.protein_sequence[0:4] == mt[0:4]:
            
        wt = svfusion.cc_1.transcript.protein_sequence
        
        for i in range(0, len(wt)):
            len_from_start1 = i
            if wt[i:i + 1] != mt[i:i + 1]:
                break
        
        wt_rev = wt[::-1]
        
        mt_rev = mt[::-1]
      
        for i in range(0, len(wt_rev)):
            len_from_end1 = i
            if wt_rev[i:i + 1] != mt_rev[i:i + 1]:
                break
            
        if len_from_end1>1:
            return mt[max(0, len_from_start1 - pad_len ):]
        else:
                
            mt= mt[max(0, len_from_start1 - pad_len):]
            return mt[:len(mt)-len_from_end1+1]
            
    
    elif svfusion.cc_2.transcript.protein_sequence[0:4] == mt[0:4]:
        wt = svfusion.cc_2.transcript.protein_sequence
        
        for i in range(0, len(wt)):
            len_from_start1 = i
            if wt[i:i + 1] != mt[i:i + 1]:
                break
        
        wt_rev = wt[::-1]
        mt_rev = mt[::-1]
      
        for i in range(0, len(wt_rev)):
            len_from_end1 = i
            if wt_rev[i:i + 1] != mt_rev[i:i + 1]:
                break
        if len_from_end1>0:
            return mt[max(0, len_from_start1 - pad_len):]
        else:
                
            mt= mt[max(0, len_from_start1 - pad_len):]
            
            return mt[:len(mt)-len_from_end1+1]
            
    
    else: 
        return mt
        
def generate_neoepitopes(svfusion, window_range):
    """
    :function: cut the WT and MUT protein sequence by window_range, then get the MUT specific peptides
    :param svfusion: a SVFusion class
    :param window_range: a list, specifying the range of window size, e.g. [8,9,10,11]
    :return: a list of neopeptides for svfusion
    """
    mut_peptides = []
    if svfusion.aa_sequence != None:
            
        for window in window_range:
            mut_peptides = mut_peptides + cut_sequence(svfusion.mt_altered_aa, window)
        wt_peptides = []
        for window in window_range:
            wt_peptides = wt_peptides + cut_sequence(svfusion.cc_1.transcript.protein_sequence, window)
            wt_peptides = wt_peptides + cut_sequence(svfusion.cc_2.transcript.protein_sequence, window)
        return list(set(mut_peptides)-set(wt_peptides))


def cut_sequence(aa_sequence, window):
    """
    :param aa_sequence: amino acid sequence
    :param window: the size of window, window should be smaller than len(aa_sequence)
    :return: all possible slices with length = window in this sequence
    """
    if len(aa_sequence) < window:
        return []
    else:
        return [aa_sequence[i: i+window] for i in range(len(aa_sequence)-window+1)]
