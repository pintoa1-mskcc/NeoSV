import subprocess
import re


def netmhc_pep_prep(filepath, svfusions):
    """
    :param svfusions: a list of svfusion classes
    :param filepath: the output file
    :return: None
    """
    print("Preparing peptides for netMHC.")
    with open(filepath, 'w') as f:
        for svfusion in svfusions:
            if svfusion.neoepitopes:
                for neoepitope in svfusion.neoepitopes:
                    f.write(neoepitope + '\n')


def netmhc_pep_prep_fasta(filepath_WT,filepath_MUT, svfusions):
    """
    :param svfusions: a list of svfusion classes
    :param filepath: the output file
    :return: None
    """
    print("Preparing FASTAs for netMHC.")
    with open(filepath_MUT, 'w') as f_mut, open(filepath_WT,'w') as f_wt:
        for svfusion in svfusions:
            if svfusion.aa_sequence:
                id = svfusion.sv.id + svfusion.sv.svtype 
                f_mut.write('>' + id + '\n')
                f_mut.write(svfusion.mt_altered_aa + '\n')
                
                if svfusion.wt_altered_aa2 != None: 
                        
                    f_wt.write('>'+ id + 'W1\n')
                    f_wt.write(svfusion.wt_altered_aa1 + '\n')
                    f_wt.write('>' + id + 'W2\n')
                    f_wt.write(svfusion.wt_altered_aa2 + '\n')
                    
                else:
                    f_wt.write('>' + id + 'W\n')
                    f_wt.write(svfusion.wt_altered_aa1 + '\n')
                
def netmhc_run(netmhcpath, peppath, alleles, outpath):
    """
    :param netmhcpath: absolute path of the netmhc execution file
    :param peppath: input peptides file
    :param alleles: HLA alleles separated by ,
    :param outpath: outfile
    :return: None
    """
    print("Running netMHCpan.")
    cmd = netmhcpath \
        + ' -a ' + alleles \
        + ' -f ' + peppath \
        + ' -inptype 1' \
        + ' -BA' \
        + ' > ' + outpath
    p = subprocess.Popen(cmd, shell=True)
    p.wait()


def netmhc_reload(resultpath):
    """
    :param resultpath: the result of netmhc
    :return: a dictionary {neoepitope: {allele: [affinity, rank, FILTER]}}
    """
    pep_dic = {}
    with open(resultpath, 'r') as f:
        for line in f:
            if re.match(r'^\s+1\s+H', line):
                tmpline = re.split(r'\s+', line)
                allele = tmpline[2]
                pep = tmpline[3]
                affinity = float(tmpline[16])
                ba_rank = float(tmpline[15])
                el_rank = float(tmpline[13])
                if pep not in pep_dic:
                    pep_dic[pep] = {}
                pep_dic[pep][allele] = [affinity, ba_rank, el_rank, 'FILTER']
    return pep_dic


def netmhc_filter(pep_dic, aff_thre, ba_rank_thre, el_rank_thre):
    """
    :param pep_dic: the dictionary for neoepitopes generated by netmhc_load
    :param aff_thre: threshold for binding affinity
    :param ba_rank_thre: threshold for BA rank percentile
    :param el_rank_thre: threshold for EL rank percentile
    :return: a dictionary in which passed neoepitopes are labeled PASS
    """
    print("Filtering NetMHCpan results")
    for pep in pep_dic:
        for allele in pep_dic[pep]:
            if pep_dic[pep][allele][0] <= aff_thre and pep_dic[pep][allele][1] <= ba_rank_thre and pep_dic[pep][allele][2] <= el_rank_thre:
                pep_dic[pep][allele][3] = 'PASS'
    return pep_dic


