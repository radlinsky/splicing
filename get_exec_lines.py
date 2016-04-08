import itertools

def get_SRR_list(settings_file,condition):
    
    fd = open(settings_file,'r')
    find = 0
    for line in fd:
        if line.startswith(condition+'='): 
            hit = line.strip()
            find += 1
        if line.startswith('[readlen]'): break
    if find > 1: raise RuntimeError('WARNING: condition (%s) was not unique in the settings file!!!'%condition)
    if find == 0: raise RuntimeError('WARNING: condition (%s) was not found in the settings file!!!'%condition)
    SRR_list = hit.split('=')[1].split(',')
    return SRR_list

def get_majiqs_and_SGs(settings_file,cond1,build_dir):
    
    cond1_list = get_SRR_list(settings_file,cond1)
    cond1_majiqs=''
    cond1_SGs=''
    for m in cond1_list:
        cond1_majiqs += ' %s/%s.majiq'%(build_dir,m)
        cond1_SGs += ' %s/%s.splicegraph'%(build_dir,m)
    return cond1_majiqs, cond1_SGs
    
def get_EXECdpsi_line(settings_file,cond1,cond2,build_dir,majiq_outdir,voila_outdir,show_all=True,nthreads=16,dpsi_thresh=0.2,default_prior=False):
    
    cond1_majiqs, cond1_SGs = get_majiqs_and_SGs(settings_file,cond1,build_dir)
    cond2_majiqs, cond2_SGs = get_majiqs_and_SGs(settings_file,cond2,build_dir)    
    
    if show_all == True:
        write_thresh = '--show-all'
    else:
        write_thresh = '--threshold %s'%dpsi_thresh
    
    majiq_line='majiq deltapsi -grp1%s -grp2%s --nthreads %s --output %s --name %s %s'%(cond1_majiqs, cond2_majiqs, nthreads, majiq_outdir, cond1, cond2)
    if default_prior == True: majiq_line += ' --default_prior\n'
    else: majiq_line += '\n'
    voila_line='voila deltapsi %s/%s_%s.deltapsi_quantify.pickle -splice-graphs1%s -splice-graphs2%s %s -o %s\n'%(majiq_outdir,cond1,cond2,cond1_SGs,cond2_SGs,write_thresh,voila_outdir)
    return majiq_line,voila_line    

def get_EXECpsi_line(settings_file,cond1,build_dir,majiq_outdir,voila_outdir,nthreads=16):
    
    cond1_majiqs, cond1_SGs = get_majiqs_and_SGs(settings_file,cond1,build_dir)
    majiq_line='majiq psi%s --nthreads %s --output %s --name %s\n'%(cond1_majiqs, nthreads, majiq_outdir, cond1)
    voila_line='voila psi %s/%s_psigroup.pickle -splice-graphs1%s -o %s\n'%(majiq_outdir,cond1,cond1_SGs,voila_outdir)
    return majiq_line, voila_line

def all_dpsi_combs_bash(bash_outfile,cond_list,settings_file,build_dir,majiq_outdir,voila_outdir,nthreads=16):
    
    fw = open(bash_outfile,'w')
    pairs = list(itertools.combinations(cond_list,2))
    print('There are %s combinations of the list of %s conditions'%(len(pairs),len(cond_list)))
    for n in range(len(pairs)):
        c1 = pairs[n][0]
        c2 = pairs[n][1]
        m1,v1 = get_EXECdpsi_line(settings_file,c1,c2,build_dir,majiq_outdir+'/%s_%s/'%(c1,c2),voila_outdir+'/%s_%s/'%(c1,c2),nthreads=nthreads)    
        fw.write(m1+v1)
    fw.close()