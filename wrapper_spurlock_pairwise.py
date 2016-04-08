import sys

folder = "/home/cradens/script/"

sys.path.append(folder)

from get_exec_lines import *

# bash file that will run majiq is being written to (and will be run from)  this location:
fw = open('/data/THelper/run_MAJIQ_all_pairwise.bash','w')

fw.write('majiq build /data/DB/hg19/ensembl.hg19.gff3 -conf /data/THelper/settings.txt --output ./build --nthreads 16\n')

ll = '''Th1_p_A=SRR1817386
Th1_e_A=SRR1817389
Th2_p_A=SRR1817387
Th2_e_A=SRR1817390
Th17_p_A=SRR1817388
Th17_e_A=SRR1817391
Th1_p_r=SRR1817392
Th1_e_r=SRR1817395
Th2_p_r=SRR1817393
Th2_e_r=SRR1817396
Th17_p_r=SRR1817394
Th17_e_r=SRR1817397'''.split()

samps = [x.split('=')[0] for x in ll]

print "Effector: TH1<->TH2, TH1<->TH17, TH2<->TH17"

m1, v1 = get_EXECdpsi_line('./settings.txt','Th1_e_A','Th2_e_A','./build','./majiq/dpsi/TH1e_TH2e','./voila/dpsi/TH1e_TH2e',show_all=False)

fw.write(m1+v1)

m1, v1 = get_EXECdpsi_line('./settings.txt','Th1_e_A','Th17_e_A','./build','./majiq/dpsi/TH1e_TH17e','./voila/dpsi/TH1e_TH17e',show_all=False)

fw.write(m1+v1)

m1, v1 = get_EXECdpsi_line('./settings.txt','Th2_e_A','Th17_e_A','./build','./majiq/dpsi/TH2e_TH17e','./voila/dpsi/TH2e_TH17e',show_all=False)

fw.write(m1+v1)

print "Primary: TH1<->TH2, TH1<->TH17, TH2<->TH17"

m1, v1 = get_EXECdpsi_line('./settings.txt','Th1_p_A','Th2_p_A','./build','./majiq/dpsi/TH1p_TH2p','./voila/dpsi/TH1p_TH2p',show_all=False)

fw.write(m1+v1)

m1, v1 = get_EXECdpsi_line('./settings.txt','Th1_p_A','Th17_p_A','./build','./majiq/dpsi/TH1p_TH17p','./voila/dpsi/TH1p_TH17p',show_all=False)

fw.write(m1+v1)

m1, v1 = get_EXECdpsi_line('./settings.txt','Th2_p_A','Th17_p_A','./build','./majiq/dpsi/TH2p_TH17p','./voila/dpsi/TH2p_TH17p',show_all=False)

fw.write(m1+v1)

print "Primary <-> Effector: TH1, TH2, TH17"

m1, v1 = get_EXECdpsi_line('./settings.txt','Th1_p_A','Th1_e_A','./build','./majiq/dpsi/TH1p_TH1e','./voila/dpsi/TH1p_TH1e',show_all=False)

fw.write(m1+v1)

m1, v1 = get_EXECdpsi_line('./settings.txt','Th2_p_A','Th2_e_A','./build','./majiq/dpsi/TH2p_TH2e','./voila/dpsi/TH2p_TH2e',show_all=False)

fw.write(m1+v1)

m1, v1 = get_EXECdpsi_line('./settings.txt','Th17_p_A','Th17_e_A','./build','./majiq/dpsi/TH17p_TH17e','./voila/dpsi/TH17p_TH17e',show_all=False)

fw.write(m1+v1)

print "PolyA <-> Total: primary <-> primary"

m1, v1 = get_EXECdpsi_line('./settings.txt','Th1_p_A','Th1_p_r','./build','./majiq/dpsi/TH1p_PolyAvTotal','./voila/dpsi/TH1p_PolyAvTotal',show_all=False)

fw.write(m1+v1)

m1, v1 = get_EXECdpsi_line('./settings.txt','Th2_p_A','Th2_p_r','./build','./majiq/dpsi/TH2p_PolyAvTotal','./voila/dpsi/TH2p_PolyAvTotal',show_all=False)

fw.write(m1+v1)

m1, v1 = get_EXECdpsi_line('./settings.txt','Th17_p_A','Th17_p_r','./build','./majiq/dpsi/TH17p_PolyAvTotal','./voila/dpsi/TH17p_PolyAvTotal',show_all=False)

fw.write(m1+v1)

print "PolyA <-> Total: effector <-> effector"

m1, v1 = get_EXECdpsi_line('./settings.txt','Th1_e_A','Th1_e_r','./build','./majiq/dpsi/TH1e_PolyAvTotal','./voila/dpsi/TH1e_PolyAvTotal',show_all=False)

fw.write(m1+v1)

m1, v1 = get_EXECdpsi_line('./settings.txt','Th2_e_A','Th2_e_r','./build','./majiq/dpsi/TH2e_PolyAvTotal','./voila/dpsi/TH2e_PolyAvTotal',show_all=False)

fw.write(m1+v1)

m1, v1 = get_EXECdpsi_line('./settings.txt','Th17_e_A','Th17_e_r','./build','./majiq/dpsi/TH17e_PolyAvTotal','./voila/dpsi/TH17e_PolyAvTotal',show_all=False)

fw.write(m1+v1)

# Psi calculations
for s in samps:
    m1,v1 = get_EXECpsi_line('./settings.txt',s,'./build/','./majiq/psi/%s'%s,'./voila/psi/%s'%s)
    fw.write(m1+v1)
    

fw.close()
