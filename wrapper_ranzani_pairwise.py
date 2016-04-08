import sys

folder = "/home/cradens/script/"

sys.path.append(folder)

from get_exec_lines import *

# bash file that will run majiq is being written to (and will be run from)  this location:
fw = open('/data/THelper/ranzani_2015/run_MAJIQ_all_pairwise.bash','w')

fw.write('majiq build /data/DB/hg19/ensembl.hg19.gff3 -conf /data/THelper/ranzani_2015/settings.txt --output ./build --nthreads 16\n')


m1, v1 = get_EXECdpsi_line('/data/THelper/ranzani_2015/settings.txt','T4_N','T4_1','./build','./majiq/dpsi/T4_NvsTH1','./voila/dpsi/T4_NvsTH1',show_all=False)

fw.write(m1+v1)

m1, v1 = get_EXECdpsi_line('/data/THelper/ranzani_2015/settings.txt','T4_N','T4_2','./build','./majiq/dpsi/T4_NvsTH2','./voila/dpsi/T4_NvsTH2',show_all=False)

fw.write(m1+v1)

m1, v1 = get_EXECdpsi_line('/data/THelper/ranzani_2015/settings.txt','T4_N','T4_17','./build','./majiq/dpsi/T4_NvsTH17','./voila/dpsi/T4_NvsTH17',show_all=False)

fw.write(m1+v1)

m1, v1 = get_EXECdpsi_line('/data/THelper/ranzani_2015/settings.txt','T4_N','T4_R','./build','./majiq/dpsi/T4_NvsR','./voila/dpsi/T4_NvsR',show_all=False)

fw.write(m1+v1)

m1, v1 = get_EXECdpsi_line('/data/THelper/ranzani_2015/settings.txt','T4_N','T4_CM','./build','./majiq/dpsi/T4_NvsCM','./voila/dpsi/T4_NvsCM',show_all=False)

fw.write(m1+v1)

m1, v1 = get_EXECdpsi_line('/data/THelper/ranzani_2015/settings.txt','T4_N','T4_EM','./build','./majiq/dpsi/T4_NvsEM','./voila/dpsi/T4_NvsEM',show_all=False)

fw.write(m1+v1)

samps = ["T4_N","T4_1","T4_2","T4_17","T4_R","T4_CM","T4_EM"]

# Psi calculations
for s in samps:
    m1,v1 = get_EXECpsi_line('/data/THelper/ranzani_2015/settings.txt',s,'./build/','./majiq/psi/%s'%s,'./voila/psi/%s'%s)
    fw.write(m1+v1)
    

fw.close()
