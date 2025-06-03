import MDAnalysis as mda
import MDAnalysis.analysis.rdf
import numpy as np
import sys 
import scipy.integrate as sci
#import scipy.optimize as optimize
#import scipy.stats
import matplotlib.pyplot as plt
import math
#import solvation_analysis as solvation

if '-h' in sys.argv:
    print('This code calculates the coordination with respect to a central reference')
    print('python *.py tpr_file trj_file step sample "selection 1" "selection 2" cutoff fileout')
    exit()
                                                                                                  
tpr, trj, step, sample, sel1, sel2, cutoff, fileout = sys.argv[1:]
print(sys.argv)

def read_in_data(u,g1,g2,g3,cutoff,step,sample,i_rad,o_rad,sel1):
    t = []
    pos = []
    rdf_arr = []
    for i in range(0,len(u.trajectory),int(step)):
        ts = u.trajectory[i]
        pos.append(g1.positions[:,2]-g3.center_of_geometry()[2])
        t.append(ts.time)
    t=np.array(t)
    pos = np.array(pos)
    zeros = np.where(t == 0.0)[0]
    print(zeros)
    if len(zeros) > 1:
        for i in range(0,len(zeros)-1):
            if zeros[i+1]-zeros[i] > 1:
                rdf_arr.append((np.average(pos[zeros[i]:zeros[i+1]]),rdf(u,g1,20.0,1,100,g2,zeros[i],zeros[i+1]-1)))
                #print(pos[zeros[i]:zeros[i+1]])
    else:
        rdf_arr.append((np.average(pos),rdf(u,g1,20.0,1,100,g2,0,len(u.trajectory))))
    rdf_arr = np.array(rdf_arr)
    rdf_arr = rdf_arr[rdf_arr[:,0].argsort()]
    coord = []
    #print(rdf_arr[0][1][:,0]/10-i_rad[sel1.split()[1]])
    for i in range(len(rdf_arr)):
        ind=np.argmin(abs(rdf_arr[i][1][:,0]/10-i_rad[sel1.split()[1]]))
        indo=np.argmin(abs(rdf_arr[i][1][:,0]/10-o_rad[sel1.split()[1]]))
        coord.append((rdf_arr[i,0],rdf_arr[i][1][ind,2],rdf_arr[i][1][indo,2],rdf_arr[i][1][:,0],rdf_arr[i][1][:,1]))
        #print(ind,coord)
    return np.array(coord)


def rdf(u,g1,dist,step,bins,g2,start,stop): 
    rdfsum = np.zeros(bins)
    rdffinal = np.zeros((bins,3))
#    try:
#        g2 = u.select_atoms('{}'.format(gname),updating=True)
#    except:
#        g2 = u.select_atoms('{}'.format(gname),updating=True)
    #print(len(g2))
    if g1 == g2:
        rdf = mda.analysis.rdf.InterRDF(g1,g2,nbins=bins,range=[0.0,dist],exclusion_block=[1,1])
        rdf.run(start=start, stop=stop,step=step)
    else:    
        rdf = mda.analysis.rdf.InterRDF(g1,g2,nbins=bins,range=[0.0,dist])
        rdf.run(start=start, stop=stop,step=step)
    V=np.cumprod(u.dimensions)[2]        
    vol = 4*np.pi*rdf.results.bins**2
    rdf_vol = np.where(np.isinf(rdf.results.rdf),0,rdf.results.rdf) 
    rdf_vol = np.where(np.isnan(rdf_vol),0,rdf_vol)
    rho=np.cumsum(rdf.results.count)/np.cumprod(u.dimensions)[2]*vol
    coord=sci.cumulative_trapezoid(vol*rdf_vol,rdf.results.bins,initial=0)*len(g2)/V
    #print(g2,rdf_vol,coord,rdf.results.count)
    rdffinal[:,0],rdffinal[:,1],rdffinal[:,2] = rdf.results.bins,rdf_vol,coord
    return rdffinal

def write_data(dz,header,fileout):
    c = open('{}.dat'.format(fileout),'w')
    c.write(",".join(header)+'\n')
    for i in range(len(dz[:,0])):
        c.write("{},{},{}\n".format(dz[i,0],dz[i,1],dz[i,2]))
    c.close()
    c = open('{}rdf.dat'.format(fileout),'w')
    c.write(",".join(header)+'\n')
    c.write('x,{}\n'.format(','.join(str(dz[0,3])[1:-1].split())))
    print(','.join(str(dz[0,3])[1:-1].split()))
    print(','.join(str(dz[0,4])[1:-1].split()))
    for i in range(len(dz[:,0])):
        c.write("{},{}\n".format(dz[i,0],','.join(str(dz[i,4])[1:-1].split())))
    c.close()
    return
    

def main(tpr,trj,step,sample,sel1,sel2,cutoff,fileout):
#K	0.392206937	0.601185962
#Na	0.359636354	0.576397103
#Li	0.311436292	0.5412059
#Ca	0.360576931	0.578186312
#Mg	0.322654642	0.548322165
#La	0.375424976	0.585809656
#Nd	0.376844245	0.592589563
#Eu	0.364447957	0.575950879
#Tb	0.367208409	0.583484511
#Yb	0.354718215	0.568010183

    i_rad = {'K':3.92206937e-1,'NA':0.359636354,'LI':0.311436292,'CA':0.360576931,'MG':0.322654642,'La':0.375424976,'Nd':0.376844245,'Eu':0.364447957,'Tb':0.367208409,'Yb':0.354718215}
    o_rad = {'K':0.601185962,'NA':0.576397103,'LI':0.5412059,'CA':0.578186312,'MG':0.548322165,'La':0.585809656,'Nd':0.592589563,'Eu':0.575950879,'Tb':0.583484511,'Yb':0.568010183}
    if trj[-3:] in ['xtc','trr','tpr']: 
        print((trj.split('\n')))
        u=mda.Universe(tpr,trj.split('\n'),continuous=False)
        print(len(u.trajectory))
        g1 = u.select_atoms('%s' % sel1)
        g2 = u.select_atoms('%s' % sel2)
        g3 = u.select_atoms('resname UNL and (name O* and bonded name C*)')
        print(g1,g2,g3)
        coord_prof=read_in_data(u,g1,g2,g3,cutoff,step,sample,i_rad,o_rad,sel1)

    write_data(coord_prof,('#z','CN {}'.format(sel2),'CNout {}'.format(sel2)),fileout)

if __name__=="__main__":
    main(tpr,trj,step,sample,sel1,sel2,cutoff,fileout)




