import MDAnalysis as mda
import numpy as np
import sys 
import scipy.integrate as integrate
import scipy.optimize as optimize
import scipy.stats
import matplotlib.pyplot as plt
import math

if '-h' in sys.argv:
    print('This code calculates the autocorrelation function')
    print('python *.py tpr_file trj_file step sample "selection 1" "selection 2" cutoff fileout')
    exit()

tpr, trj, step, sample, sel1, sel2, cutoff, fileout = sys.argv[1:]
sample=float(sample)

def stretched_e(t,tao,a):
    return np.exp(-(t/tao)**a)
def read_in_data(u,g1,g2,g3,cutoff,step,sample):
    t = []
    pos = []
    pos_av = []
    for i in range(0,len(u.trajectory),int(step)):
        ts = u.trajectory[i]
        pos.append(g1.positions-g3.center_of_geometry())
        t.append(ts.time)
    t=np.array(t)
    pos = np.array(pos)
    zeros = np.where(t == 0.0)[0]
    return t, pos, zeros


def read_in_files(filename):
    c = open(filename,'r')
    t=[]
    pos = []
    for line in c:
        #print(line)
        try:
            d = np.genfromtxt(line[:-1],skip_header=1)
        except:
            print(line[:-1])
        try:
            if abs(np.average(d[0,2]) - np.average(d[0,3])) < 0.5:  
                t.append(d[:,0])
                pos.append(d[:,2])
            else:
                print('Harmonic restraint outside of bounds: %s' % line)
        except:
            print(line[:-1],' Did not meet minimum requirements as a file')
        #d = open(line[:-1],'r') 
        #for dat in d:
#       #     print(dat)
        #    t.append(float(dat.split()[0]))
        #    pos.append(float(dat.split()[2]))
        #d.close()
    c.close()
    t=np.concatenate(t)
    pos = np.concatenate(pos)
    zeros = np.where(t == 0.0)[0]
    print(np.shape(pos))
    return t, pos, zeros




def self_var_tau(t,pos,cutoff,step,sample,ax,zeros):
#    sample = int(zeros[1]-zeros[0])/2.0
#    pos_block = [[(np.average(pos[zeros[j]:zeros[j]+sample]),np.average(pos[zeros[j]+sample:zeros[j]+2*sample]))]  for j in range(len(zeros)-1)]
#    print(pos_block)
    def stretched_e(x,a,b):
        return np.exp(-(x/a)**b)
    sample_ref=sample
    print(zeros)
    if len(zeros) == 1:
        print('Zeros len is 1')
        zeros = [0,len(pos)]
    dz = []
    z = []
    ct_collect = []
    cutoff=0.01
    for j in range(len(zeros)-1):
        if int((zeros[j+1]-zeros[j])/int(sample_ref/(t[zeros[j]+1]-t[zeros[j]]))) > 1:
#            n = 50
#            sample = int((zeros[j+1]-zeros[j])/n)
            traj_l = (zeros[j+1]-zeros[j])
            delt = t[zeros[j]+1]-t[zeros[j]]
            sample = int(sample_ref/(delt)) 
#            print(sample)
            n = int(traj_l/sample)
            print(len(range(0,sample*n,sample)))
#            print(n)
            n_calc=traj_l#/delt
#            delt = sample*delt
            
#            varz = np.var([np.average(pos[zeros[j]:zeros[j+1]]) 
#            pos_block = [np.average(pos[zeros[j]:zeros[j]+sample]),np.average(pos[zeros[j]+sample:zeros[j]+2*sample])]
#            pos_block = 0.5*(pos[zeros[j]:zeros[j]+sample]+pos[zeros[j]+sample:zeros[j]+2*sample])
#            pos_block = [np.average(pos[zeros[j]+i:zeros[j]+i+sample-1]) for i in range(sample)]
#            try:
##                varz = np.var([np.average(pos[zeros[j]+i:zeros[j]+1+i,:,2]) for i in range(0,traj_l-1,2)])            
#                varz = np.var(pos[zeros[j]:zeros[j+1],:,2])
#                pos_block = [np.average(pos[zeros[j]+i:zeros[j]+i+sample,:,2]) for i in range(0,sample*n,sample)]
#                ct = [np.where(abs(pos[zeros[j]+i:zeros[j]+i+sample,:,2]-pos[zeros[j]]) < np.var(pos_block),1,0) for i in range(0,sample*n,sample)]
#                ct = np.average(ct,axis=0)
#                t_x = t[zeros[j]:zeros[j]+sample]
#                avz = np.average(pos[zeros[j]:zeros[j+1],:,2])
#                tau = (n_calc*np.var(pos_block)/(varz) - 1)*delt/2 # ps
#                dz.append((avz,varz/tau*1e-8**2/1e-12,varz,tau))
#                ax[1].plot(t_x,ct)

#            except:
#                varz = np.var([np.average(pos[zeros[j]+i:zeros[j]+1+i]) for i in range(0,traj_l-1,2)])            
            varz = np.var(pos[zeros[j]:zeros[j+1]])
            pos_block = [np.average(pos[zeros[j]+i:zeros[j]+i+sample]) for i in range(0,sample*n,sample)]
#            varz_block = np.average([np.var(pos[zeros[j]+i:zeros[j]+i+sample]) for i in range(0,sample*n,sample)])
            varz_block = np.var([np.average(pos[zeros[j]+i:zeros[j]+i+sample]) for i in range(0,sample*n,sample)])

            ct = [np.where(abs(pos[zeros[j]+i:zeros[j]+i+sample]-pos[zeros[j]+i]) < cutoff,1,0) for i in range(0,sample*n,sample)]
            ct = np.average(ct,axis=0)
            ct_collect.append(ct)
            t_x = t[zeros[j]:zeros[j]+sample]
            avz = np.average(pos[zeros[j]:zeros[j+1]])
            tau = (n_calc*np.var(pos_block)/(varz) - 1)*delt/2 # ps

            try:
                parameters=optimize.curve_fit(stretched_e,t_x,ct,p0=[25,1])
                a=parameters[0][0]; b=parameters[0][1]
            except:
                a=1; b=1
#            tau_st = a/b*math.factorial(abs(int(1/b)))
#            print(tau,tau_st)
            dz.append((avz,varz/tau*1e-7**2/1e-12,varz/integrate.trapezoid(ct,t_x)*1e-7**2/1e-12,varz))

            #print('z_av: ',np.average(pos_block),'var(z_): ',np.var(pos_block),'var(z): ',varz,'n: ',n,'delt: ',delt,'traj_l: ',traj_l,'sample: ',sample)
        else:
            print('Length of trajectory insufficient')

    t_x[0] = t_x[1]/100
    ax[1].errorbar(t_x,np.average(np.array(ct_collect),axis=0),yerr=np.std(np.array(ct_collect),axis=0),capsize=2,errorevery=10,label='%s cutoff' % cutoff); ax[1].legend(); ax[1].set_xlabel('Time (ps)',fontsize=28); ax[1].set_ylabel('S(t)',fontsize=28); ax[1].set_ylim(0,1) 
    ax[0].set_ylabel('$D_{ion}$',fontsize=28); ax[0].set_xlabel('Distance from Channel Center (nm)',fontsize=28); ax[1].set_xscale('log')
    dz = np.array(dz)
    dz = dz[dz[:,0].argsort()]
    ax[0].plot(dz[:,0],dz[:,1],linestyle='none',marker='*'); ax[0].set_yscale('log')
    print(dz)
    return dz



def self_auto(t,pos,cutoff,step,sample,ax,zeros):
#    zeros = np.where(t == 0.0)[0]
    print(zeros)
    if len(zeros) == 1:
        print('Zeros len is 1')
        zeros = [0,len(pos)]
#    print(pos[0,:])
#    print(pos[:,0])
    dz = []
    z = []
    cut = []
#    print(len(zeros))
    for j in range(len(zeros)-1):
        sample = zeros[j+1]-zeros[j]
        x_t = np.array([np.linalg.norm(pos[i:zeros[j+1]]-pos[i],axis=2) for i in range(zeros[j],zeros[j+1])])
#        print(x_t)
        c_t = np.zeros((sample,sample,len(g1)))
        c_t[:,:] = np.nan
#        print(t)
        for i in range(0,sample):
            c_t[i,:sample-i,:] = np.where(x_t[i] <= float(cutoff),1,0) 
#        print(c_t)
        c_t = np.nanmean(c_t,axis=0)
        c_t = np.concatenate(c_t)
#        print([np.where(x_t[i] <=float(cutoff),1,0) for i in range(len(pos))])
#        print(np.shape(c_t[:,0]),np.shape(t[zeros[j]:zeros[j+1]]))
#        print(zeros[j],zeros[j+1])
        z_av = np.average(pos[zeros[j]:zeros[j+1],:,2]) 
        ax[0].plot(t[:sample],c_t,label='C(t)')#;ax.plot(t[:sample],stretched_e(t[:sample],10,0.5));plt.savefig('ct.png') #;ax.plot(t[:sample],stretched_e(t[:sample],100,2),label='e^(-t/10000)^0.1'); plt.savefig('c_t_{:.2f}.png'.format(z_av))
#        print(np.var(pos[zeros[j]:zeros[j+1],:,2])**2,integrate.trapezoid(c_t,t[zeros[j]:zeros[j+1]]))
#        print('z= ',np.average(pos[zeros[j]:zeros[j+1],:,2]),'d_z= ',np.var(pos[zeros[j]:zeros[j+1],:,2])**2/integrate.trapezoid(c_t,t[zeros[j]:zeros[j+1]])*1e-8**2/1e-12,'c_t= ',np.transpose(c_t),'cutoff= ',cutoff)
        z.append(np.average(pos[zeros[j]:zeros[j+1],:,2]))
        tau=integrate.trapezoid(c_t,t[zeros[j]:zeros[j+1]])/1e-12
#        var_z=np.var(pos[ze0ros[j]:zeros[j+1]])
#        print(z[-1],var_z/tau)
#        fig, ax = plt.subplot
#        ct = np.convolve(np.roll(c_t,-1),c_t,'same')
#        try:
        try:
            popt, pcov = optimize.curve_fit(stretched_e,t[:sample],c_t,p0=[10,0.5],bounds=[0,1000])
            print(popt,pcov) 
            tau_k, a = popt[:]
        except:
            tau_k, a = dz[np.argmin(abs(np.array(dz)[:,0]-z_av))][2:]

            #try:
            #    slope, intercept, rvalue, pvalue, stderr = scipy.stats.linregress(np.log(t[:sample]),np.log(-np.log(c_t)))
            #    print(slope,intercept)
            #    a, tau_k = np.exp(-intercept/slope), slope
            #except:
            #    try:
            #        tau_k, a = dz[np.argmin(abs(np.concatenate(dz)[:,0]-z_av))][2:]
            #    except:
            #        tau_k, a = 1, 1
            
            
#            ax.plot(t[:sample])
#        except:
#            popt, pcov = [1,1],0
#        print(popt,pcov)
        print(a,tau_k)
        ax[0].plot(t[:sample],stretched_e(t[:sample],tau_k,a)) #; plt.savefig('c_t_{:.2f}_{:.2f}.png'.format(z_av,cutoff))
 #       tau = tau_k/a*math.factorial(int(1/a-1))
        var=np.var(pos[zeros[j]:zeros[j+1],:,2])
        dz.append((z_av,var*1e-7**2/tau,tau_k,a))
        print(cutoff,dz)
        cut.append(cutoff)
        results = np.array([z,dz,cut])

    dz = np.array(dz)
    dz = dz[dz[:,0].argsort()]
    ax[1].plot(dz[:,0],dz[:,1],linestyle='None',marker='*');ax[1].set_yscale('log') #;fig2.savefig('dz.png')
    print(dz[:,0],dz[:,1],np.shape(dz))
#    histo = np.histogram2d(dz[:,0],dz[:,1],bins=50)
#    print(histo)
    return dz

def write_data(dz,header,fileout):
    c = open('{}.dat'.format(fileout),'w')
    c.write(",".join(header)+'\n')
    for i in range(len(dz[:,0])):
        c.write("{},{},{},{}\n".format(dz[i,0],dz[i,1],dz[i,2],dz[i,3]))
    c.close()
    return
    

def main(tpr,trj,step,sample,sel1,sel2,cutoff,fileout):
    if trj[-3:] in ['xtc','trr','tpr']: 
        u=mda.Universe(tpr,trj)
        g1 = u.select_atoms('%s' % sel1)
        g2 = u.select_atoms('%s' % sel2)
        g3 = u.select_atoms('resname UNL and (name O* and bonded name C*)')
        print(g1,g2,g3)
        t,pos,zeros=read_in_data(u,g1,g2,g3,cutoff,step,sample)
    else:
        print('Loading from data file')
        t,pos,zeros=read_in_files(trj)
        print(t,pos,zeros)
    fig, ax = plt.subplots(2,figsize=(11,8.5))
    if float(cutoff) == 0.0:
        print('Iterating Cutoff')
        for cutoff in np.linspace(1,5,9):
            dz=self_auto(t,pos,cutoff,step,sample,ax,zeros)
    elif float(cutoff) > 0.0:
        dz=self_auto(t,pos,cutoff,step,sample,ax,zeros)

    elif float(cutoff) < 0.0:
        dz=self_var_tau(t,pos,cutoff,step,sample,ax,zeros)

    plt.subplots_adjust(left=0.14, right=0.95, top=0.95, bottom=0.14)
#plt.savefig('figures/frp_stl.pdf', format = 'pdf', dpi = 1200)
    plt.savefig('c_t_%s.png' % sample,dpi=1200)
    plt.savefig('c_t_%s.pdf' % sample,dpi=1200,format='pdf')
    write_data(dz,('z','Diffusivity','tau_k','alpha'),fileout)

if __name__=="__main__":
    main(tpr,trj,step,sample,sel1,sel2,cutoff,fileout)



