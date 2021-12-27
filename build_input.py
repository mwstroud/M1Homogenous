import numpy as np
import sys
from bmtools.cli.plugins.util import util
from bmtk.utils.reports.spike_trains import PoissonSpikeGenerator
import random
import csv

#################################################################
####### IMPORTANT: the code here doesnt scale properly! #########
#################################################################

def lognorm_fr_list(n,m,s):
    mean = np.log(m) - 0.5 * np.log((s/m)**2+1)
    std = np.sqrt(np.log((s/m)**2 + 1))
    ranlog= np.random.lognormal(mean,std)
    #print(ranlog)
    return [ranlog for i in range(n)]

#def build_input(t_sim, numPN_A = 640, numPN_C=260, numBask = 100):
def build_poisson_input(population,node_ids,mean,std,output_h5,tstart,tend,t_sim=15000):
    print('Building input for ' + population + "[" + str(len(node_ids)) + " cells at " + str(mean) + "(" + str(std) + ") Hz]")
    psg = PoissonSpikeGenerator(population=population)
    psg.add(node_ids=node_ids,  
    firing_rate=lognorm_fr_list(len(node_ids),mean,std),
    times=(tstart, tend))  
    psg.to_sonata(output_h5)

def populations(config):
    nodes = util.load_nodes_from_config(config)
    CP_nodes = []
    CS_nodes = []
    FSI_nodes = []
    LTS_nodes = []
    for j in nodes:
        node=nodes[j]
        num_cells=len(node['node_type_id'])
        if node['model_type'][0]=='virtual':
            num_cells=len(node['node_type_id'])
            for i in range(num_cells):
                if(node['pop_name'][i]=='CP'):
                    CP_nodes.append(i)
                if(node['pop_name'][i]=='CS'):
                    CS_nodes.append(i)
                if(node['pop_name'][i]=='FSI'):
                    FSI_nodes.append(i)
                if(node['pop_name'][i]=='LTS'):
                    LTS_nodes.append(i)
    return CP_nodes, CS_nodes, FSI_nodes, LTS_nodes 

def thalnodes(config):
    nodes = util.load_nodes_from_config(config)
    Thal_nodes = []
    for j in nodes:
        node=nodes[j]
        num_cells=len(node['node_type_id'])
        if node['model_type'][0]=='virtual':
            num_cells=len(node['node_type_id'])
            for i in range(num_cells):
                if(node['pop_name'][i]=='thal'):
                    Thal_nodes.append(i)
    return Thal_nodes 

def Intnodes(config):
    nodes = util.load_nodes_from_config(config)
    Int_nodes = []
    for j in nodes:
        node=nodes[j]
        num_cells=len(node['node_type_id'])
        if node['model_type'][0]=='virtual':
            num_cells=len(node['node_type_id'])
            for i in range(num_cells):
                if(node['pop_name'][i]=='Int'):
                    Int_nodes.append(i)
    return Int_nodes

def build_input(t_sim, num_thal = 800, num_CS=200, num_CTH = 200, num_CC=200, num_FSI=120, num_LTS=80, scale=1):

    print("Building all input spike trains")
    # M2
    #build_poisson_input(population='M2',
                        #node_ids=range((numPN_A+numPN_C+numPV)*scale),
                        #mean=50,std=1,
                        #output_h5='M2.h5',
                        #t_sim=t_sim)

                       
    Int_nodes=Intnodes("config.json")
    num_int=len(Int_nodes)
    Thal_nodes=thalnodes("config.json") 
    num_thal = len(Thal_nodes)
    
    # Get the thalamic node id that will connect to CP vs CS cells
    CP_nodes = []
    CS_nodes = []
    for i in range(num_thal):
        if  Thal_nodes[i] < 20:
            CP_nodes.append(Thal_nodes[i])
        elif Thal_nodes[i] >= 20 and Thal_nodes[i] < 400:
            CS_nodes.append(Thal_nodes[i])
        elif Thal_nodes[i] >= 400 and Thal_nodes[i] < 780:
            CP_nodes.append(Thal_nodes[i])
        elif Thal_nodes[i] >= 780 and Thal_nodes[i] < 800:
            CS_nodes.append(Thal_nodes[i])
    
    num_CP = len(CP_nodes)
    num_CS = len(CS_nodes)
    #print (num_CP+num_CS)
    # Divide 
    assemblies = 8
    num_group = round(num_thal/assemblies)
    # Shuffle the nodes
    #random.shuffle(Thal_nodes)
    random.shuffle(CP_nodes)
    random.shuffle(CS_nodes)
    #print(Thal_nodes)
    thal1=[]
    thal2=[]
    thal3=[]
    thal4=[]
    thal5=[]
    thal6=[]
    thal7=[]
    thal8=[]
    # Split the CP nodes into 8 assemblies
    for i in range(num_CP):
        if i < round(1*(num_CP/assemblies)):
            thal1.append(CP_nodes[i])
        if i >= round(1*(num_CP/assemblies)) and i < round(2*(num_CP/assemblies)):
            thal2.append(CP_nodes[i])
        if i >= round(2*(num_CP/assemblies)) and i < round(3*(num_CP/assemblies)):
            thal3.append(CP_nodes[i])
        if i >= round(3*(num_CP/assemblies)) and i < round(4*(num_CP/assemblies)):
            thal4.append(CP_nodes[i])
        if i >= round(4*(num_CP/assemblies)) and i < round(5*(num_CP/assemblies)):
            thal5.append(CP_nodes[i])
        if i >= round(5*(num_CP/assemblies)) and i < round(6*(num_CP/assemblies)):
            thal6.append(CP_nodes[i])
        if i >= round(6*(num_CP/assemblies)) and i < round(7*(num_CP/assemblies)):
            thal7.append(CP_nodes[i])
        if i >= round(7*(num_CP/assemblies)) and i < round(8*(num_CP/assemblies)):
            thal8.append(CP_nodes[i])
            
    # Split the CS nodes into 8 assemblies. Append to the same thal lists 
    for i in range(num_CS):
        if i < round(1*(num_CS/assemblies)):
            thal1.append(CS_nodes[i])
        if i >= round(1*(num_CS/assemblies)) and i < round(2*(num_CS/assemblies)):
            thal2.append(CS_nodes[i])
        if i >= round(2*(num_CS/assemblies)) and i < round(3*(num_CS/assemblies)):
            thal3.append(CS_nodes[i])
        if i >= round(3*(num_CS/assemblies)) and i < round(4*(num_CS/assemblies)):
            thal4.append(CS_nodes[i])
        if i >= round(4*(num_CS/assemblies)) and i < round(5*(num_CS/assemblies)):
            thal5.append(CS_nodes[i])
        if i >= round(5*(num_CS/assemblies)) and i < round(6*(num_CS/assemblies)):
            thal6.append(CS_nodes[i])
        if i >= round(6*(num_CS/assemblies)) and i < round(7*(num_CS/assemblies)):
            thal7.append(CS_nodes[i])
        if i >= round(7*(num_CS/assemblies)) and i < round(8*(num_CS/assemblies)):
            thal8.append(CS_nodes[i])


    
    # Short burst input for 1000 ms with 100 ms subbursts followed by 500 ms silence
    Thal=[thal1,thal2,thal3,thal4,thal5,thal6,thal7,thal8]
    #print (len(Thal))
    

    with open("./input/Assembly_ids.csv", "w", newline="") as f:
        writer = csv.writer(f)
        for row in Thal:
            writer.writerow(row)
        
    psgs = PoissonSpikeGenerator(population='thalamus')
    for i in range(8):
        for j in range(assemblies):
            time1=i*1.5 + j*0.125
            time2=i*1.5 + j*0.125 +0.125
            psgs.add(node_ids=Thal[j],  
            firing_rate=lognorm_fr_list(num_group*scale,50,0),
            times=(time1, time2)) 
    # Long burst input for 1000 ms followed by 500 ms silence
    psgl = PoissonSpikeGenerator(population='thalamus') 
    for i in range(assemblies):
        time1=i*1.5
        time2=i*1.5+1
        psgl.add(node_ids=Thal[i],  
            firing_rate=lognorm_fr_list(num_group*scale,50,0),
            times=(time1, time2)) 
        
    #print(len(Thal_nodes))
    # Thalamus test constant 50 Hz input
    psgc = PoissonSpikeGenerator(population='thalamus')
    psgc.add(node_ids=Thal_nodes,  
    firing_rate=lognorm_fr_list(num_thal*scale,50,0),
    times=(0, 10))  
    
    # Thalamus baseline
    psgb = PoissonSpikeGenerator(population='thalamus')
    psgb.add(node_ids=Thal_nodes,  
    firing_rate=lognorm_fr_list(num_thal*scale,2,0),
    times=(0, 10))  
    
    # Interneuron baseline
    psgi = PoissonSpikeGenerator(population='Intbase')
    psgi.add(node_ids=Int_nodes,  
    firing_rate=lognorm_fr_list(num_int*scale,2,0),
    times=(0, 10))

    """
    def random_thal(p=0.1,num_cells=num_thal,nodes=Thal_nodes):
        prob = np.random.uniform(0,1,num_cells)
        pick = []
        i=0
        for choice in prob:
            if choice < p:
                pick.append(nodes[i])
                i=i+1
            else:
                i=i+1
        num_chosen = len(pick)
        return pick , num_chosen

    pick, num_chosen = random_thal()
    # Long burst 1000 ms duration followed by 500 ms silence
    psgl = PoissonSpikeGenerator(population='thalamus')
    psgl.add(node_ids=pick,  
    firing_rate=lognorm_fr_list(num_chosen*scale,50,0),
    times=(0, 1))  

    pick, num_chosen = random_thal()
    
    psgl.add(node_ids=pick,  
    firing_rate=lognorm_fr_list(num_chosen*scale,50,0),
    times=(1.5, 2.5))
    
    pick, num_chosen = random_thal()
    
    psgl.add(node_ids=pick,  
    firing_rate=lognorm_fr_list(num_chosen*scale,50,0),
    times=(3, 4))
    
    pick, num_chosen = random_thal()
    
    psgl.add(node_ids=pick,  
    firing_rate=lognorm_fr_list(num_chosen*scale,50,0),
    times=(4.5, 5.5))
    
    pick, num_chosen = random_thal()
    
    psgl.add(node_ids=pick,  
    firing_rate=lognorm_fr_list(num_chosen*scale,50,0),
    times=(6, 7))
    
    pick, num_chosen = random_thal()
    
    psgl.add(node_ids=pick,  
    firing_rate=lognorm_fr_list(num_chosen*scale,50,0),
    times=(7.5, 8.5))
    
    pick, num_chosen = random_thal()
    
    psgl.add(node_ids=pick,  
    firing_rate=lognorm_fr_list(num_chosen*scale,50,0),
    times=(9, 10))
    
    pick, num_chosen = random_thal()
    
    psgl.add(node_ids=pick,  
    firing_rate=lognorm_fr_list(num_chosen*scale,50,0),
    times=(10.5, 11.5))
    """ 
    
    psgc.to_sonata('./input/thalamus_const.h5')
    psgb.to_sonata('./input/thalamus_base.h5')
    psgs.to_sonata('./input/thalamus_short.h5')
    psgl.to_sonata('./input/thalamus_long.h5')
    psgi.to_sonata('./input/Intbase.h5')

    # These inputs are for the baseline firing rates of the cells in the shell.
    ############ Short Burst Input ##############
    CP_nodes, CS_nodes, FSI_nodes, LTS_nodes = populations("config.json")
    num_CP = len(CP_nodes)
    num_CS = len(CS_nodes)
    num_FSI = len(FSI_nodes)
    num_LTS = len(LTS_nodes)
    # CP Average
    psgl = PoissonSpikeGenerator(population='shell')
    psgl.add(node_ids=CP_nodes,  
    firing_rate=lognorm_fr_list(num_CP*scale,0.5,0.2),
    times=(0, 11.5))  
    psgl.to_sonata('./input/CP_shell_short.h5')

    # CS Average
    psgl = PoissonSpikeGenerator(population='shell')
    psgl.add(node_ids=CS_nodes,  
    firing_rate=lognorm_fr_list(num_CS*scale,0.5,0.2),
    times=(0, 11.5))  
    psgl.to_sonata('./input/CS_shell_short.h5')

    # FSI Average
    psgl = PoissonSpikeGenerator(population='shell')
    psgl.add(node_ids=FSI_nodes,  
    firing_rate=lognorm_fr_list(num_FSI*scale,5,3),
    times=(0, 11.5))
    psgl.to_sonata('./input/FSI_shell_short.h5')

    # LTS Average
    psgl = PoissonSpikeGenerator(population='shell')
    psgl.add(node_ids=LTS_nodes,  
    firing_rate=lognorm_fr_list(num_LTS*scale,1,0.1),
    times=(0, 11.5))
    psgl.to_sonata('./input/LTS_shell_short.h5')

    ######### Long Burst Shell input ##########
    # CP Average
    psgl = PoissonSpikeGenerator(population='shell')
    psgl.add(node_ids=CP_nodes,  
    firing_rate=lognorm_fr_list(num_CP*scale,0.5,0.2),
    times=(0, 11.5))  
    psgl.to_sonata('./input/CP_shell_long.h5')

    # CS Average
    psgl = PoissonSpikeGenerator(population='shell')
    psgl.add(node_ids=CS_nodes,  
    firing_rate=lognorm_fr_list(num_CS*scale,0.5,0.2),
    times=(0, 11.5))  
    psgl.to_sonata('./input/CS_shell_long.h5')

    # FSI Average
    psgl = PoissonSpikeGenerator(population='shell')
    psgl.add(node_ids=FSI_nodes,  
    firing_rate=lognorm_fr_list(num_FSI*scale,5,3),
    times=(0, 11.5))
    psgl.to_sonata('./input/FSI_shell_long.h5')

    # LTS Average
    psgl = PoissonSpikeGenerator(population='shell')
    psgl.add(node_ids=LTS_nodes,  
    firing_rate=lognorm_fr_list(num_LTS*scale,1.5,1),
    times=(0, 11.5))
    psgl.to_sonata('./input/LTS_shell_long.h5')


    print("Done")


if __name__ == '__main__':
    if __file__ != sys.argv[-1]:
        build_input(int(sys.argv[-1]))
    else:
        build_input(10)
