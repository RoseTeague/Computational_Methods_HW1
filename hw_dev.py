"""M3C 2017 Homework 1
Rosemary Teague
00828351
"""
import numpy as np
import matplotlib.pyplot as plt
import networkx as nx
import time

def rw2d(Nt,M,a=0,b=0):
    """Input variables
    Nt: number of time steps
    M: Number of simulations
    a,b: bias parameters
    """

    start = time.time()

    assert a>=-1, \
        "a should not be less than -1"
    assert b<=1, \
        "b should not be greater than 1"

    X=0;Y=0;X2=0;Y2=0;XY=0
    m=np.arange(1,M+1)

    for i in m:
        Choice=np.random.randint(0,2,Nt)
        Choice2=np.random.randint(0,2,Nt)#Makes a random choice of 0 or 1 for each of Nt time steps
        #Assign the number 0 to a step of 1+a and the number 1 to a step of 1 in X
        # and similarly for steps in Y. Cumulatively add up these steps to create random paths
        x=np.cumsum(Choice*(-1)+(1-Choice)*(1+a))
        y=np.cumsum(Choice2*(-1)+(1-Choice2)*(1-b))
        #Adds each successive path, M, at each timestep. All in the format of Nt dimensional arrays
        x=np.insert(x,0,0)
        y=np.insert(y,0,0)
        X2=np.multiply(x,x)
        Y2=np.multiply(y,y)
        XY=np.multiply(x,y)

    #converts each averaged output to a list, ensures the starting position is at zero and divides by
    #the total number of paths traversed. Stores the values(lists) to be returned in one list.
    output=[x,y,X2,Y2,XY]
    for i in range(5):
        output[i]=(output[i]/M).tolist()
        output[i][:0]=[0]

    end = time.time()
    #print(end-start)

    return output


def test(n,l,d,y,x,X,Y,BRK):
    if BRK==0:
        y1=-1
        for i in range(n):
            y1=y.index(Y+l,y1+1)
            if (np.abs(X-x[y1]))<=d:
                BRK=1
    return BRK


def rwnet1(H,Hf,a=0,display=False,analyze=False):
    """Input variables
    H: Height at which new nodes are initially introduced
    Hf: Final network height
    a: horizontal bias parameter, should be -1, 0, or 1
    display: figure displaying the network is created when true
    Output variables
    X,Y: Final node coordinates
    output: a tuple containing any other information you would
    like the function to return, may be left empty
    """
    #start = time.time()
    assert a in [-1,0,1], \
        "a should be -1, 0 or 1"

    M=H*10
    D=np.sqrt(1+(1+a)**2)
    G=nx.Graph()
    G.add_node(0,pos=(0,0))
    pos=nx.get_node_attributes(G,'pos')
    l=len(pos)
    x=[0]
    y=[0]
    ymax=[]
    tmax=[]

    for m in range(M):
        poslist=rw2d(M,1,a,1)

        X=poslist[0]
        Y=[yi+H for yi in poslist[1]]
        BRK=0

        for i in range(M+1):
            if Y[i]==0:
                BRK=1
            else:
                if a==0:
                    n_1,n0,n1=y.count(Y[i]-1),y.count(Y[i]),y.count(Y[i]+1)
                    BRK=test(n_1,-1,1,y,x,X[i],Y[i],BRK)
                    BRK=test(n0,0,1,y,x,X[i],Y[i],BRK)
                    BRK=test(n1,1,1,y,x,X[i],Y[i],BRK)
                elif a==-1:
                    n_1,n0,n1=y.count(Y[i]-1),y.count(Y[i]),y.count(Y[i]+1)
                    BRK=test(n_1,-1,0,y,x,X[i],Y[i],BRK)
                    BRK=test(n0,0,1,y,x,X[i],Y[i],BRK)
                    BRK=test(n1,1,0,y,x,X[i],Y[i],BRK)
                elif a==1:
                    n_2,n_1,n0,n1,n2=y.count(Y[i]-2),y.count(Y[i]-1),y.count(Y[i]),y.count(Y[i]+1),y.count(Y[i]+2)
                    BRK=test(n_2,-2,1,y,x,X[i],Y[i],BRK)
                    BRK=test(n_1,-1,2,y,x,X[i],Y[i],BRK)
                    BRK=test(n0,0,2,y,x,X[i],Y[i],BRK)
                    BRK=test(n1,1,2,y,x,X[i],Y[i],BRK)
                    BRK=test(n2,2,1,y,x,X[i],Y[i],BRK)

            if BRK==1:
                G.add_node(l,pos=(X[i],Y[i]))
                if analyze==True:
                    if Y[i]>max(y):
                        ymax.append(Y[i])
                        tmax.append(l)
                x.append(X[i])
                y.append(Y[i])
                break
            else:
                continue

        if y[-1]>=Hf:
            G.remove_node(l)
            break

        pos=nx.get_node_attributes(G,'pos')
        l=len(pos)

    if display==True:
        n='H='+str(H)+',Hf='+str(Hf)+',a='+str(a)
        nx.draw_networkx(G,pos,node_size=6,with_labels=False)
        plt.title(n)
        plt.savefig(n)
        plt.close()

    #end = time.time()
    #print(end-start)

    return G,pos,x,y,ymax,tmax#,output

def rwnet2(L,H,Hf,a=0,display=False,analyze=False):
    """Input variables
    L: Walls are placed at X = +/- L
    H: Height at which new nodes are initially introduced
    Hf: Final network height
    a: horizontal bias parameter, should be -1, 0, or 1
    display: figure displaying the network is created when true
    Output variables
    X,Y: Final node coordinates
    output: a tuple containing any other information you would
    like the function to return, may be left empty
    """
    #start = time.time()
    assert a in [-1,0,1], \
        "a should be -1, 0 or 1"

    M=H*10
    D=np.sqrt(1+(1+a)**2)
    G=nx.Graph()
    G.add_node(0,pos=(0,0))
    pos=nx.get_node_attributes(G,'pos')
    l=len(pos)
    x=[0]
    y=[0]
    ymax=[]
    tmax=[]

    for m in range(M):
        poslist=rw2d(M,1,a,1)

        X=poslist[0]
        Y=[yi+H for yi in poslist[1]]
        BRK=0

        for i in range(M+1):
            if Y[i]<=0:
                BRK=1
            elif np.abs(X[i])>=L:
                BRK=1
            else:
                if a==0:
                    n_1,n0,n1=y.count(Y[i]-1),y.count(Y[i]),y.count(Y[i]+1)
                    BRK=test(n_1,-1,1,y,x,X[i],Y[i],BRK)
                    BRK=test(n0,0,1,y,x,X[i],Y[i],BRK)
                    BRK=test(n1,1,1,y,x,X[i],Y[i],BRK)
                elif a==-1:
                    n_1,n0,n1=y.count(Y[i]-1),y.count(Y[i]),y.count(Y[i]+1)
                    BRK=test(n_1,-1,0,y,x,X[i],Y[i],BRK)
                    BRK=test(n0,0,1,y,x,X[i],Y[i],BRK)
                    BRK=test(n1,1,0,y,x,X[i],Y[i],BRK)
                elif a==1:
                    n_2,n_1,n0,n1,n2=y.count(Y[i]-2),y.count(Y[i]-1),y.count(Y[i]),y.count(Y[i]+1),y.count(Y[i]+2)
                    BRK=test(n_2,-2,1,y,x,X[i],Y[i],BRK)
                    BRK=test(n_1,-1,2,y,x,X[i],Y[i],BRK)
                    BRK=test(n0,0,2,y,x,X[i],Y[i],BRK)
                    BRK=test(n1,1,2,y,x,X[i],Y[i],BRK)
                    BRK=test(n2,2,1,y,x,X[i],Y[i],BRK)

            if BRK==1:
                G.add_node(l,pos=(X[i],Y[i]))
                if analyze==True:
                    if Y[i]>max(y):
                        ymax.append(Y[i])
                        tmax.append(l)
                x.append(X[i])
                y.append(Y[i])
                break
            else:
                continue

        if y[-1]>=Hf:
            G.remove_node(l)
            break

        pos=nx.get_node_attributes(G,'pos')
        l=len(pos)

    if display==True:
        n='H='+str(H)+',Hf='+str(Hf)+',a='+str(a)+',L='+str(L)
        nx.draw_networkx(G,pos,node_size=6,with_labels=False)
        plt.title(n)
        plt.savefig(n)
        plt.close()

    #end = time.time()
    #print(end-start)

    return G,pos,x,y,ymax,tmax

def plotting(x1,y1,x2,y2,L=0):

    if L==0:
        n='L=Infinity'
    else:
        n='L='+str(L)

    one,=plt.plot(x1,y1,'b-')
    two,=plt.plot(x2,y2,'g-')
    one.set_label('alpha=0')
    two.set_label('alpha=1')
    plt.xlabel('time')
    plt.ylabel('network height')
    plt.title(n)
    plt.legend()


def analyze():
    """ No inmput variables required for analasis.
    Y of highest node plotted against the iteration it corresponds to for the following cases:
    Figure 1: L=infinity uses input from rwnet1, plots the trends for when a=0 and a=1 as comparison
    Figure 2: L=30 uses input from rwnet2, plots the trends for when a=0 and a=1 as comparison.
        Note that in this case, L<<Hf and so most of the graphs height will be a result of nodes 'sticking' to the Walls.
        For the case when a=0 we don't expect any height as the
    Figure 3: L=150 uses input from rwnet2, plots the trends for when a=0 and a-1 as comparison.
        Note in this case L=Hf and so most of the height from the graph will be a result of nodes sticking to other nodes.
    """

    G,pos,x,y,ymax0i,tmax0i=rwnet1(200,150,0,display=True,analyze=True)
    G,pos,x,y,ymax1i,tmax1i=rwnet1(200,150,1,display=True,analyze=True)

    fig1=plotting(tmax0i,ymax0i,tmax1i,ymax1i)
    plt.savefig('hw11.png')
    plt.close()

    G,pos,x,y,ymax0j,tmax0j=rwnet2(30,200,150,0,display=True,analyze = True)
    G,pos,x,y,ymax1j,tmax1j=rwnet2(30,200,150,1,display=True,analyze = True)

    fig2=plotting(tmax0j,ymax0j,tmax1j,ymax1j,30)
    plt.savefig('hw12.png')
    plt.close()

    G,pos,x,y,ymax0k,tmax0k=rwnet2(150,200,150,0,display=True,analyze = True)
    G,pos,x,y,ymax1k,tmax1k=rwnet2(150,200,150,1,display=True,analyze = True)

    fig3=plotting(tmax0k,ymax0k,tmax1k,ymax1k,150)
    plt.savefig('hw13.png')
    plt.close()

    return


def network(x,y,dstar,display=False,degree=False):
    """ Input variables
    X,Y: Numpy arrays containing network node coordinates
    dstar2: Links are placed between nodes within a distance, d<=dstar of each other
        and dstar2 = dstar*dstar
    display: Draw graph when true
    degree: Compute, display and return degree distribution for graph when true
    Output variables:
    G: NetworkX graph corresponding to X,Y,dstar
    D: degree distribution, only returned when degree is true
    """
    l=[]
    g=nx.Graph()

    dstar=np.sqrt(2)
    for j in range(len(y)):
        links=[i for i in np.arange(j+1,len(y)) if np.sqrt(np.abs(x[i]-x[j])**2+np.abs(y[i]-y[j])**2)<=dstar]
        l.append(links)
    #print(l)

    for j in range(len(l)):
        if len(l[j])>0:
            g.add_node(j)
            for i in range(len(l[j])):
                #print('j=',j,'lj=',l[j],'i=',i,'lij=',l[j][i])
                g.add_edge(j,l[j][i])
    #print(len(g.node))
    #plt.plot(nx.degree_histogram(g))
    nx.draw_networkx(g,node_size=6,with_labels=False)
    plt.title('Connections Network')
    plt.savefig('Network of points')
    plt.show()



if __name__ == '__main__':
    #The code here should call analyze and generate the
    #figures that you are submitting with your code
    analyze()
