"""M3C 2017 Homework 1
Rosemary Teague
00828351
"""
import numpy as np
import matplotlib.pyplot as plt
import networkx as nx
import time
import matplotlib.animation as animation

def rw2d(Nt,M,a=0,b=0):
    """Input variables
    Nt: number of time steps
    M: Number of simulations
    a,b: bias parameters
    """

    start = time.time()

    #assert conditions on the values of a and b
    assert a>=-1, \
        "a should not be less than -1"
    assert b<=1, \
        "b should not be greater than 1"

    #set range of simulations
    m=np.arange(1,M)

    # m random walks, each of length Nt, created by setting up 2 arrays of random choices, 0 or 1,
    # which are then weighted by the values a and b and added cumulatively. Each time step is added
    # to the previous simulations in order that a division by M will provide an average
    x0=np.cumsum(np.random.randint(0,2,Nt)*(-2-a)+(1+a))
    y0=y=np.cumsum(np.random.randint(0,2,Nt)*(-2+b)+(1-b))
    for i in m:
        x=np.cumsum(np.random.randint(0,2,Nt)*(-2-a)+(1+a))
        y=np.cumsum(np.random.randint(0,2,Nt)*(-2+b)+(1-b))
        x0=x0+x
        y0=y0+y

    # all outputs are calculated and collected into a list,
    # whereby they are each averaged (divided through by M), converted
    # from arrays to lists, and initial positions at (0,0) are prepended.
    X2=np.multiply(x0,x0)
    Y2=np.multiply(y0,y0)
    XY=np.multiply(x0,y0)
    output=[x,y,X2,Y2,XY]
    for i in range(5):
        output[i]=(output[i]/M).tolist()
        output[i][:0]=[0]

    end = time.time()
    #print(end-start)

    return output


def test(n,l,d,x,y,X,Y,BRK):
    """
    A testing function to be used within network generation to compare distances
        between a new node and all existing points.
    Input variable:
    n : number of occurences of Y+l within the existing network, y.
    l : disatnce from test node, Y, an existing node could be in order for Y to be
        added to the network.
    d : distance from test node, X, an existing node could be in order for X to be
        added to the network, given the value of l.

        e.g if the overall disatnce between a new node,(X,Y), and any existing
        node (x,y), must be less than sqrt(1) to be added to the network, the
        following combinations of l and d are allowed:
        l=-1, d=0; l=0, d=1, l=1, d=0.

    x,y : list of existing node coordinates
    X,Y : new node coordinates to be tested.
    BRK : Flag to check whether node should be added to network

    Output variables
    BRK : Flag will be set to 1 if the node is to be added to the network and 0 otherwise
    """
    if BRK==0:
        y1=-1
        for i in range(n):
            y1=y.index(Y+l,y1+1)
            if (np.abs(X-x[y1]))<=d:
                BRK=1
                break
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
    # start = time.clock()
    #assertions on the value of a
    assert a in [-1,0,1], \
        "a should be -1, 0 or 1"

    #prepare figure and set up initial parameters
    plt.close()
    x=[0];y=[0];t=0;X=0;Y=H;BRK=0
    ymax=[]
    tmax=[]

    while max(y)<Hf:
        X=X+np.random.randint(0,2)*(-2-a)+(1+a)
        Y=Y-np.random.randint(0,2)
        if Y<=0:
            BRK=1
        else:
            if a==0:
                n_1,n0,n1=y.count(Y-1),y.count(Y),y.count(Y+1)
                BRK=test(n_1,-1,1,x,y,X,Y,BRK)
                BRK=test(n0,0,1,x,y,X,Y,BRK)
                BRK=test(n1,1,1,x,y,X,Y,BRK)
            elif a==-1:
                n_1,n0,n1=y.count(Y-1),y.count(Y),y.count(Y+1)
                BRK=test(n_1,-1,0,x,y,X,Y,BRK)
                BRK=test(n0,0,1,x,y,X,Y,BRK)
                BRK=test(n1,1,0,x,y,X,Y,BRK)
            elif a==1:
                n_2,n_1,n0,n1,n2=y.count(Y-2),y.count(Y-1),y.count(Y),y.count(Y+1),y.count(Y+2)
                BRK=test(n_2,-2,1,x,y,X,Y,BRK)
                BRK=test(n_1,-1,2,x,y,X,Y,BRK)
                BRK=test(n0,0,2,x,y,X,Y,BRK)
                BRK=test(n1,1,2,x,y,X,Y,BRK)
                BRK=test(n2,2,1,x,y,X,Y,BRK)
        if BRK==1:
            if analyze==True:
                t=t+1
                if Y>max(y):
                    ymax.append(Y)
                    tmax.append(t)
            x.append(X);y.append(Y)
            X=0;Y=H;BRK=0

    if display==True:
        n='H='+str(H)+',Hf='+str(Hf)+',a='+str(a)
        plt.plot(x,y,'ro',markersize=1)
        plt.title(n)
        plt.savefig(n,dpi=700)
        plt.show()

    # end = time.clock()
    # print(end-start)

    return x,y,ymax,tmax

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
    # start = time.clock()
    assert a in [-1,0,1], \
        "a should be -1, 0 or 1"

    plt.close()
    x=[0];y=[0];t=0;X=0;Y=H;BRK=0
    ymax=[]
    tmax=[]

    while max(y)<Hf:
        X=X+np.random.randint(0,2)*(-2-a)+(1+a)
        Y=Y-np.random.randint(0,2)
        if X>=L:
            X=L
        elif X<=-L:
            X=-L

        if Y<=0:
            BRK=1
        else:
            if a==0:
                n_1,n0,n1=y.count(Y-1),y.count(Y),y.count(Y+1)
                BRK=test(n_1,-1,1,x,y,X,Y,BRK)
                BRK=test(n0,0,1,x,y,X,Y,BRK)
                BRK=test(n1,1,1,x,y,X,Y,BRK)
            elif a==-1:
                n_1,n0,n1=y.count(Y-1),y.count(Y),y.count(Y+1)
                BRK=test(n_1,-1,0,x,y,X,Y,BRK)
                BRK=test(n0,0,1,x,y,X,Y,BRK)
                BRK=test(n1,1,0,x,y,X,Y,BRK)
            elif a==1:
                n_2,n_1,n0,n1,n2=y.count(Y-2),y.count(Y-1),y.count(Y),y.count(Y+1),y.count(Y+2)
                BRK=test(n_2,-2,1,x,y,X,Y,BRK)
                BRK=test(n_1,-1,2,x,y,X,Y,BRK)
                BRK=test(n0,0,2,x,y,X,Y,BRK)
                BRK=test(n1,1,2,x,y,X,Y,BRK)
                BRK=test(n2,2,1,x,y,X,Y,BRK)
        if BRK==1:
            if analyze==True:
                t=t+1
                if Y>max(y):
                    ymax.append(Y)
                    tmax.append(t)
            x.append(X);y.append(Y)
            X=0;Y=H;BRK=0

    if display==True:
        n='H='+str(H)+',Hf='+str(Hf)+',a='+str(a)+',L='+str(L)
        plt.plot(x,y,'ro',markersize=1)
        plt.title(n)
        plt.savefig(n,dpi=700)
        plt.show()

    # end = time.clock()
    # print(end-start)

    return x,y,ymax,tmax

def plotting(x1,y1,x2,y2,L=0):
    """Input variables
    x1,y1 : list data for first plot
    x2,y2 : list data for second plot to be superimposed on the first
    L : Parameter for wall distance to be included in plot title
    """

    if L==0:
        n='L=Infinity'
    else:
        n='L='+str(L)

    one,=plt.plot(x1,y1,'b-')
    two,=plt.plot(x2,y2,'g-')
    one.set_label('alpha=0')
    two.set_label('alpha=1')
    plt.xlabel('Time (Iteration Number)')
    plt.ylabel('Network height')
    plt.title(n)
    plt.legend()


def analyze():
    """ No inmput variables required for analasis.
    Y of highest node plotted against the iteration it corresponds to for the following cases:
    Figure 1: L=infinity uses input from rwnet1, plots the trends for when a=0 and a=1 as comparison
    Figure 2: L=30 uses input from rwnet2, plots the trends for when a=0 and a=1 as comparison.
        Note that in this case, L<<Hf and so most of the graphs height will be a result of nodes 'sticking' to the Walls.
        For the case when a=1, we therefore expect the height to build up very rapidly.
        For the case when a=0, we expect a slower increase as the height is centered around x=0
    Figure 3: L=150 uses input from rwnet2, plots the trends for when a=0 and a-1 as comparison.
        Note in this case L=Hf and so most of the height from the graph will be a result of nodes sticking to other nodes.
    """

    x,y,ymax0i,tmax0i=rwnet1(200,150,0,analyze=True)
    x,y,ymax1i,tmax1i=rwnet1(200,150,1,analyze=True)

    fig1=plotting(tmax0i,ymax0i,tmax1i,ymax1i)
    plt.savefig('hw11.png',dpi=700)
    plt.close()

    x,y,ymax0j,tmax0j=rwnet2(30,200,150,0,analyze = True)
    x,y,ymax1j,tmax1j=rwnet2(30,200,150,1,analyze = True)

    fig2=plotting(tmax0j,ymax0j,tmax1j,ymax1j,30)
    plt.savefig('hw12.png',dpi=700)
    plt.close()

    x,y,ymax0k,tmax0k=rwnet2(150,200,150,0,analyze = True)
    x,y,ymax1k,tmax1k=rwnet2(150,200,150,1,analyze = True)

    fig3=plotting(tmax0k,ymax0k,tmax1k,ymax1k,150)
    plt.savefig('hw13.png',dpi=700)
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

    Compares each set of node positions (x[j],y[j]) with every other node (x[i],y[i])
        if i>j. If the distance between these any two such nodes is smaller than or
        equal to dstar, the nodes will be added to a list, l. This list will then a sub-list
        for each possible node position containing indices of all other nodes it will
        link to. If this list is empty, it signifies the node is not linked to any others and
        will thus be disregarded. Else it will be added to the graph g, and edges will be added
        to all connecting nodes.

    If display is True, this network will be drawn
    If Degree is true, the degree of the network will be calculated and plotted as a bar chart
        to best illustrate the number of nodes with each degree.
    """
    plt.close()

    l=[]
    g=nx.Graph()

    for j in range(len(y)):
        links=[i for i in np.arange(j+1,len(y)) if np.sqrt(np.abs(x[i]-x[j])**2+np.abs(y[i]-y[j])**2)<=dstar]
        l.append(links)

    for j in range(len(l)):
        if len(l[j])>0:
            g.add_node(j)
            for i in range(len(l[j])):
                g.add_edge(j,l[j][i])

    if display:
        nx.draw_networkx(g,node_size=4,with_labels=False)
        plt.title('Connections Network')
        plt.savefig('Network of points',dpi=700)
        plt.show()

    if degree:
        D=nx.degree(g)
        plt.plot(nx.degree_histogram(g),'k-')
        plt.bar(np.arange(len(nx.degree_histogram(g))),nx.degree_histogram(g))
        plt.savefig('Degree of nodes', dpi=700)
        plt.show()
    else:
        D=[]

    return g,D


if __name__ == '__main__':
    #The code here should call analyze and generate the
    #figures that you are submitting with your code
    analyze()
