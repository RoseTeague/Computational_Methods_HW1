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
        l=-1, d=0; l=0, d=1; l=1, d=0.

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
    analyze: data corresponding to rate of growth of network collected when true

    Output variables
    X,Y: Final node coordinates
    output variables
    X : list of x-coordinated of nodes in the network
    Y : list of x-coordinated of nodes in the network
    Ym : list of y coordinates corresponding to nodes which increased the overall
        height of the network
    Tm : list of iteration numbers corresponding to points the overall height of
        the network was increased
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

    # create random step weighted by b and a=1 and starting at (0,H)
    # network growth will stop after reaching a height of Hf
    while max(y)<Hf:
        X=X+np.random.randint(0,2)*(-2-a)+(1+a)
        Y=Y-np.random.randint(0,2)

        # Nested loops to check if Y has reached 0, else if (X,Y) is within d=sqrt(1+(1+a)^2)
        # of any existing node. If either of these conditions are met then a BRK flag is set to one.
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

        # When BRK flag triggered, (X,Y) to be appended to the lists, x,y, of existing nodes.
        # If analyze is set to true, the new node will be tested to check if it is higher than
        # any other existing node, and if it is, it will be appended to a new list, ymax and
        # its iteration number will be appended to tmax.
        if BRK==1:
            if analyze:
                t=t+1
                if Y>max(y):
                    ymax.append(Y)
                    tmax.append(t)
            x.append(X);y.append(Y)
            X=0;Y=H;BRK=0

    # displays and saves a plot of node positions if display is set to True
    if display:
        n='H='+str(H)+',Hf='+str(Hf)+',a='+str(a)
        plt.plot(x,y,'ro',markersize=1)
        plt.title('Rosemary Teague, rwnet1 \n'+n)
        plt.savefig(n,dpi=700)
        plt.show()

    # end = time.clock()
    # print(end-start)

    return x,y,ymax,tmax

def rwnet2(L,H,Hf,a=0,display=False,analyze=False):
    """Input variables
    L: Walls are placed at X=+/-L
    H: Height at which new nodes are initially introduced
    Hf: Final network height
    a: horizontal bias parameter, should be -1, 0, or 1
    display: figure displaying the network is created when true
    analyze: data corresponding to rate of growth of network collected when true

    Output variables
    X,Y: Final node coordinates
    output variables
    X : list of x-coordinated of nodes in the network
    Y : list of x-coordinated of nodes in the network
    Ym : list of y coordinates corresponding to nodes which increased the overall
        height of the network
    Tm : list of iteration numbers corresponding to points the overall height of
        the network was increased
    """
    # start = time.clock()
    assert a in [-1,0,1], \
        "a should be -1, 0 or 1"

    #prepare figure and set up initial parameters
    plt.close()
    x=[0];y=[0];t=0;X=0;Y=H;BRK=0
    ymax=[]
    tmax=[]

    # create random step weighted by b and a=1 and starting at (0,H)
    # network growth will stop after reaching a height of Hf
    while max(y)<Hf:
        X=X+np.random.randint(0,2)*(-2-a)+(1+a)
        Y=Y-np.random.randint(0,2)
        #If X is greater than or equal to the wall position, it will be altered to be set to L.
        if X>=L:
            X=L
        elif X<=-L:
            X=-L

        # Nested loops to check if Y has reached 0, else if (X,Y) is within d=sqrt(1+(1+a)^2)
        # of any existing node. If either of these conditions are met then a BRK flag is set to one.
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

        # When BRK flag triggered, (X,Y) to be appended to the lists, x,y, of existing nodes.
        # If analyze is set to true, the new node will be tested to check if it is higher than
        # any other existing node, and if it is, it will be appended to a new list, ymax and
        # its iteration number will be appended to tmax.
        if BRK==1:
            if analyze==True:
                t=t+1
                if Y>max(y):
                    ymax.append(Y)
                    tmax.append(t)
            x.append(X);y.append(Y)
            X=0;Y=H;BRK=0

    # displays and saves a plot of node positions if display is set to True
    if display==True:
        n='H='+str(H)+',Hf='+str(Hf)+',a='+str(a)+',L='+str(L)
        plt.plot(x,y,'ro',markersize=1)
        plt.title('Rosemary Teague, rwnet2 \n'+n)
        plt.savefig(n,dpi=700)
        plt.show()

    # end = time.clock()
    # print(end-start)

    return x,y,ymax,tmax

def plotting(x1,y1,x2,y2,H,Hf,L=0):
    """Input variables
    x1,y1 : list data for first plot
    x2,y2 : list data for second plot to be superimposed on the first
    L : Parameter for wall distance to be included in plot title
    """

    if L==0:
        n='H='+str(H)+',Hf='+str(Hf)+',L=Infinity'
    else:
        n='H='+str(H)+',Hf='+str(Hf)+',L='+str(L)

    one,=plt.plot(x1,y1,'b-')
    two,=plt.plot(x2,y2,'g-')
    one.set_label('alpha=0')
    two.set_label('alpha=1')
    plt.xlabel('Time (Iteration Number)')
    plt.ylabel('Network height')
    plt.title('Rosemary Teague, Analyze \n'+n)
    plt.legend()


def analyze(H,Hf,L1=30,L2=150):
    """ Input variables
    H : Height at which new nodes are initially introduced
    Hf : Final network height
    L1 : Position of walls for first comparison
    L2 : Position of walls for second comparison

    Y of highest node plotted against the iteration it corresponds to for the following cases:

    Figure 1: L=infinity, uses input from rwnet1, plots the trends for when a=0 and a=1 as comparison.
        The case for a=1 should show a more rapid growth (steeper curve) as new nodes will practically
        always arrive from the same side, causing the network to have an angle. In contrast, when a=0
        there is an equal chance for nodes to be added either side of the starting growth and so the
        network will grow outwards more than when a=1 and hence the height will not increase so rapidly.

    Figure 2: L1=30, uses input from rwnet2, plots the trends for when a=0 and a=1 as comparison.
        In this case, L<<Hf and so most of the graphs height will be a result of nodes
        'sticking' to the walls. When a=0, we will still expect the majority of the growth to be upwards
        as paths are made in a general downwards direction, however it will become truncated for low
        enough Hf. When a=1 however, the paths have a strong bias to the left and so very rapidly encounter
        the walls. Due to this strong bias they stay fairly close to the wall until reaching y=0 or
        falling within d of an existing node. This causes the height of the network to increase very
        rapidly and growth to stop after a comparitively short time.

    Figure 3: L2=150, uses input from rwnet2, plots the trends for when a=0 and a-1 as comparison.
        In this case, as in the previous cases, when a=0 the path has no left or right bias and will
        grow vertically at a fairly steady rate. However, when a=1, the nodes will tend to gather along
        the right hand wall and an initial, sharp increase in network height is expected. However,
        as L=Hf, there is a higher chance of a path falling more before reaching the wall, and hence
        more of an outward (towards the centre) growth is expected. As this builds up the nodes will
        encounter the existing network before reaching the wall and the growth rate will tend towards
        that seen for L=infinity. i.e, a change in rate will be observed from more steep to more shallow.
    """

    #generates x and y coordinates for distributions to be tested using rwnet1
    x,y,ymax0i,tmax0i=rwnet1(H,Hf,0,analyze=True)
    x,y,ymax1i,tmax1i=rwnet1(H,Hf,1,analyze=True)

    #plots the height vs time graph for the network and saves the figure
    fig1=plotting(tmax0i,ymax0i,tmax1i,ymax1i,H,Hf)
    plt.savefig('hw11.png',dpi=700)
    plt.close()

    #generates x and y coordinates for distributions to be tested using rwnet2
    x,y,ymax0j,tmax0j=rwnet2(L1,H,Hf,0,analyze = True)
    x,y,ymax1j,tmax1j=rwnet2(L1,H,Hf,1,analyze = True)

    #plots the height vs time graph for the network and saves the figure
    fig2=plotting(tmax0j,ymax0j,tmax1j,ymax1j,H,Hf,L1)
    plt.savefig('hw12.png',dpi=700)
    plt.close()

    #generates x and y coordinates for distributions to be tested using rwnet2
    x,y,ymax0k,tmax0k=rwnet2(L2,H,Hf,0,analyze = True)
    x,y,ymax1k,tmax1k=rwnet2(L2,H,Hf,1,analyze = True)

    #plots the height vs time graph for the network and saves the figure
    fig3=plotting(tmax0k,ymax0k,tmax1k,ymax1k,H,Hf,L2)
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

    This network has a limit on the maximum degree as only neighbouring nodes (within dstar)
        will link to form edges. If a network with a significantly larger maximum degree was
        required, the network could be modified such that any new node will link to all
        previous successively linked nodes. This would not require any modifications to
        the rwnet functions but would require modifications to this function such that edges
        are added not only to each value in l[j] but to all values in l.
        Alternatively, simply increasing the value of d would increase the number of possible
        connections between nodes and would require no significant modifications to any
        part of the code

    If display is True, this network will be drawn
    If Degree is true, the degree of the network will be calculated and plotted as a bar chart
        to best illustrate the number of nodes with each degree.
    """

    # initialises and empty network
    plt.close()
    l=[]
    g=nx.Graph()

    # compares (x,y) coordinates of nodes. For each pair, if any coordinate at a higher position
    # in the list lies within a distance of dstar, the corresponding location in the coordinates
    # list is appended to 'l'
    for j in range(len(y)):
        links=[i for i in np.arange(j+1,len(y)) if np.sqrt(np.abs(x[i]-x[j])**2+np.abs(y[i]-y[j])**2)<=dstar]
        l.append(links)

    # For each node, if it does not lie within dstar it will be disregarded, otherwise it will be
    # added to the network and linked to all nodes listed in l.
    for j in range(len(l)):
        if len(l[j])>0:
            g.add_node(j)
            for i in range(len(l[j])):
                g.add_edge(j,l[j][i])

    if display:
        nx.draw_networkx(g,node_size=4,with_labels=False)
        plt.title('Rosemary Teague, Network \n Network of Points')
        plt.savefig('Network of points',dpi=700)
        plt.show()

    if degree:
        D=nx.degree(g)
        plt.plot(nx.degree_histogram(g),'k-')
        plt.bar(np.arange(len(nx.degree_histogram(g))),nx.degree_histogram(g))
        plt.title('Rosemary Teague, Network \n Network')
        plt.xlabel('Degree')
        plt.ylabel('Number of nodes')
        plt.savefig('Degree of nodes', dpi=700)
        plt.show()
    else:
        D=[]

    return g,D


if __name__ == '__main__':
    analyze(200,150)
