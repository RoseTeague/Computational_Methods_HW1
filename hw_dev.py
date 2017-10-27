"""M3C 2017 Homework 1
Rosemary Teague
00828351
"""
import numpy as np
import matplotlib.pyplot as plt
import networkx as nx

def rw2d(Nt,M,a=0,b=0):
    """Input variables
    Nt: number of time steps
    M: Number of simulations
    a,b: bias parameters
    """
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
        X=X+x
        Y=Y+y
        X2=X2+np.multiply(x,x)
        Y2=Y2+np.multiply(y,y)
        XY=XY+np.multiply(x,y)

    #converts each averaged output to a list, ensures the starting position is at zero and divides by
    #the total number of paths traversed. Stores the values(lists) to be returned in one list.
    output=[X,Y,X2,Y2,XY]
    for i in range(5):
        output[i]=(output[i]/M).tolist()
        output[i][:0]=[0]

    return output


def rwnet1(H,Hf,a=0,display=False):
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

    #assert a==-1,0,1 \
    #    "a should be -1, 0 or +1"
    M=20
    D=np.sqrt(1+(1+a)**2)
    G=nx.Graph()
    G.add_node(0,pos=(0,0))
    #G.add_node(1,pos=(0,H))
    pos=nx.get_node_attributes(G,'pos')
    l=len(pos)

    for m in range(M):
        poslist=rw2d(50,1,a,1)

        X=poslist[0]
        Y=poslist[1]+(np.zeros(51)+H)
        BRK=0

        for i in range(51):
            new_node_pos=(X[i],Y[i])
            for j in np.arange(l-1,-1,-1):
                d=np.sqrt((new_node_pos[0]-pos[j][0])**2+(new_node_pos[1]-pos[j][1])**2)
                if d<=D:
                    G.add_node(l,pos=new_node_pos)
                    BRK=1
                    break
                elif new_node_pos[1]==0:
                    G.add_node(l,pos=new_node_pos)
                    BRK=1
                    break
                else:
                    continue

            if BRK == 0:
                continue
            else:
                break

        pos=nx.get_node_attributes(G,'pos')
        l=len(pos)

    return G,pos,X,Y#,output

def rwnet2(L,H,Hf,a=0,display=False):
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

    return X,Y,output

def analyze():
    """ Add input variables as needed
    """

def network(X,Y,dstar,display=False,degree=False):
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



if __name__ == '__main__':
    #The code here should call analyze and generate the
    #figures that you are submitting with your code
    analyze()
