"""Lecture 3 example, computing sqrt with Newton's method
Will be developed further during Lecture 4
"""

def sqrt2(a,x0=1):
    """ Use Newton's method to calculate sqrt(a)
    """

    #Check input
    assert type(a) is int or type(a) is float, "error, input must be numeric"
    assert a>=0, "a must be non-negative"

    tol = 1.0e-12 #convergence criteria
    imax = 10000 #maximum number of iterations

    #Newton's method
    for i in range(imax):
        x1 = x0/2 + a/(2*x0)
        print("iteration %d, x = %16.14f" %(i+1,x1))
        del_x = abs(x1-x0)
        if del_x < tol:
            print("converged!")
            break
        x0 = x1
        

    return x1
