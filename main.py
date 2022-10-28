
import json
# from scipy.sparse import csr_matrix
import numpy as np

"""
LINELEM     Array of Linear elements (each row corresponds to one element)
NLNELEM     Nonlinear elements (in our case, it will be the list of MOSFETs)
INFO        Auxiliary information on the simulation
NODES       Actual numbers ("names") of the nodes as given in the circuit
LNAME       Names of the linear elements as given in the circuit
NNAME       Names of the nonlinear elements as given in the circuit
PRINTNV     Node voltages to be printed
PRINTBV     Branch voltages to be printed
PRINTBI     Branch currents to be printed
PLOTNV      Node voltages to be plotted
PLOTBV      Branch voltages to be plotted
PLOTBI      Branch currents to be plotted
"""

class linelem():
    ElemType_C = 67
    ElemType_L = 76
    ElemType_I = 73
    ElemType_R = 82
    ElemType_V = 86
    V_TYPE_DC = 0
    def __init__(self, LINELEM_ARRAY):
        self.ElemType = LINELEM_ARRAY[0]
        self.Init_Value = LINELEM_ARRAY[1]
        self.Node1 = LINELEM_ARRAY[2]
        self.Node2 = LINELEM_ARRAY[3]
        if self.ElemType==linelem.ElemType_R:
            pass
        elif self.ElemType==linelem.ElemType_V:
            self.V_TYPE_ = LINELEM_ARRAY[4]
        elif self.ElemType==linelem.ElemType_I:
            self.V_TYPE_ = LINELEM_ARRAY[4]
            if not self.V_TYPE_==linelem.V_TYPE_DC:
                print("Error, V_TYPE_ of {} not supported".format(self.V_TYPE_))
                exit(1)
        elif self.ElemType==linelem.ElemType_C:
            self.Init_Voltage = LINELEM_ARRAY[4]
        elif self.ElemType==linelem.ElemType_L:
            self.Init_Current = LINELEM_ARRAY[4]
        else:
            print("Error, ElemType of {} not supported".format(self.ElemType))
            exit(1)
    def __repr__(self):
        return "ElemType: %s\nInit_Value: %s\nNode1: %s\nNode2: %s\n" % (self.ElemType, self.Init_Value, self.Node1, self.Node2)


def circuit_from_filename(filename):
    # Opening JSON file
    f = open(filename)
    data = json.load(f)
    f.close()
    return data

# Follow the general flow for nonlinear transient analysis discussed in class.
# Use the multi-dimensional Newton-Raphson method in the inner loop to do the nonlinear iterations.
# To get the correct solution at t=0, you will need to do a nonlinear DC analysis.

def mna_from_circuit(circuit):
    # https://lpsa.swarthmore.edu/Systems/Electrical/mna/MNA3.html
    linear_elements = None
    try:
        linear_elements = [linelem(circuit["LINELEM"][i]) for i in range(len(circuit["LINELEM"]))]
    except:
        linear_elements = [ linelem(circuit["LINELEM"]) ]
    n = None
    try:
        n = len(circuit["NODES"])
    except:
        n = 1
    m = 0
    for e in linear_elements:
        if e.ElemType==linelem.ElemType_V:
            m += 1
        elif e.ElemType==linelem.ElemType_C:
            m += 1
    G = np.zeros((n, n))
    B = np.zeros((n, m))
    I = np.zeros((n, 1))
    E = np.zeros((m, 1))
    E_i = 0
    for e in linear_elements:
        if e.ElemType==linelem.ElemType_R:
            if not e.Node1==-1:
                G[e.Node1-1, e.Node1-1] += 1.0 / e.Init_Value
            if not e.Node2==-1:
                G[e.Node2-1, e.Node2-1] += 1.0 / e.Init_Value
            if not (e.Node1==-1 or e.Node2==-1):
                G[e.Node1-1, e.Node2-1] -= 1.0 / e.Init_Value
                G[e.Node2-1, e.Node1-1] -= 1.0 / e.Init_Value
        elif e.ElemType==linelem.ElemType_V:
            if not e.Node1==-1:
                B[e.Node1-1, E_i] = -1
            if not e.Node2==-1:
                B[e.Node2-1, E_i] = 1
            E[E_i, 0] = -e.Init_Value
            E_i += 1
        elif e.ElemType==linelem.ElemType_C:
            if not e.Node1==-1:
                B[e.Node1-1, E_i] = -1
            if not e.Node2==-1:
                B[e.Node2-1, E_i] = 1
            E[E_i, 0] = -e.Init_Voltage
            E_i += 1
        elif e.ElemType==linelem.ElemType_I:
            if not e.Node1==-1:
                I[e.Node1-1, 0] += e.Init_Value
        elif e.ElemType==linelem.ElemType_L:
            if not e.Node1==-1:
                I[e.Node1-1, 0] += e.Init_Current
        else:
            print("Error, ElemType of {} not supported".format(e.ElemType))
            exit(1)

    C = B.T
    D = np.zeros((m, m))
    GC = np.concatenate((G, C), axis=0)
    BD = np.concatenate((B, D), axis=0)
    A = np.concatenate((GC, BD), axis=1)

    z = np.concatenate((I, E))

    return [A, z]


if __name__=='__main__':
    circuit = circuit_from_filename('circuits_json/dc_r.json')
    A, z = mna_from_circuit(circuit)
    x = np.linalg.solve(A, z)
    print(A)
    print(x)
    print(z)
