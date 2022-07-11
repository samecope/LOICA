from .impactor import *

class Degrader(Impactor):
    """
    A class that represents an enzyme that degrades a substrate (or substrates).
    ...
    
    Attributes
    ----------
    enzyme : Regulator | GeneProduct
        The enzyme that degrades substrate
    substrate : Regulator | Reporter | GeneProduct | List [Regulator | Reporter | GeneProduct]
        The substrate that is degraded by enzyme
    k1 : int | float | List [int | float]
        Binding constant to produce enzyme-substrate complex
    k1r : int | float | List [int | float]
        Dissociation constant to reduce the enzyme-substrate complex
    k2 : int | float | List [int | float]
        Degradation constant to reduce enzyme-substrate complex into enzyme only
    name : str, optional
        Name of the operator displayed on the network representation
    color: str, optional
        Color displayed on the network representation

    """

    def __init__(self, enzyme, substrate, k1, k1r, k2, name=None, uri=None, sbol_comp=None, color='skyblue'):
        super().__init__(enzyme, name, color)
        self.enzyme = enzyme

        if type(k1)==list:
            self.k1 = k1
        else:
            self.k1 = [k1]

        if type(k1r)==list:
            self.k1r = k1r
        else:
            self.k1r = [k1r]

        if type(k2)==list:
            self.k2 = k2
        else:
            self.k2 = [k2]

        if type(substrate)==list:
            self.substrate = substrate
        else:
            self.substrate = [substrate]

        self.es_complex = []
        for o in self.substrate:
            self.es_complex.append(0)

    def __str__(self):
        if self.name == None:
            return 'DEGRADER'
        else: return self.name
        
    def degradation_rate(self, dt):
        degradation_rate = []
        enzyme = []
        
        for i, substrate in enumerate(self.substrate):
            # test
            print(f'{(self.enzyme.concentration / len(self.substrate) - self.es_complex[i])} split enzyme conc')
            # if enzyme has multiple substrates, enzyme is split equally between all substrates (account for enzymes in complex)
            enzyme.append(self.enzyme.concentration / len(self.substrate) - self.es_complex[i])

            # calculate substrate degradation rate
            substrate_change_rate = -self.k1[i] * enzyme[i] * substrate.concentration + self.k1r[i] * self.es_complex[i] - self.k2[i] * self.es_complex[i]
            # test
            print(f'degr_rate is {substrate_change_rate}')
            degradation_rate.append(-substrate_change_rate)

            # update enzyme-substrate complex concentrations
            es_change = self.k1[i] * enzyme[i] * substrate.concentration - self.k1r[i] * self.es_complex[i] - self.k2[i] * self.es_complex[i]
            new_es = self.es_complex[i] + es_change * dt
            self.es_complex[i] = new_es

        return degradation_rate