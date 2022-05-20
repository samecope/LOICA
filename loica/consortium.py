from .sample import Sample


class Consortium(Sample):
    """
    Representation of microbial consortium (a child of Sample)
    Consortium contains GeneticNetworks (which equal strains) and their respective 
    Metabolism. The strains are connected through extracellular space/medium.
    
    Attributes
    ----------
    strains : List[GeneticNetwork]
        strains that form consortium
    metabolism : List[Metabolism]
        metabolism that drives the genetic network of each strain
    media : str
        Name of the media in the consortium
    strain : List[str]
        Name of the strains in the consortium
    
    """

    """ 
    I need to have something similar to metabolism here
    So there will be either nested list or dictionary containing wires and their 
    extracelllar concentration

    For example, extracellular_conc = [[wire1, conc1], [wire2, conc2]]
    And each extracellular concentration will be set to 0 in the beginning
    extracellular_conc = [[wire1, 0], [wire2, 0]]

    Since supplement can be a wire, I need to change wire concentration to be equal 
    to supplement

    so let's say supplement = lc.Supplement(name="wire1", concentration=5)
    and wire1.name = "wire1"
    since supplement.name == wire1.name
    extracellular_conc = [[wire1, 5], [wire2, 0] at time 0
    """
    def __init__(self,
            strains=None, 
            metabolism=None, 
            media=None,
            strain=None,
            ):
        super().__init__(metabolism, media)
        self.strains = strains
        self.strain = strain
        #self.vector = self.genetic_network.vector

        if self.strains:
            # code that pulls reporters from all strains, and removes repeating
            # add code for latter (not sure if I need to remove repeating)
            self.reporters = []
            for cons_strain in self.strains:
                self.reporters.append(cons_strain.reporters)
        if metabolism:
            self.biomass = []
            for strain_metabolism in self.metabolism:
                self.biomass.append(strain_metabolism.biomass)
        self.supplements = {} 


    # function that checks whether supplement is the same as wire
    # if yes, set the extracellular concentration of the wire 
    def supplement_is_wire(self, supplement, wire):
        if supplement.name == wire.name:
            extracellular_conc[wire]=supplement.concentration
