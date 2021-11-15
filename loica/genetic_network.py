import sbol3
import networkx as nx
from .geneproduct import Regulator, Reporter

class GeneticNetwork():
    def __init__(self, vector=None):
        self.operators = []
        self.regulators = []
        self.reporters = []
        self.vector = vector

    def initialize(self):
        for regulator in self.regulators:
            regulator.initialize()
        for reporter in self.reporters:
            reporter.initialize()

    def add_operator(self, op):
        self.operators.append(op)

    def add_regulator(self, reg):
        self.regulators.append(reg)

    def add_reporter(self, rep):
        self.reporters.append(rep)

    def add_operators(self, ops):
        for op in ops:
            self.operators.append(op)

    def add_regulators(self, regs):
        for reg in regs:
            self.regulators.append(reg)

    def add_reporters(self, reps):
        for rep in reps:
            self.reporters.append(rep)

    def step(self, growth_rate=1, t=0, dt=0.1):
        for op in self.operators:
            expression_rate = op.expression_rate(t, dt)
            op.output.express(expression_rate)

        for regulator in self.regulators:
            regulator.step(growth_rate, dt)

        for reporter in self.reporters:
            reporter.step(growth_rate, dt)

    def to_graph(self):
        g = nx.DiGraph()
        for op in self.operators:
            #if type(op.output)==Regulator:
            if type(op.input)==list:
                for i in op.input:
                    g.add_edge(i, op)
            else:
                g.add_edge(op.input, op)
            g.add_edge(op, op.output)
        return g

    def draw(
        self,
        node_shape='o',
        node_size=200,
        linewidths=2,
        alpha=0.75,
        arrowsize=5,
        font_size=6,
        font_family='Tahoma',
        font_weight='bold',
        pos=nx.kamada_kawai_layout
        ):
        g = self.to_graph()
        pos = pos(g)
        nx.draw_networkx_nodes(
            g, 
            pos=pos, 
            node_color=[n.color for n in g.nodes], 
            node_shape=node_shape, 
            node_size=node_size,
            linewidths=linewidths,
            alpha=alpha
            )
        nx.draw_networkx_edges(
            g, 
            pos=pos, 
            width=1, 
            node_shape=node_shape, 
            node_size=node_size,
            arrowsize=arrowsize
            )
        nx.draw_networkx_labels(
            g,
            pos=pos,
            font_size=font_size,
            font_family=font_family,
            font_weight=font_weight
            )

    def to_sbol(self):
        # prototype for not gate interaction
        doc = sbol3.Document()
        geneticnetwork = sbol3.Component('geneticnetwork', sbol3.SBO_DNA)
        geneticnetwork.roles.append(sbol3.SO_ENGINEERED_REGION)
        loica_set = set()
        for op in self.operators:
            operator_comp = op.sbol_comp
            output_comp = op.output.sbol_comp
            input_comp = op.input.sbol_comp
            print(f'Input: {input_comp}, Operator {operator_comp}, Output:{output_comp}')
            operator_sc = sbol3.SubComponent(operator_comp)
            output_sc = sbol3.SubComponent(output_comp)
            input_sc = sbol3.SubComponent(input_comp)
            #if not in doc
            loica_set.add(operator_comp)
            loica_set.add(output_comp)
            loica_set.add(input_comp)
            # TU Component
            tu = sbol3.Component(f'TU_{op.input.name}_{op}_{op.output.name}', sbol3.SBO_DNA) #generalize to multi input/output TUs
            tu.roles.append(sbol3.SO_ENGINEERED_REGION)
            tu.features = [operator_sc, output_sc, input_sc]                  
            tu.constraints = [sbol3.Constraint(sbol3.SBOL_PRECEDES, operator_sc, output_sc)]
            # generate a sequence for the TU assuming assembly by type IIS REsnf both parts will have the fusion sites.
            # compare last 4 bp with thefirst 4 bp of the next part, given the preceds constraint.
            # if they are the same then delete one of them and concatenate the rest.
            # else error or comment TU sequence can not be generated, provide ways to add it.

            # Protein Component
            protein = sbol3.Component(f'protein_{op.input.name}_{op}_{op.output.name}', sbol3.SBO_PROTEIN)
            protein.roles.append(sbol3.SO_TRANSCRIPTION_FACTOR)
            protein_sc = sbol3.SubComponent(protein)
            tu.features.append(protein_sc)
            # obtain TU subcomponents sequences, specially CDS and flanking parts sequences
            # look for ATG on the CDS and upstream part sequences (in the case of MoClo the ATG is in the fusion sites)
            # look for stop codons on frame with the ATG.
            # add translated the sequence between the ATG and the stop codon as protein sequence.
            #protein.sequence = tu.cds.sequence

            # Production Interaction (change TU for CDS)
            output_participation = sbol3.Participation(roles=[sbol3.SBO_TEMPLATE], participant=output_sc)
            protein_participation = sbol3.Participation(roles=[sbol3.SBO_PRODUCT], participant=protein_sc)
            production = sbol3.Interaction(types=[sbol3.SBO_GENETIC_PRODUCTION], participations=[output_participation, protein_participation])
            tu.interactions.append(production)
            # Inhibition Interaction
            input_participation = sbol3.Participation(roles=[sbol3.SBO_INHIBITOR], participant=input_sc)
            op_participation = sbol3.Participation(roles=[sbol3.SBO_INHIBITED], participant=operator_sc)
            interaction = sbol3.Interaction(types=[sbol3.SBO_INHIBITION], participations=[input_participation, op_participation])
            tu.interactions.append(interaction)
            # Model
            #op_model_attachment = sbol3.Attachment()
            op_model = sbol3.Model(f'LOICA_{op.input.name}_{op}_{op.output.name}_model', 
                            source='https://github.com/SynBioUC/LOICA/blob/master/loica/operators/not_.py',
                            language='http://identifiers.org/EDAM:format_3996', # http://edamontology.org/format_3996
                            framework='http://identifiers.org/SBO:0000062' #,attachments=[op]
                            )
            doc.add(op_model)
            tu.models.append(op_model)
            #tu.derived_from = []
            doc.add(tu)
            tu_sc = sbol3.SubComponent(tu)#(instance_of=tu.identity )
            geneticnetwork.features.append(tu_sc)
        loica_list = list(loica_set)
        #doc.add(loica_list) # waiting release
        for obj in loica_list:
            doc.add(obj)
        if len(geneticnetwork.features) > 1:
            for i in range(len(geneticnetwork.features)-1):
                geneticnetwork.constraints = [sbol3.Constraint(sbol3.SBOL_PRECEDES, geneticnetwork.features[i], geneticnetwork.features[i+1])]
        else: pass
        doc.add(geneticnetwork)
        return doc

        
