import subprocess
import os

ESPRESO_TESTS = os.path.dirname(os.path.abspath(__file__))
ESPRESO_ROOT = os.path.dirname(ESPRESO_TESTS)
EXAMPLES = os.path.join(ESPRESO_TESTS, "examples")

if not os.path.isdir(EXAMPLES):
    os.mkdir(EXAMPLES)

class Mesh:

    def __init__(self, path, subdomains, fixpoints=0):
        self.path = path
        self.subdomains = subdomains
        self.fixpoints = fixpoints

    def get_path(self):
        return self.path

class Generator(Mesh):

    def __init__(self, path, subdomains, fixpoints):
        Mesh.__init__(self, path, subdomains, fixpoints)

class CubeGenerator(Generator):

    def __init__(self, etype, subdomains, elements):
        self.config = os.path.join(EXAMPLES, "cube.txt")
        Generator.__init__(self, self.config, subdomains[0] * subdomains[1] * subdomains[2], 8)
        self.etype = etype
        self.subdomains = subdomains
        self.elements = elements
        self.assemble()

    def assemble(self):
        def set(parameter, value):
            file.write("{0} = {1}\n".format(parameter, value))

        file = open(self.config, "w")

        set("SHAPE", 0)
        set("ELEMENT_TYPE", self.etype)
        set("SUBDOMAINS_X", self.subdomains[0])
        set("SUBDOMAINS_Y", self.subdomains[1])
        set("SUBDOMAINS_Z", self.subdomains[2])
        set("ELEMENTS_X", self.elements[0])
        set("ELEMENTS_Y", self.elements[1])
        set("ELEMENTS_Z", self.elements[2])
        set("DIRICHLET_BOTTOM_X", 0)
        set("DIRICHLET_BOTTOM_Y", 0)
        set("DIRICHLET_BOTTOM_Z", 0)

        file.close()

class Assembler:

    def __init__(self, discretization, assembler):
        self.discretization = discretization
        self.assembler = assembler

class Solver:

    def __init__(self, iterations, epsilon):
        self.iteration = iterations
        self.epsilon = epsilon

class Espreso:

    def __init__(self, mesh, assembler=None, solver=None):
        self.mesh = mesh
        self.assembler = assembler
        self.solver = solver

    def run(self, processes):
        program = [ "mpirun" ]
        program.append("-np={0}".format(processes))
        program.append(os.path.join(ESPRESO_ROOT, "espreso"))
        program.append(self.mesh.get_path())
        result = subprocess.Popen(program, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        output, error = result.communicate()
        print output
        print error








