
from collections import defaultdict


class LockedModule(Exception):
    pass


class Atom(object):

    atomic_weight = {
        'H': 1.0,
        'B': 10.8,
        'C': 12.0,
        'N': 14.0,
        'O': 16.0,
        'F': 19.0,
        'Mg': 24.3,
        'P': 31.0,
        'S': 32.1,
        'Cl': 35.5,
        'Br': 80.0,
    }

    def __init__(self, elt, id_):
        self.element = elt
        self.id = id_
        self.weight = self.atomic_weight[elt]

    def __hash__(self):
        return self.id

    def __eq__(self, other):
        return self.id == other.id


class Molecule(object):

    def __init__(self, name=''):
        self.name = name
        self.formula_dict = defaultdict(lambda: 0)
        self._molecular_weight = 0.0
        self.atoms = []
        self.branches = {}
        self.locked = False

    @property
    def formula(self):
        # C, H, O, then other elements in alphabetic order
        return ''

    @property
    def molecular_weight(self):
        return 0.0

    def brancher(self, *args):
        for num_of_carbons in args:
            self.branches[len(self.branches) + 1] = []
            for _ in range(num_of_carbons):
                carbon = Atom(elt='C', id_=len(self.atoms) + 1)
                self.atoms.append(carbon)
                self.branches[len(self.branches)].append(carbon)
                self._molecular_weight += carbon.weight
                self.formula_dict['C'] += 1
