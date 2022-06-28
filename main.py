
from collections import defaultdict


class UnlockedMolecule(Exception):
    pass


class LockedMolecule(Exception):
    pass


class InvalidBond(Exception):
    pass


class Atom(object):

    atomic_weights = {
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

    valence_numbers = {
        'H': 1,
        'B': 3,
        'C': 4,
        'N': 3,
        'O': 2,
        'F': 1,
        'Mg': 2,
        'P': 3,
        'S': 2,
        'Cl': 1,
        'Br': 1,
    }

    def __init__(self, elt, id_):
        self.element = elt
        self.id = id_
        self.weight = self.atomic_weights[elt]
        self.valence_number = self.valence_numbers[elt]
        self.bonded_atoms = []

    def __hash__(self):
        return self.id

    def __eq__(self, other):
        return self.id == other.id

    def __str__(self):
        s = f'Atom({self.element}.{self.id}'
        if self.bonded_atoms:
            # TODO: implement __lt__ to sort self.bonded_atoms
            pass
        return s + ')'

    @staticmethod
    def bond(atom1, atom2):
        try:
            assert atom1 != atom2
            assert len(atom1.bonded_atoms) + 1 <= atom1.valence_number
            assert len(atom2.bonded_atoms) + 1 <= atom2.valence_number
        except AssertionError:
            raise InvalidBond
        else:
            atom1.bonded_atoms.append(atom2)
            atom2.bonded_atoms.append(atom1)


class Molecule(object):

    formula_order = ['C', 'H', 'O', 'B', 'N', 'F', 'Mg',
                     'P', 'S', 'Cl', 'Br']

    def __init__(self, name=''):
        self.name = name
        self.formula_dict = defaultdict(lambda: 0)
        self._molecular_weight = 0.0
        self.atoms = []
        self.branches = {}
        self.locked = False

    @property
    def formula(self):
        if self.locked:
            formula = ''
            for elt in self.formula_order:
                amount = self.formula_dict[elt]
                if amount == 1:
                    formula += elt
                elif amount > 1:
                    formula += elt + str(amount)
            return formula
        else:
            raise UnlockedMolecule

    @property
    def molecular_weight(self):
        if self.locked:
            return self._molecular_weight
        else:
            raise UnlockedMolecule

    def brancher(self, *args):
        if not self.locked:
            for num_of_carbons in args:
                self.branches[len(self.branches) + 1] = []
                for _ in range(num_of_carbons):
                    carbon = Atom(elt='C', id_=len(self.atoms) + 1)
                    self.atoms.append(carbon)
                    self.branches[len(self.branches)].append(carbon)
                    self._molecular_weight += carbon.weight
                    self.formula_dict['C'] += 1
        else:
            raise LockedMolecule

    def bounder(self, *args):
        if not self.locked:
            for carbon1, branch1, carbon2, branch2 in args:
                atom1 = self.branches[branch1][carbon1 - 1]
                atom2 = self.branches[branch2][carbon2 - 1]
                Atom.bond(atom1, atom2)
        else:
            raise LockedMolecule
