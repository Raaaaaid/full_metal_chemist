
from collections import defaultdict
from copy import copy


class UnlockedMolecule(Exception):
    pass


class LockedMolecule(Exception):
    pass


class InvalidBond(Exception):
    pass


class EmptyMolecule(Exception):
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

    def mutate(self, elt):
        new_valence_number = self.valence_numbers[elt]
        try:
            assert len(self.bonded_atoms) <= new_valence_number
        except AssertionError:
            raise InvalidBond
        else:
            old_element, old_weight = self.element, self.weight
            self.element = elt
            self.weight = self.atomic_weights[elt]
            self.valence_number = new_valence_number
            return old_element, old_weight

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

    @staticmethod
    def delete_bond(atom1, atom2):
        try:
            assert atom1 != atom2
            assert atom1 in atom2.bonded_atoms
            assert atom2 in atom1.bonded_atoms
        except AssertionError:
            raise InvalidBond
        else:
            atom1.bonded_atoms.remove(atom2)
            atom2.bonded_atoms.remove(atom1)


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
            return self
        else:
            raise LockedMolecule

    def bounder(self, *args):
        if not self.locked:
            for atom_ind1, branch1, atom_ind2, branch2 in args:
                atom1 = self.branches[branch1][atom_ind1 - 1]
                atom2 = self.branches[branch2][atom_ind2 - 1]
                Atom.bond(atom1, atom2)
            return self
        else:
            raise LockedMolecule

    def mutate(self, *args):
        if not self.locked:
            for atom_ind, branch, elt in args:
                atom = self.branches[branch][atom_ind - 1]
                old_element, old_weight = atom.mutate(elt)
                self.formula_dict[old_element] -= 1
                self.formula_dict[elt] += 1
                self._molecular_weight -= old_weight
                self._molecular_weight += atom.weight
            return self
        else:
            raise LockedMolecule

    def add(self, *args):
        if not self.locked:
            for atom_ind, branch, elt in args:
                branch_atom = self.branches[branch][atom_ind - 1]
                new_atom = Atom(elt=elt, id_=len(self.atoms) + 1)
                Atom.bond(branch_atom, new_atom)
                self.formula_dict[elt] += 1
                self._molecular_weight += new_atom.weight
                self.atoms.append(new_atom)
            return self
        else:
            raise LockedMolecule

    def add_chaining(self, atom_ind, branch, elt, *elts):
        if not self.locked:
            bonded_atoms = []
            try:
                self.add((atom_ind, branch, elt))
                bonded_atoms.append((self.branches[branch][atom_ind - 1], self.atoms[-1]))
                for elt in elts:
                    new_atom = Atom(elt=elt, id_=len(self.atoms) + 1)
                    Atom.bond(self.atoms[-1], new_atom)
                    self.formula_dict[elt] += 1
                    self._molecular_weight += new_atom.weight
                    self.atoms.append(new_atom)
                    bonded_atoms.append((bonded_atoms[-1][1], new_atom))
                return self
            except InvalidBond:
                for atom1, atom2 in bonded_atoms[::-1]:
                    Atom.delete_bond(atom1, atom2)
                    self.atoms.remove(atom2)
                    self.formula_dict[atom2.element] -= 1
                    self._molecular_weight -= atom2.weight
                raise InvalidBond
        else:
            raise LockedMolecule

    def closer(self):
        if not self.locked:
            self.locked = True
            atoms_to_check = copy(self.atoms)
            for atom in atoms_to_check:
                while len(atom.bonded_atoms) < atom.valence_number:
                    hydrogen = Atom(elt='H', id_=len(self.atoms) + 1)
                    Atom.bond(atom, hydrogen)
                    self.formula_dict['H'] += 1
                    self._molecular_weight += hydrogen.weight
                    self.atoms.append(hydrogen)
            return self
        else:
            raise LockedMolecule

    def unlock(self):
        if self.locked:
            self.locked = False
            atoms_to_check = copy(self.atoms)
            for atom in atoms_to_check:
                if atom.element == 'H':
                    Atom.delete_bond(atom, atom.bonded_atoms[0])
                    self.atoms.remove(atom)
                    self.formula_dict['H'] -= 1
                    self._molecular_weight -= atom.weight
            for i, branch in self.branches.items():
                self.branches[i] = [atom for atom in branch if atom.element != 'H']
            for i, atom in enumerate(self.atoms):
                atom.id = i + 1
            self.branches = {i: branch for i, branch in self.branches.items() if branch}
            if not self.branches:
                raise EmptyMolecule
            self.branches = {i + 1: branch for i, branch in enumerate(self.branches.values())}
            return self
        else:
            raise UnlockedMolecule
