import abc


class MoleculeSystem(abc.ABC):

    @abc.abstractmethod
    def __init__(self):
        ...

    @abc.abstractmethod
    def select_atoms(self, selection: str) -> None:
        ...
