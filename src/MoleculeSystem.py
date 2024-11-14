import abc


class MoleculeSystem(abc.ABC):

    @abc.abstractmethod
    def __init__(self, top_path: str, trj_path=None):
        ...

    @abc.abstractmethod
    def select_atoms(self, selection: str) -> None:
        ...
