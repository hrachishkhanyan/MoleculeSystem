import MDAnalysis as mda


def strip_water(top: str, trj=None, file_name="no_water", save_path=".", water_resname="TIP3") -> None:
    """
    Strips a residue type from the trajectory and saves it. The default residue type is TIP3 water.

    :param top: Path to topology
    :param trj: Path to trajectory
    :param file_name: Name of the new trajectory/topology file
    :param save_path: Path to the save directory
    :param water_resname: Residue name to strip from the trajectory
    :return: None
    """

    u = mda.Universe(top, trj)
    ag = u.select_atoms(f"not resname {water_resname}")

    print("Writing topology...")
    with mda.Writer(f"{save_path}/{file_name}.gro", ag.n_atoms) as W:
        W.write(ag)
    print("Writing trajectory...")

    if trj:
        with mda.Writer(f"{save_path}/{file_name}.xtc", ag.n_atoms) as W:
            for _ in u.trajectory:
                W.write(ag)

    print("Finished!")


