# All functions needed for post processing CIDMD simulations #


import re
import numpy as np
import pandas as pd
import os
import sys
import glob
import matplotlib
import matplotlib.pyplot as plt
from pathlib import Path
from collections import Counter
from datetime import date

today = date.today()


def get_molinfo(mol_info_infile):

    with open(mol_info_infile) as f:
        mol_name = f.readline().split(" ")[-1].strip()
        mol_nameid = f.readline().split(" ")[-1].strip()
        mol_formula = f.readline().split(" ")[-1].strip()
        mol_weight = f.readline().split(" ")[-1].strip()
        mol_exactmass = f.readline().split(" ")[-1].strip()

    cwd = os.getcwd()
    mol_id = cwd.split("/")[-3]
    mol_name2 = cwd.split("/")[-5]
    mol_protcase = cwd.split("/")[-4]
    return (
        mol_name,
        mol_nameid,
        mol_id,
        mol_formula,
        mol_weight,
        mol_exactmass,
        mol_protcase,
    )


def parse_molecules(infile, collect_mol):

    keep_molecules = []
    keep_molecules_found = []
    keep_molecules_known = []
    keep_molecules_transient = []
    string1 = " Found"
    string2 = "Known"
    string3 = "Transient"

    with open(infile, "r") as f:
        for line in f:
            if collect_mol == "all":
                if string1 in line or string2 in line or string3 in line:
                    keep_molecules.append(line.split()[1])
            elif collect_mol == "found":
                if string1 in line:
                    keep_molecules.append(line.split()[1])
            elif collect_mol == "known":
                if string2 in line:
                    keep_molecules.append(line.split()[1])
            elif collect_mol == "transient":
                if string3 in line:
                    keep_molecules.append(line.split()[1])
    return keep_molecules


def sorting(mols):

    sorted_mols = sorted(mols)
    sorted_count = {}
    for mol in sorted_mols:
        sorted_count[mol] = mols[mol]
    return sorted_count


def count_molecules(molecules):

    mols = Counter(molecules)
    return sorting(mols)


mass_table = {
    "H": 1.00783,
    "C": 12.00000,
    "N": 14.0031,
    "O": 15.9949,
    "P": 30.9738,
    "S": 32.9715,
    "Ar": 35.967546,
}


def load_frags(infilename):
    frags = []
    with open(infilename, "r") as f:
        for line in f:
            l = line.split()
            frags.append(l[-1])
    return frags


def load_count(infilename):
    count = []
    with open(infilename, "r") as f:
        for line in f:
            l = line.split()
            count.append(float(l[0]))
    # print('count for frags =')
    # print(count)
    return count


def get_mass_take2(my_frags_list):

    frags_map = []
    my_frags_mass_list = []
    for frag in my_frags_list:
        frag_map = re.findall(r"([A-Z]?)(\d*)", frag)
        frags_map.append(frag_map)
    for frag_map in frags_map:
        collect_atom_mass = []
        for i in frag_map:
            print(i)
            if len(i[-1]) != 0:
                atom_base_mass = mass_table.get(i[0]) * float(i[-1])
                collect_atom_mass.append(atom_base_mass)
            if len(i[-1]) == 0:
                atom_base_mass = mass_table.get(i[0])
                collect_atom_mass.append(atom_base_mass)
        real_collect_atom_mass = []
        for val in collect_atom_mass:
            if val != None:
                real_collect_atom_mass.append(val)
        tot_mass = 0.0
        for i in real_collect_atom_mass:
            tot_mass += float(i)
        my_frags_mass_list.append(tot_mass)
    return my_frags_mass_list


def write_cidmd_jdx(
    today,
    mol_name,
    mol_formula,
    mol_weight,
    mol_id,
    collect_rxntype,
    my_mass_list,
    my_frags_list,
    my_count_list,
):

    avgcharge_outfile = "../results/pop_" + collect_rxntype + ".out"
    countrxn_outfile = "../results/rxn_" + collect_rxntype + ".out"
    a = np.genfromtxt(avgcharge_outfile, usecols=0, dtype=str)
    b = np.genfromtxt(avgcharge_outfile, usecols=1, dtype=float)
    c = np.genfromtxt(countrxn_outfile, usecols=0, dtype=int)
    my_frags = []
    my_emzie = []
    my_abundance = []
    for i in range(len(b)):
        if b[i] > 0.5:
            my_frags.append(a[i])
            my_emzie.append(my_mass_list[i])
            my_abundance.append(my_count_list[i])
    cidmd_jdx_outfile = "../results/cidmd.jdx"

    with open(cidmd_jdx_outfile, "w") as cidmdout:
        for out in sys.stdout, cidmdout:
            print(f"##TITLE= CIDMD in-silico spectrum      ", file=out)
            print(f"##JCAMP-DX=Revision 4.10               ", file=out)
            print(f"##DATA TYPE=MASS SPECTRUM              ", file=out)
            print(f"##SAMPLE DESCRIPTION= {mol_id}         ", file=out)
            print(f"##NAMES= Jesi Lee                      ", file=out)
            print(f"##CAS NAME= {mol_name}                 ", file=out)
            print(f"##MOLFORM = {mol_formula}              ", file=out)
            print(f"##CIDMD spec created= {today}          ", file=out)
            print(f"##MP= -300                             ", file=out)
            print(f"##BP= -300                             ", file=out)
            print(f"##MW=  {mol_weight}                    ", file=out)
            print(f"##$RETENTION INDEX=0                   ", file=out)
            print(f"##$CONDENSED SPECTRUM=NO               ", file=out)
            print(f"##NPOINTS=  {len(my_frags_list)}       ", file=out)
            print(f"##XYDATA=(XY..XY) 1                    ", file=out)
            for i in range(len(my_emzie)):
                print(f"{my_emzie[i]:.4f} {my_abundance[i]:.2f} ", file=out)
    print(f"##%     Wrote {cidmd_jdx_outfile}     ##%")


def write_cidmd_msp(
    today,
    mol_name,
    mol_formula,
    mol_weight,
    mol_id,
    collect_rxntype,
    my_mass_list,
    my_frags_list,
    my_count_list,
):

    cidmd_msp_outfile = "../results/cidmd.msp"
    avgcharge_outfile = "../results/pop_" + collect_rxntype + ".out"
    countrxn_outfile = "../results/rxn_" + collect_rxntype + ".out"
    a = np.genfromtxt(avgcharge_outfile, usecols=0, dtype=str)
    b = np.genfromtxt(avgcharge_outfile, usecols=1, dtype=float)
    c = np.genfromtxt(countrxn_outfile, usecols=0, dtype=int)
    my_frags = []
    my_emzie = []
    my_abundance = []
    for i in range(len(b)):
        if b[i] > 0.5:
            my_frags.append(a[i])
            my_emzie.append(my_mass_list[i])
            my_abundance.append(my_count_list[i])
    cidmd_jdx_outfile = "../results/cidmd.jdx"
    with open(cidmd_msp_outfile, "w") as cidmdout:
        for out in sys.stdout, cidmdout:
            print(f"Name: {mol_name}                                 ", file=out)
            print(f"Notes: CIDMD_mol_id {mol_id}                     ", file=out)
            print(f"Precursor_type: [M+H]+                           ", file=out)
            print(f"Spectrum_type: MS2                               ", file=out)
            # print(f'PrecursorMZ: {mol_weight+1.}                    ', file=out)
            print(f"Instrument_type: in-silico CID                   ", file=out)
            print(f"Instrument: TeraChem CIDMD                       ", file=out)
            print(f"Ionization: manual                               ", file=out)
            print(f"Collision_gas: Ar                                ", file=out)
            print(f"Ion_mode: Positive                               ", file=out)
            print(f"Formula: {mol_formula}                           ", file=out)
            print(f"MW: {mol_weight}                                 ", file=out)
            print(f"Comments: CIDMD mass spectrum                    ", file=out)
            print(f"Num Peaks: {len(my_emzie)}                  ", file=out)

            for i in range(len(my_emzie)):
                print(f"{my_emzie[i]:.1f} {my_abundance[i]:.1f} ", file=out)
    print(f"##%     Wrote {cidmd_msp_outfile}     ##%")


def get_charge(popfilename):

    tot_frame = 0
    frag_charge = []
    frag_name = []

    popfile = open(popfilename, "r")
    for line in popfile:
        if "frame" in line:
            tot_frame += 1
            charge = line.split()[6]
            frag_charge.append(float(charge))

            frag = line.split()[0]
            frag_name.append(frag)
            # print(frag, charge)
    return frag_name[-1], frag_charge[-1]


def read_jdx_csv_pd(infilename, normalize_to=0):

    df = pd.read_csv(
        infilename,
        comment="#",
        header=None,
        names=("mz", "height"),
        dtype=({"mz": np.float64, "height": np.float64}),
        skipinitialspace=True,
        delim_whitespace=True,
    )
    if normalize_to != 0:
        df.height /= df.height.max()
        df.height *= normalize_to
    print("df shape is", df.shape)
    return df


def plot_cidmdjdx(qcjdx, mol_id, saveto=None):
    """
    plots cidmd given a jdx file.
    """

    fig = plt.figure(figsize=(20, 14), dpi=350)
    matplotlib.rc("font", **{"size": 24})
    ax = fig.add_subplot(111)

    qc_spec = ax.bar(
        qcjdx.mz,
        qcjdx.height,
        width=0.6,
        label="Theoretical spectrum",
        color="magenta",
    )

    ax.spines["right"].set_color("none")
    ax.spines["bottom"].set_position("zero")
    ax.spines["top"].set_color("none")
    ax.xaxis.set_ticks_position("bottom")
    ax.yaxis.set_ticks_position("left")
    ax.set_ylabel("Intensity")
    ax.set_xlabel("m/z")

    ax.legend(fontsize=20, loc="best")
    ax.set_title("Molecule " + mol_id, fontsize=16)
    ax.tick_params(axis="y", labelsize=18)
    ax.tick_params(axis="x", labelsize=18)

    label_peaks(ax, qc_spec)

    ticks = ax.get_yticks()
    ax.set_yticklabels([int(abs(tick)) for tick in ticks])
    ticks = ax.get_xticks()
    ax.set_xticklabels([int(abs(tick)) for tick in ticks])

    trim_xaxis = False
    if trim_xaxis:
        min_mz = 0.0
        # defined as the mz with the maximum intensity
        max_mz = qcjdx.mz.max()
        ax.set_xlim(min_mz, 1.1 * max_mz)
    fig.tight_layout()
    if saveto:
        fig.savefig(saveto)


def peaks(x, y, separation=0.5, ratio=1.01, threshold=0.0):
    """
    separation means the peak should be separated in mz by +- this number
    ratio means the peak should be this many times greater than the other points
    threshhold means the height should be at least this number to be a peak
    """
    assert separation > 0 and ratio > 0
    x = np.asarray(x)
    y = np.asarray(y)

    # sort them to make sure we always prefer the largest peaks first
    srt = np.abs(y).argsort()[::-1]

    xp = []
    yp = []

    for i in range(len(x)):
        i = srt[i]
        mask = np.abs(x - x[i])
        mask = mask < separation
        if not mask.any():
            print("skipping because empty mask")
            continue
        group = y[mask]

        # quick way to just get the 2 largest values
        if len(group) > 1:
            max_indices = np.argpartition(np.abs(group), len(group) - 2)[-2:]
            xval = x[mask][max_indices[1]]
            yval = group[max_indices[1]]
            second_max = group[max_indices[0]]
        else:
            idx = group.argmax()
            xval = x[mask][idx]
            yval = group[idx]
            second_max = 0

        print("Considering xval=", xval, "yval=", yval, "secondmax=", second_max)
        if second_max == 0:
            ratio_i = ratio * 2
        else:
            ratio_i = yval / second_max
        uniq = [abs(xpi - xval) > separation for xpi in xp]
        print("******ADD", xval, yval)

        xp.append(xval)
        yp.append(yval)

    # print(" DONE")
    print("xp = ", xp)
    print("yp = ", yp)
    return xp, yp


def label_peaks(ax, rects, peaks_only=True):
    """
    Attach a text label above each bar displaying its height
    """

    dat = [[rect.get_x() + rect.get_width() / 2.0, rect.get_height()] for rect in rects]
    if peaks_only:
        x = [v[0] for v in dat]
        y = [v[1] for v in dat]
        xp, yp = peaks(x, y)
    else:
        xp = [v[0] for v in dat if v[1] > 0.0]
        yp = [v[1] for v in dat if v[1] > 0.0]
    for mz, height in zip(xp, yp):
        ax.annotate(
            "%.2f" % mz,
            (mz, np.sign(height) * 4 + height),
            ha="center",
            va="center",
            fontsize=18,
        )


##% ##% ##% ##% ##%


def lets_parse_and_clean(collect_rxntype):
    base = os.path.join("..", "calcs")
    assert os.path.exists(base)
    molecules = []
    logs = glob.glob(os.path.join(base, "*", "gathered", "learn_rxn.log"))
    N = len(logs)
    for infile in logs:
        # molecules.extend(parse_molecules(infile))
        molecules.extend(parse_molecules(infile, collect_rxntype))
    counts = count_molecules(molecules)
    # del counts['Ar']

    for key in list(counts.keys()):
        if "Ar" in key:
            del counts[key]
    return base, counts


def write_sy_outfile(sy, collect_rxntype):

    cidmd_sy_outfile = "../results/sy_" + collect_rxntype + ".out"
    with open(cidmd_sy_outfile, "w") as sy_val:
        for out in sys.stdout, sy_val:
            print(f"{sy:.6f}", file=out)
    print(f"##%     Wrote {cidmd_sy_outfile}     ##%")
    return cidmd_sy_outfile


def lets_write_rxn_outfile(counts, collect_rxntype):
    rxn_outfile = "../results/rxn_" + collect_rxntype + ".out"
    with open(rxn_outfile, "w") as rxns:
        for out in sys.stdout, rxns:
            # only prints
            for mol, count in counts.items():
                print(f"{count:10d} {mol:s}", file=out)
            # print(f"Total: {N}", file=out)
    print()
    print(f"##%     Wrote {rxn_outfile}     ##%")
    print()
    return rxn_outfile


def lets_get_my_frags_count_mass_list(rxn_outfile):
    ##% get_mass
    # infilename = sys.argv[1]
    infilename_toGetMass = rxn_outfile
    my_frags_list = load_frags(infilename_toGetMass)
    my_count_list = load_count(infilename_toGetMass)
    # my_mass_list = mass_list(my_frags_list)
    # spec_out(infilename_toGetMass, my_frags_list, my_mass_list, my_count_list)
    return my_frags_list, my_count_list


def lets_get_spec(rxn_outfile, my_frags_list, my_count_list, my_mass_list):
    ##% get spec
    infilename_toGetSpec = rxn_outfile
    my_count_array = np.array(my_count_list, dtype=np.float16)
    my_mass_array = np.array(my_mass_list, dtype=np.int32)
    # jdx_output(infilename_toGetSpec, my_frags_list, my_mass_array, my_count_array)
    return my_count_array, my_mass_array


def lets_get_charges(base, collect_rxntype, my_frags_list, my_count_list, my_mass_list):
    ##% get charges
    collect_avg_charges = []
    charges = []
    poplogs = glob.glob(
        os.path.join(base, "*", "gathered", "molecules", "molecule_*.pop")
    )
    M = len(poplogs)

    collect_frags = []
    collect_charges = []

    for poplog in poplogs:
        frag, charge = get_charge(poplog)
        collect_frags.append(frag)
        collect_charges.append(charge)
    collect_avg_charges = []

    print()
    report_outfile = "../results/report_" + collect_rxntype + ".out"
    print()
    print(my_mass_list)
    # print(my_mass_array)
    with open(report_outfile, "w") as reportout:
        for out in sys.stdout, reportout:
            # print('my_frag      n_frags   mean        std        min        max', file = out )
            # print('my_frag       my_mass   n_frags   mean        std        min        max', file = out )
            print(
                "my_frag       my_mass    counts   mean        std        min        max",
                file=out,
            )

            for j in range(len(my_frags_list)):
                sum_collect_charges = 0.0
                n_frags = 0
                collect_charges_for_frags = []
                print("j       = ", j)
                for i in range(len(collect_frags)):
                    if collect_frags[i] == my_frags_list[j]:
                        collect_charges_for_frags.append(collect_charges[i])
                        print("i = ", i)
                        print(collect_frags[i], collect_charges[i])
                        n_frags += 1
                        sum_collect_charges += collect_charges[i]
                x = np.array(collect_charges_for_frags)
                avg_charge = sum_collect_charges / n_frags
                collect_avg_charges.append(avg_charge)
                mean = x.mean()
                std = x.std()
                minimum = x.min()
                maximum = x.max()
                print(
                    f"{my_frags_list[j]:10s} {my_mass_list[j]:10.4f} {int(my_count_list[j]):6d} {mean:10.4f} {std:10.4f} {minimum:10.4f} {maximum:10.4f}",
                    file=out,
                )
    print()
    print(f"##%     Wrote {report_outfile}     ##%")
    print()

    avgcharge_outfile = "../results/pop_" + collect_rxntype + ".out"
    with open(avgcharge_outfile, "w") as popout:
        for out in sys.stdout, popout:
            for i in range(len(my_frags_list)):
                print(
                    f"{my_frags_list[i]:12s} {collect_avg_charges[i]:6f}",
                    file=out,
                )
    print()
    print(f"##%     Wrote {avgcharge_outfile}     ##%")
    print()
    return collect_avg_charges


def find_mol_ion_info(my_frags_list, my_count_list, my_mass_list, mol_weight):

    mol_ion_mass = int(mol_weight) + 1.0
    for i in range(len(my_mass_list)):
        # print(i, my_mass_list[i])
        diff = float(my_mass_list[i]) - float(mol_ion_mass)
        mol_ion_exmass = my_mass_list[i]
        mol_ion_form = my_frags_list[i]
        mol_ion_count = my_count_list[i]
        mol_ion_index = i
        if abs(diff) < 0.7:
            mol_ion_exmass = my_mass_list[i]
            mol_ion_form = my_frags_list[i]
            mol_ion_count = my_count_list[i]
            mol_ion_index = i
    return mol_ion_exmass, mol_ion_form, mol_ion_count, mol_ion_index


def calc_survivalyield(mol_ion_count, my_count_list):
    sy = 0.0
    I_p = 0.0
    I_f = 0.0

    I_p = mol_ion_count
    I_f = sum(my_count_list)
    sy = I_p / I_f
    return sy


def get_gputime(infile_gputimelog):

    infile = "../results/gputime.log"
    D = np.loadtxt(infile_gputimelog)
    totaltime_sec = np.sum(D)
    totaltime_hr = totaltime_sec / 3600.0

    outfilename = "../results/gputime.txt"
    with open(outfilename, "w") as cidmdtimeout:
        for out in sys.stdout, cidmdtimeout:
            print(
                f"total time in seconds = {totaltime_sec}                                 ",
                file=out,
            )
            print(
                f"total time in hours = {totaltime_hr}                                 ",
                file=out,
            )
    print(f"##%     Wrote {outfilename}     ##%")


def plot_ar_vel(infile):
    """
    This file calculates and plots terachem raw velocity or converted KE data.
    make sure to check do_conv
    """
    delta_i = []
    min_i = []
    Ar_mass = 0.039948  # in Kg/mol
    conversion = 10 ** (-10)  # A to m

    to_Jmol = 0.5 * Ar_mass * conversion**2 / (4.888821e-14) ** 2

    with open("plot_ar_vel.in", "r") as fin:
        for line in fin:
            print(line)
            data = np.loadtxt(
                line.replace("\n", "")
            )  # shape (1001,2)=(samples, dimensions) where dim1=timestep dim2=vel
            data = data.T  # (2,1001)

            if do_conv:
                data[1] *= data[1]
                data[1] *= to_Jmol
                data[1] *= 0.000239  # to kcal/mol
            plt.plot(*data)

    outfile = ""
    if do_conv:
        outfile = "ar_vel_KE.png"
    else:
        outfile = "ar_vel_raw.png"

    plt.xlabel("time step")
    plt.ylabel("KE")
    plt.title(label=outfile.split(".")[0])
    plt.savefig("../results/" + outfile)


def main():

    print("Today's date:", today)
    collect_rxntype = (
        "found"  # choose from the options: found, known, transient, and all #
    )

    mol_info_infile = "../../../../mol_info.in"
    (
        mol_name,
        mol_name,
        mol_id,
        mol_nameid,
        mol_formula,
        mol_weight,
        mol_protcase,
    ) = get_molinfo(mol_info_infile)
    base, counts = lets_parse_and_clean(collect_rxntype)

    rxn_outfile = lets_write_rxn_outfile(counts, collect_rxntype)

    my_frags_list, my_count_list = lets_get_my_frags_count_mass_list(rxn_outfile)
    my_mass_list = get_mass_take2(my_frags_list)
    my_count_array, my_mass_array = lets_get_spec(
        rxn_outfile, my_frags_list, my_count_list, my_mass_list
    )

    collect_avg_charges = lets_get_charges(
        base, collect_rxntype, my_frags_list, my_count_list, my_mass_list
    )

    write_cidmd_jdx(
        today,
        mol_name,
        mol_formula,
        mol_weight,
        mol_id,
        collect_rxntype,
        my_mass_list,
        my_frags_list,
        my_count_list,
    )

    write_cidmd_msp(
        today,
        mol_name,
        mol_formula,
        mol_weight,
        mol_id,
        collect_rxntype,
        my_mass_list,
        my_frags_list,
        my_count_list,
    )

    (
        mol_ion_exmass,
        mol_ion_form,
        mol_ion_count,
        mol_ion_index,
    ) = find_mol_ion_info(my_frags_list, my_count_list, my_mass_list, mol_weight)
    # sy = calc_survivalyield(mol_ion_count, my_count_list)
    # cidmd_sy_outfile = write_sy_outfile(sy, collect_rxntype)

    cidmd_jdx_infile = "../results/cidmd.jdx"

    jdx_pd = read_jdx_csv_pd(cidmd_jdx_infile, normalize_to=0)
    fig = plot_cidmdjdx(jdx_pd, mol_id, saveto="../results/cidmd." + mol_id + ".png")

    jdx_pd = read_jdx_csv_pd(cidmd_jdx_infile, normalize_to=1)
    fig = plot_cidmdjdx(jdx_pd, mol_id, saveto="../results/cidmd." + mol_id + "norm.png")

    infile_gputimelog = "../results/gputime.log"
    get_gputime(infile_gputimelog)

    infile = "plot_ar_vel.in"
    plot_ar_vel(infile)

    ##%##% jesiplot_cidmd_spec_slide.py is not included ##%
    print(
        " CIDMD post-processing done. Please check ../results/ and run ../compare_"
        + mol_name
    )


if __name__ == "__main__":
    main()
