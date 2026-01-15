import os
import polars as pl
import sys
from itertools import product

subs = pl.read_ipc(sys.argv[1]+'subs.ipc')


### from peptide and charge make m/z with M/C mods
# Monoisotopic *residue* masses (no H2O)
AA_MONO = {
    "A": 71.03711,
    "R": 156.10111,
    "N": 114.04293,
    "D": 115.02694,
    "C": 103.00919,
    "E": 129.04259,
    "Q": 128.05858,
    "G": 57.02146,
    "H": 137.05891,
    "I": 113.08406,
    "L": 113.08406,
    "K": 128.09496,
    "M": 131.04049,
    "F": 147.06841,
    "P": 97.05276,
    "S": 87.03203,
    "T": 101.04768,
    "W": 186.07931,
    "Y": 163.06333,
    "V": 99.06841,
}

WATER_MASS  = 18.01056
OX_MASS     = 15.994915   # Met oxidation
CAM_MASS    = 57.021464   # Carbamidomethyl on Cys
PROTON_MASS = 1.007276466812


def peptide_mz_all_states(seq: str, charge: int) -> dict[str, float]:
    """
    Calculate monoisotopic m/z for all combinations of:
      - Met unoxidized / oxidized (+15.994915)
      - Cys unmodified / CAM (+57.021464)

    Returns:
        dict mapping "M<pos>_...C<pos>_..." -> m/z, e.g.
          "M0C0"      : m/z with no Ox / no CAM
          "M3C1"      : m/z with Ox on M3, CAM on C1
          "M3_7C1_5"  : m/z with Ox on M3 & M7, CAM on C1 & C5
    """
    if charge <= 0:
        raise ValueError("Charge must be a positive integer")

    # Base neutral peptide mass (no Ox / no CAM)
    try:
        base_neutral = WATER_MASS + sum(AA_MONO[aa] for aa in seq)
    except KeyError as e:
        raise ValueError(f"Unknown amino acid in sequence: {e.args[0]}")

    # 0-based indices of M and C in the sequence
    m_positions = [i for i, aa in enumerate(seq) if aa == "M"]
    c_positions = [i for i, aa in enumerate(seq) if aa == "C"]

    nM = len(m_positions)
    nC = len(c_positions)

    result: dict[str, float] = {}

    # Each M can be 0 (unoxidized) or 1 (oxidized)
    # Each C can be 0 (unmodified) or 1 (CAM)
    for m_state in product([0, 1], repeat=nM):
        for c_state in product([0, 1], repeat=nC):
            neutral_mass = base_neutral
            neutral_mass += sum(s * OX_MASS for s in m_state)
            neutral_mass += sum(s * CAM_MASS for s in c_state)

            # 1-based positions of modified residues
            ox_m_pos = [m_positions[i] + 1 for i, s in enumerate(m_state) if s]
            cam_c_pos = [c_positions[i] + 1 for i, s in enumerate(c_state) if s]

            # Build key: M<positions> C<positions>, zeros if none
            if ox_m_pos:
                m_part = "M" + "_".join(str(p) for p in ox_m_pos)
            else:
                m_part = "M0"

            if cam_c_pos:
                c_part = "C" + "_".join(str(p) for p in cam_c_pos)
            else:
                c_part = "C0"

            key = f"{m_part}{c_part}"

            mz = (neutral_mass + charge * PROTON_MASS) / charge
            result[key] = mz

    # If there are no M and no C, ensure we still return M0C0
    if not m_positions and not c_positions:
        neutral_mass = base_neutral
        key = "M0C0"
        mz = (neutral_mass + charge * PROTON_MASS) / charge
        result[key] = mz

    return result


def canon_title(t: str, use_base: bool) -> str:
    """Normalize titles for matching."""
    t = t.strip()
    if use_base and ":" in t:
        # Strip trailing :something (e.g., :3 or :scan)
        return t.rsplit(":", 1)[0]
    return t

### make a fake mgf for the bases
def expand_mgf_with_mod_states_mz(
    df: pl.DataFrame,
    mgf_in_path: str,
    mgf_out_path: str,
    use_base_title: bool = False,
    debug_unmatched: bool = True,
):
    """
    Expand an MGF with mod-state variants for base peptides.

    Parameters
    ----------
    df : pl.DataFrame
        Must have columns:
          - 'Spectrum'  : string, must match (or canonically match) TITLE in MGF
          - 'base' : string, peptide sequence for that spectrum

    mgf_in_path : str
        Input MGF file (can be large; processed in streaming fashion).

    mgf_out_path : str
        Output MGF file; will contain one spectrum per (Spectrum, mod_state).

    use_base_title : bool
        If True, match on TITLE with trailing ':something' stripped in both
        DF and MGF. (e.g. 'kjkjljds_:xxx:4' -> 'kjkjljds_:xxx').

    debug_unmatched : bool
        If True, print how many DF titles were not found in the MGF and
        how many were skipped due to missing/invalid CHARGE.
    """
    if not {"Spectrum", "base"}.issubset(df.columns):
        raise ValueError("DataFrame must have 'Spectrum' and 'base' columns.")

    # --- 1. Build mapping from (canonical) Spectrum to peptide list ---
    title_to_peps: dict[str, list[str]] = {}
    for row in df.iter_rows(named=True):
        raw_title = str(row["Spectrum"])
        key = canon_title(raw_title, use_base=use_base_title)
        pep = str(row["base"])
        title_to_peps.setdefault(key, []).append(pep)

    wanted_titles = set(title_to_peps.keys())

    # Track what we actually matched and what was skipped due to charge
    matched_titles: set[str] = set()
    charge_missing_for: set[str] = set()

    # --- 2. Stream through MGF ---
    with open(mgf_in_path, "r") as fin, open(mgf_out_path, "w") as fout:
        in_block = False
        block_lines: list[str] = []
        current_title_raw: str | None = None
        current_title_key: str | None = None
        pepmass_idx: int | None = None
        title_idx: int | None = None
        charge: int | None = None

        for line in fin:
            stripped = line.strip()

            if not in_block:
                if stripped == "BEGIN IONS":
                    in_block = True
                    block_lines = [line]
                    current_title_raw = None
                    current_title_key = None
                    pepmass_idx = None
                    title_idx = None
                    charge = None
                continue

            # Inside BEGIN IONS ... END IONS
            block_lines.append(line)

            if stripped.startswith("TITLE="):
                current_title_raw = stripped[len("TITLE="):].strip()
                current_title_key = canon_title(current_title_raw, use_base=use_base_title)
                title_idx = len(block_lines) - 1

            elif stripped.startswith("PEPMASS="):
                pepmass_idx = len(block_lines) - 1

            elif stripped.startswith("CHARGE="):
                charge_str = stripped.split("=", 1)[1].strip()
                charge_str = charge_str.rstrip("+")
                try:
                    charge = int(charge_str)
                except ValueError:
                    charge = None

            elif stripped == "END IONS":
                # End of spectrum block: see if we want this one
                if current_title_key in wanted_titles:
                    if charge is None:
                        # We wanted it but couldn't parse CHARGE
                        charge_missing_for.add(current_title_key)
                    elif pepmass_idx is not None and title_idx is not None:
                        matched_titles.add(current_title_key)

                        for pep in title_to_peps[current_title_key]:
                            mod_mz_dict = peptide_mz_all_states(pep, charge)

                            # Base title for writing (optionally stripped)
                            base_title_for_write = (
                                canon_title(current_title_raw, use_base=True)
                                if use_base_title
                                else current_title_raw
                            )

                            for mod_str, theo_mz in mod_mz_dict.items():
                                new_block = list(block_lines)

                                # Patch TITLE
                                new_title_line = f"TITLE={base_title_for_write}_{mod_str}\n"
                                new_block[title_idx] = new_title_line

                                # Patch PEPMASS (overwrite mass, drop old intensity)
                                new_pepmass_line = f"PEPMASS={theo_mz:.6f}\n"
                                new_block[pepmass_idx] = new_pepmass_line

                                fout.writelines(new_block)
                    # else: malformed block (no TITLE or PEPMASS)

                # Reset for next block
                in_block = False
                block_lines = []
                current_title_raw = None
                current_title_key = None
                pepmass_idx = None
                title_idx = None
                charge = None

    # --- 3. Debugging info ---
    if debug_unmatched:
        unmatched = wanted_titles - matched_titles
        print(f"Total DF titles: {len(wanted_titles)}")
        print(f"Matched in MGF : {len(matched_titles)}")
        print(f"Missing charge : {len(charge_missing_for)}")
        print(f"Not found in MGF (by key): {len(unmatched)}")

        if unmatched:
            print("Example unmatched titles (up to 10):")
            for t in list(unmatched)[:10]:
                print("  ", t)


def mod_suffixes_for_seq(seq: str) -> list[str]:
    """
    Return all Cys/Met modification suffixes for a peptide sequence.

    Suffix format:
      - 'M0C0'       : no oxidized M, no CAM C
      - 'M3C1'       : Ox on M at pos 3, CAM on C at pos 1
      - 'M3_7C1_5'   : Ox on M3 & M7, CAM on C1 & C5
    """
    # 0-based indices of M and C
    m_positions = [i for i, aa in enumerate(seq) if aa == "M"]
    c_positions = [i for i, aa in enumerate(seq) if aa == "C"]

    nM = len(m_positions)
    nC = len(c_positions)

    # No M, no C → single state M0C0
    if nM == 0 and nC == 0:
        return ["M0C0"]

    suffixes: list[str] = []

    for m_state in product([0, 1], repeat=nM):
        for c_state in product([0, 1], repeat=nC):
            # 1‑based positions of modified residues
            ox_m_pos = [m_positions[i] + 1 for i, s in enumerate(m_state) if s]
            cam_c_pos = [c_positions[i] + 1 for i, s in enumerate(c_state) if s]

            if ox_m_pos:
                m_part = "M" + "_".join(str(p) for p in ox_m_pos)
            else:
                m_part = "M0"

            if cam_c_pos:
                c_part = "C" + "_".join(str(p) for p in cam_c_pos)
            else:
                c_part = "C0"

            suffixes.append(f"{m_part}{c_part}")

    return suffixes


def write_dep_spectrum_chunks(
    df: pl.DataFrame,
    n: int,
    out_prefix: str = "dep_spec_chunk",
    out_dir: str = ".",
):
    """
    From a Polars DataFrame with columns 'base' and 'Spectrum',
    create many text files in `out_dir` with tab-separated 'base' and 'Spectrum'
    values, expanded over all Met/Cys mod combinations.

    Each line:
        <base>\\t<Spectrum>_<suffix>

    Files are split into chunks of size `n` lines:
        <out_dir>/<out_prefix>_1.txt, <out_prefix>_2.txt, ...

    Parameters
    ----------
    df : pl.DataFrame
        Must contain columns 'base' and 'Spectrum'.
    n : int
        Maximum number of lines per output file.
    out_prefix : str
        Prefix for output file names.
    out_dir : str
        Directory where output files will be created (will be created if needed).
    """
    if not {"base", "Spectrum"}.issubset(df.columns):
        raise ValueError("DataFrame must have 'base' and 'Spectrum' columns.")

    if n <= 0:
        raise ValueError("n must be a positive integer.")

    # Ensure output directory exists
    os.makedirs(out_dir, exist_ok=True)

    file_index = 1
    line_count = 0
    out_f = None

    def open_new_file():
        nonlocal out_f, file_index, line_count
        if out_f is not None:
            out_f.close()
        fname = f"{out_prefix}_{file_index}.txt"
        out_path = os.path.join(out_dir, fname)
        out_f = open(out_path, "w")
        file_index += 1
        line_count = 0

    opened = False

    for row in df.iter_rows(named=True):
        pep = str(row["base"])
        spec = str(row["Spectrum"])

        suffixes = mod_suffixes_for_seq(pep)

        for suffix in suffixes:
            if not opened or line_count >= n:
                open_new_file()
                opened = True

            out_f.write(f"{pep}\t{spec}_{suffix}\n")
            line_count += 1

    if out_f is not None:
        out_f.close()


expand_mgf_with_mod_states_mz(subs, sys.argv[1]+'mgfall.mgf',sys.argv[1]+'mgfall_fakes.mgf')
write_dep_spectrum_chunks(subs, n=1, out_prefix="variants", out_dir=sys.argv[1]+'pepquery_input_base')