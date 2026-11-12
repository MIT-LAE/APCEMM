"""Migrate an APCEMM input file from v1.2 to v1.3.
You should still double check the output yourself.

Usage: python input-migration.py <old_input.yaml> <new_input.yaml>

The script reads the input file line by line, and applies four changes:

1. Delete every removed key. If the key opens a submenu, also delete its whole block.
   Also delete the comment lines directly above a removed key (no blank line between),
   if they are at least as indented as the key.
2. Rename the legacy time step keys in the TRANSPORT MENU to their current names.
   If the file already has the current name, delete the legacy key instead.
3. Set 'Turn on solid coagulation (T/F)' to F, because the YAML reader now rejects T.
4. Convert 'Ice growth timestep [min]' in the AEROSOL MENU to
   'Inner physics time step [s]' in the TRANSPORT MENU, to keep the v1.2 step size.

If a flag in the METEOROLOGICAL INPUT SUBMENU turns off the met file (entirely, or for
temperature, RH, wind shear, or vertical velocity), the script cannot convert the file
as v1.3 relies entirely on the netCDF file. Conversion fails.

Stale comments will still be in the new input file.
"""

import re
import sys
from pathlib import Path

# Copied from checkRemovedKey() in Code.v05-00/src/YamlInputReader/YamlInputReader.cpp.
REMOVED_SWEEP_KEYS = {
    "PARAM SWEEP SUBMENU",
    "Parameter sweep (T/F)",
    "Run Monte Carlo (T/F)",
    "Num Monte Carlo runs (int)",
}

REMOVED_UNUSED_KEYS = {
    # SIMULATION MENU: spectral solver leftovers and never-implemented run modes
    "Use threaded FFT (T/F)",
    "FFTW WISDOM SUBMENU",
    "Use FFTW WISDOM (T/F)",
    "Dir w/ write permission (string)",
    "SAVE FORWARD RESULTS SUBMENU",
    "Save forward results (T/F)",
    "ADJOINT OPTIMIZATION SUBMENU",
    "Turn on adjoint optim. (T/F)",
    "BOX MODEL SUBMENU",
    "Run box model (T/F)",
    "netCDF filename format (string)",
    # PARAMETER MENU: emission indices that reach Emission and stop
    "NOx [g(NO2)/kg_fuel] (double)",
    "CO [g/kg_fuel] (double)",
    "UHC [g/kg_fuel] (double)",
    # TRANSPORT MENU: spectral solver leftover, and updraft now comes from the met file
    "Fill Negative Values (T/F)",
    "PLUME UPDRAFT SUBMENU",
    "Turn on plume updraft (T/F)",
    "Updraft timescale [s] (double)",
    "Updraft veloc. [cm/s] (double)",
    # The two legacy time step keys are also in RENAMED_KEYS. The rename runs first.
    "Transport Timestep [min] (double)",
    "Transport and ice growth substep [s] (double)",
    # CHEMISTRY MENU: chemistry is no longer supported
    "CHEMISTRY MENU",
    "Turn on Chemistry (T/F)",
    "Chemistry Timestep [min] (double)",
    "Perform hetero. chem. (T/F)",
    "Photolysis rates folder (string)",
    # AEROSOL MENU: no coagulation runs, so its timestep never applied
    "Coag. timestep [min] (double)",
    "Turn on liquid coagulation (T/F)",
    # AEROSOL MENU: the outer and inner time steps of the TRANSPORT MENU set the ice growth step
    "Ice growth timestep [min] (double)",
    # METEOROLOGY MENU: met now comes from the netCDF file only
    "Use met. input (T/F)",
    "Init temp. from met. (T/F)",
    "Init RH from met. (T/F)",
    "Init wind shear from met. (T/F)",
    "Init vert. veloc. from met. data (T/F)",
    "Temp. Perturb. Timescale (min)",
    # DIAGNOSTIC MENU: species and P/L diagnostics needed chemistry
    "SPECIES TIMESERIES SUBMENU",
    "Save species timeseries (T/F)",
    "Species indices to include (list of ints)",
    "Aerosol indices to include (list of ints)",
    "PRODUCTION & LOSS SUBMENU",
    "Turn on P/L diag (T/F)",
    "Save O3 P/L (T/F)",
    # ADVANCED OPTIONS MENU: the documented scaling formula was never implemented
    "INITIAL CONTRAIL SIZE SUBMENU",
    "Base Contrail Depth [m] (double)",
    "Contrail Depth Scaling Factor [-] (double)",
    "Base Contrail Width [m] (double)",
    "Contrail Width Scaling Factor [-] (double)",
    # ADVANCED OPTIONS MENU: the EPM no longer models vortex sinking, so it needs no lapse rate
    "Ambient Lapse Rate [K/km] (double)",
}

REMOVED_DEPRECATED_KEYS = {
    # PARAMETER MENU -> METEOROLOGICAL PARAMETERS SUBMENU: met now comes from the netCDF file
    "Temperature [K] (double)",
    "R.Hum. wrt water [%] (double)",
    "Wind shear [1/s] (double)",
    # METEOROLOGY MENU: humidity, lapse rate, and diurnal options with no effect
    "HUMIDITY SCALING OPTIONS",
    "IMPOSE MOIST LAYER DEPTH SUBMENU",
    "IMPOSE LAPSE RATE SUBMENU",
    "Add diurnal variations (T/F)",
}

REMOVED_KEYS = REMOVED_SWEEP_KEYS | REMOVED_UNUSED_KEYS | REMOVED_DEPRECATED_KEYS

# Legacy key -> current key.
# These were the aliases in checkAliasedKey() and mergeYamlInputFiles().
RENAMED_KEYS = {
    ("TRANSPORT MENU",): {
        "Transport Timestep [min] (double)": "Outer time step [min] (double)",
        "Transport and ice growth substep [s] (double)": "Inner physics time step [s] (double)",
    },
}

# In v1.2 ice growth ran every 'Ice growth timestep [min]'. Now it runs every
# 'Inner physics time step [s]'. The two keys are in different menus and use different
# units, so convert_ice_growth_step() converts the value instead of a rename.
ICE_GROWTH_PATH = ("AEROSOL MENU",)
ICE_GROWTH_KEY = "Ice growth timestep [min] (double)"
ICE_GROWTH_SWITCH_KEY = "Turn on ice growth (T/F)"
TRANSPORT_PATH = ("TRANSPORT MENU",)
OUTER_STEP_KEY = "Outer time step [min] (double)"
INNER_STEP_KEY = "Inner physics time step [s] (double)"
# Accepted spellings of false, from parseBoolString() in YamlInputReader.hpp.
FALSE_STRINGS = {"F", "FALSE", "NO", "N", "0"}

# In v1.2 these flags could turn off the met file, entirely or for one variable.
# APCEMM then used the constant values and the moist layer and lapse rate options,
# which are all removed now. v1.3 reads met data only from the met file, so the script
# cannot convert such a file.
MET_INPUT_PATH = ("METEOROLOGY MENU", "METEOROLOGICAL INPUT SUBMENU")
MET_INPUT_FLAGS = (
    "Use met. input (T/F)",
    "Init temp. from met. (T/F)",
    "Init RH from met. (T/F)",
    "Init wind shear from met. (T/F)",
    "Init vert. veloc. from met. data (T/F)",
)

# Keys that we force to a specific setting.
FORCED_VALUES = {
    ("AEROSOL MENU",): {
        "Turn on solid coagulation (T/F)": "F",
    },
}

# In v1.2 the reader ignored this value and the code used 3.0e4 Pa.
# v1.3 correctly uses the value, though it does not impact contrail results.
TROPOPAUSE_KEY = "Tropopause Pressure [Pa] (double)"
TROPOPAUSE_OLD_HARDCODED_PA = 3.0e4

# A mapping line: indent, key (plain or quoted), colon, then the rest of the line.
# A plain YAML key cannot contain ": ", so the first ":" followed by a space or the
# end of the line is the separator.
KEY_LINE = re.compile(
    r"""^(?P<indent>[ ]*)
        (?P<key>"[^"]*"|'[^']*'|[^\s#][^#]*?)
        :(?P<rest>(?:\s.*)?)$""",
    re.VERBOSE,
)


class Line:
    """One line of the input file, with the YAML key it holds (if any)."""

    def __init__(self, text: str):
        self.text = text
        self.indent = len(text) - len(text.lstrip(" "))
        stripped = text.strip()
        self.is_blank = stripped == ""
        self.is_comment = stripped.startswith("#")
        self.key = None
        self.rest = None
        if not self.is_blank and not self.is_comment:
            match = KEY_LINE.match(text.rstrip("\r\n"))
            if match:
                self.key = match.group("key").strip()
                if self.key[0] in "\"'" and self.key[0] == self.key[-1]:
                    self.key = self.key[1:-1]
                self.rest = match.group("rest")
        self.path = ()  # Parent keys of this line, filled in by annotate_paths().
        self.delete = False


def annotate_paths(lines: list[Line]) -> None:
    """Give each key line the tuple of its parent keys."""
    stack: list[tuple[int, str]] = []
    for line in lines:
        if line.key is None:
            continue
        while stack and stack[-1][0] >= line.indent:
            stack.pop()
        line.path = tuple(key for _, key in stack)
        stack.append((line.indent, line.key))


def block_end(lines: list[Line], start: int) -> int:
    """Return the index after the last line of the block that starts at lines[start].

    The block holds all following lines that are more indented than lines[start].
    Blank lines and comments inside the block belong to it. Trailing blank lines and
    comments at the indent of the key or less belong to the next block.
    """
    indent = lines[start].indent
    end = start + 1
    last_child = start
    while end < len(lines):
        line = lines[end]
        if line.is_blank or line.is_comment:
            if line.is_comment and line.indent > indent:
                last_child = end
            end += 1
            continue
        if line.indent <= indent:
            break
        last_child = end
        end += 1
    return last_child + 1


def split_value(rest: str) -> tuple[str, str]:
    """Split the text after the colon into (value, trailing comment)."""
    match = re.match(r"^(\s*)(.*?)(\s+#.*)?$", rest)
    value = match.group(2)
    comment = match.group(3) or ""
    return value, comment


def find_key(lines: list[Line], path: tuple[str, ...], key: str) -> int | None:
    """Return the index of the kept line with this parent path and key, or None."""
    for i, line in enumerate(lines):
        if line.key == key and line.path == path and not line.delete:
            return i
    return None


def format_double(value: float) -> str:
    """Write a double so that it reads as a double, for example 60.0 and not 60."""
    text = f"{value:.10g}"
    if not any(c in text for c in ".eEn"):
        text += ".0"
    return text


def convert_ice_growth_step(lines: list[Line], messages: list[str]) -> None:
    """Replace the v1.2 ice growth step with the inner physics time step.

    Run this after the deletions and the renames: the ice growth line is already
    marked for deletion, and the outer time step has its current name.
    Insert the new key directly after the outer time step, or after the menu key.
    """
    growth_index = None
    for i, line in enumerate(lines):
        if line.key == ICE_GROWTH_KEY and line.path == ICE_GROWTH_PATH:
            growth_index = i
            break
    if growth_index is None:
        return
    growth_path = " -> ".join(ICE_GROWTH_PATH + (ICE_GROWTH_KEY,))
    inner_path = " -> ".join(TRANSPORT_PATH + (INNER_STEP_KEY,))
    growth_value, _ = split_value(lines[growth_index].rest)

    # If ice growth is off, the old step had no effect.
    switch_index = find_key(lines, ICE_GROWTH_PATH, ICE_GROWTH_SWITCH_KEY)
    if switch_index is not None:
        switch_value, _ = split_value(lines[switch_index].rest)
        if switch_value.strip().upper() in FALSE_STRINGS:
            messages.append(f"Did not convert '{growth_path}': ice growth is off.")
            return

    if find_key(lines, TRANSPORT_PATH, INNER_STEP_KEY) is not None:
        messages.append(
            f"Did not convert '{growth_path}' ({growth_value} min): "
            f"the file already sets '{inner_path}'."
        )
        return

    try:
        inner_s = float(growth_value) * 60.0
    except ValueError:
        messages.append(
            f"WARNING: '{growth_path}' is '{growth_value}', which is not a number. "
            f"Set '{inner_path}' yourself."
        )
        return
    if inner_s <= 0.0:
        messages.append(
            f"WARNING: '{growth_path}' is {growth_value}, but the inner time step must be "
            f"greater than 0. Set '{inner_path}' yourself."
        )
        return

    # The reader rejects an inner step that is longer than the outer step, so the
    # conversion needs the outer step of the file. Without it, write no inner step:
    # APCEMM then uses the value of defaults/input.yaml.
    outer_path = " -> ".join(TRANSPORT_PATH + (OUTER_STEP_KEY,))
    outer_index = find_key(lines, TRANSPORT_PATH, OUTER_STEP_KEY)
    if outer_index is None:
        messages.append(
            f"WARNING: the file does not set '{outer_path}', so '{growth_path}' is not converted. "
            f"APCEMM uses the default inner time step. To keep the step size of v1.2, "
            f"set '{inner_path}' to {format_double(inner_s)}."
        )
        return
    outer_value, _ = split_value(lines[outer_index].rest)
    try:
        outer_s = float(outer_value) * 60.0
    except ValueError:
        messages.append(
            f"WARNING: '{outer_path}' is '{outer_value}', which is not a number, "
            f"so '{growth_path}' is not converted. To keep the step size of v1.2, "
            f"set '{inner_path}' to {format_double(inner_s)}."
        )
        return

    if inner_s > outer_s:
        messages.append(
            f"WARNING: '{growth_path}' ({growth_value} min) is longer than '{outer_path}' "
            f"({outer_value} min), so it is not converted. APCEMM does not accept an inner "
            f"time step that is longer than the outer time step. "
            f"APCEMM uses the default inner time step. Set '{inner_path}' yourself."
        )
        return

    indent = lines[outer_index].indent
    new_line = Line(f"{' ' * indent}{INNER_STEP_KEY}: {format_double(inner_s)}\n")
    new_line.path = TRANSPORT_PATH
    lines.insert(outer_index + 1, new_line)
    messages.append(
        f"Converted '{growth_path}' ({growth_value} min) to "
        f"'{inner_path}' ({format_double(inner_s)} s), to keep the step size of v1.2. "
        f"To use the default inner time step instead, delete this line."
    )


def find_met_input_off(lines: list[Line]) -> list[str]:
    """Return the met input flags that are off, as full paths with their values."""
    off = []
    for line in lines:
        if line.key in MET_INPUT_FLAGS and line.path == MET_INPUT_PATH:
            value, _ = split_value(line.rest)
            if value.strip().upper() in FALSE_STRINGS:
                off.append(f"{' -> '.join(line.path + (line.key,))}: {value}")
    return off


def migrate(lines: list[Line]) -> list[str]:
    """Edit the lines in place and return the messages for the user."""
    messages = []
    annotate_paths(lines)

    # Keys that are present, per parent path. A rename must not create a duplicate.
    present: dict[tuple[str, ...], set[str]] = {}
    for line in lines:
        if line.key is not None:
            present.setdefault(line.path, set()).add(line.key)

    for i, line in enumerate(lines):
        if line.key is None or line.delete:
            continue
        full_path = " -> ".join(line.path + (line.key,))

        renames = RENAMED_KEYS.get(line.path, {})
        if line.key in renames:
            new_key = renames[line.key]
            if new_key in present.get(line.path, set()):
                messages.append(
                    f"Deleted '{full_path}': the file already sets '{new_key}'."
                )
            else:
                prefix = line.text[: line.indent]
                line.text = f"{prefix}{new_key}:{line.rest}\n"
                present[line.path].add(new_key)
                messages.append(f"Renamed '{full_path}' to '{new_key}'.")
                line.key = new_key
                continue

        if line.key in REMOVED_KEYS or line.key in renames:
            end = block_end(lines, i)
            for j in range(i, end):
                lines[j].delete = True
            # Delete the comments directly above the key (they describe it).
            # A comment that is less indented than the key is a section header: keep it.
            j = i - 1
            while (
                j >= 0
                and lines[j].is_comment
                and not lines[j].delete
                and lines[j].indent >= line.indent
            ):
                lines[j].delete = True
                j -= 1
            if line.key not in renames:
                messages.append(f"Deleted '{full_path}'.")
            continue

        forced = FORCED_VALUES.get(line.path, {})
        if line.key in forced:
            value, comment = split_value(line.rest)
            new_value = forced[line.key]
            if value != new_value:
                prefix = line.text[: line.indent]
                line.text = f"{prefix}{line.key}: {new_value}{comment}\n"
                messages.append(f"Set '{full_path}' to {new_value} (was {value}).")
                # Fix the comment above, if it tells the user to keep the option on.
                if (
                    i > 0
                    and lines[i - 1].is_comment
                    and lines[i - 1].text.strip() == "# Keep on"
                ):
                    prev = lines[i - 1]
                    prev.text = (
                        f"{prev.text[: prev.indent]}# Keep off, currently disabled\n"
                    )

        if line.key == TROPOPAUSE_KEY:
            value, _ = split_value(line.rest)
            try:
                pressure = float(value)
            except ValueError:
                pressure = None
            if pressure != TROPOPAUSE_OLD_HARDCODED_PA:
                messages.append(
                    f"WARNING: '{full_path}' is {value}. On v1.2 APCEMM ignored this value and "
                    f"used {TROPOPAUSE_OLD_HARDCODED_PA:.1e} Pa. Now APCEMM uses the value. "
                    f"Set it to {TROPOPAUSE_OLD_HARDCODED_PA:.1e} to be consistent with v1.2. "
                    "Contrail simulation results should not be affected by the value anyway."
                )

    convert_ice_growth_step(lines, messages)
    return messages


def drop_extra_blank_lines(texts: list[str]) -> list[str]:
    """Collapse the runs of blank lines that the deletions create."""
    out = []
    for text in texts:
        if text.strip() == "" and out and out[-1].strip() == "":
            continue
        out.append(text)
    return out


class ConversionError(Exception):
    """The script cannot convert the input file. It writes no output file."""


def run(input_file: Path, output_file: Path):
    with open(input_file, encoding="utf-8") as f:
        original = f.readlines()
    if original and not original[-1].endswith("\n"):
        original[-1] += "\n"

    lines = [Line(text) for text in original]
    annotate_paths(lines)
    met_input_off = find_met_input_off(lines)
    if met_input_off:
        flags = "\n".join(f"  {flag}" for flag in met_input_off)
        raise ConversionError(
            f"{input_file} cannot be converted. These flags turn off / overrode the met file:\n"
            f"{flags}\n"
            "In v1.2 APCEMM then used the constant temperature, RH, wind shear, moist layer, "
            "and lapse rate options. v1.3 removes these options and reads met data only from "
            "the met file."
        )
    messages = migrate(lines)

    kept = [line.text for line in lines if not line.delete]
    # Only collapse blank lines that the deletions created, not the ones of the user.
    had_double_blank = any(
        a.strip() == "" and b.strip() == "" for a, b in zip(original, original[1:])
    )
    if not had_double_blank:
        kept = drop_extra_blank_lines(kept)

    with open(output_file, "w", encoding="utf-8") as f:
        f.writelines(kept)

    for message in messages:
        print(message)
    if not messages:
        print("No change needed.")

    print(
        "\nFinal checks you need to do:\n"
        " -> Paths to engine/ambient conditions files might need to be replaced with '=DEFAULT='\n"
        " -> Check that 'Inner physics time step [s] (double)' is consistent with your previous input & new scheme\n"
    )


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print(__doc__)
        sys.exit(1)

    in_path = Path(sys.argv[1])
    out_path = Path(sys.argv[2])

    try:
        run(input_file=in_path, output_file=out_path)
    except ConversionError as error:
        print(error)
        print(f"Did not write {out_path}.")
        sys.exit(1)
    print(f"Saved to {out_path}")
