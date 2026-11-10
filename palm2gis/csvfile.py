"""CSV files palm2gis writes, in the layout the reader's Excel expects.

advanced.csv_format chooses it:

  excel_cz (default)  ';' between fields and ',' as the decimal mark, UTF-8
                      with a byte-order mark: what Excel with Czech (and most
                      continental European) regional settings splits into
                      columns and reads as numbers on a double-click.
  standard            ',' between fields and '.' as the decimal mark, plain
                      UTF-8: Excel with English settings, pandas, R, QGIS.

A comma-separated file opened in a Czech Excel lands whole in column A,
and its '19.44' is read as text or a date, so a plain csv.writer is not
enough for the files people open by hand.
"""

import csv
from numbers import Integral, Real

FORMATS = ("excel_cz", "standard")


class CsvOut:
    """Context manager: `with CsvOut(path, fmt) as out: out.row([...])`.

    Floats are written with 6 significant digits and the format's decimal
    mark; integers and strings as they are.
    """

    def __init__(self, path, fmt="excel_cz"):
        if fmt not in FORMATS:
            raise ValueError(f"csv_format must be one of {', '.join(FORMATS)}, "
                             f"got {fmt!r}.")
        self.path = path
        self.excel = fmt == "excel_cz"
        self._f = None
        self._w = None

    def __enter__(self):
        self._f = open(self.path, "w", newline="",
                       encoding="utf-8-sig" if self.excel else "utf-8")
        self._w = csv.writer(self._f, delimiter=";" if self.excel else ",")
        return self

    def __exit__(self, *exc):
        self._f.close()
        return False

    def cell(self, v):
        if isinstance(v, bool) or v is None:
            return "" if v is None else str(v)
        if isinstance(v, Integral):
            return str(int(v))
        if isinstance(v, Real):
            s = "{:.6g}".format(float(v))
            return s.replace(".", ",") if self.excel else s
        return str(v)

    def row(self, values):
        self._w.writerow([self.cell(v) for v in values])


def read_rows(path):
    """Rows of a palm2gis CSV in either format, numbers as floats where a
    cell parses as one (for tests and scripts)."""
    with open(path, encoding="utf-8-sig", newline="") as f:
        text = f.read()
    first = text.splitlines()[0] if text else ""
    delim = ";" if first.count(";") > first.count(",") else ","
    rows = list(csv.reader(text.splitlines(), delimiter=delim))

    def num(c):
        try:
            return float(c.replace(",", ".") if delim == ";" else c)
        except ValueError:
            return c
    return [rows[0]] + [[num(c) for c in r] for r in rows[1:]], delim
