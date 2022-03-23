import io


def nucleotide_complements_table():
    # A nice list of complements can be found at:
    # http://www.reverse-complement.com/ambiguity.html
    table = list(range(256))

    table[ord('A')] = "'T'"
    table[ord('a')] = "'t'"
    table[ord('C')] = "'G'"
    table[ord('c')] = "'g'"
    table[ord('G')] = "'C'"
    table[ord('g')] = "'c'"
    table[ord('T')] = "'A'"
    table[ord('t')] = "'a'"

    table[ord('U')] = "'A'"
    table[ord('u')] = "'a'"

    # R, purine (A, G) vs Y, pyrimidine (C, T)
    table[ord('R')] = "'Y'"
    table[ord('r')] = "'y'"
    table[ord('Y')] = "'R'"
    table[ord('y')] = "'r'"

    # K, keto (G, T) vs A, amino (A, C)
    table[ord('K')] = "'M'"
    table[ord('k')] = "'m'"
    table[ord('M')] = "'K'"
    table[ord('m')] = "'k'"

    # B, not A, vs V, not T
    table[ord('B')] = "'V'"
    table[ord('b')] = "'v'"
    table[ord('V')] = "'B'"
    table[ord('v')] = "'b'"

    # D, not C vs H, not G
    table[ord('D')] = "'H'"
    table[ord('d')] = "'h'"
    table[ord('H')] = "'D'"
    table[ord('h')] = "'d'"

    return table


def make_table(variable_name, table, row_size=16):
    out = io.StringIO()
    out.write(variable_name + ' = {\n    ')
    i = 0
    for i, literal in enumerate(table):
        if i % row_size == 0 and i != 0:
            out.write("\n    ")
        out.write(str(literal).rjust(3, " ") + ", ")
    out.write("\n")
    out.write("};\n")
    return out.getvalue()


def main():
    with open("src/dnaio/_conversions.h", "wt", encoding="utf-8") as out:
        out.write(make_table(
            "static const char NUCLEOTIDE_COMPLEMENTS[256]",
            nucleotide_complements_table())
        )
        out.write('\n')


if __name__ == "__main__":
    main()
