"""Test helper script that reads from stdin"""

import dnaio

if __name__ == "__main__":
    with dnaio.open("-") as f:
        records = list(f)
    print(len(records))
