"""Plot libraries for the per-method sweeps (NJ, SNJ, bpart).

These live in ``src/`` because notebooks and ``src/runners/`` import them. They used to
sit in ``scripts/``, which forced ``scripts/README.md`` to open with a warning that seven
of its ``plot_*`` files were libraries rather than scripts -- and made ``src/`` depend on
``scripts/``. Each module keeps its ``__main__`` guard, so they still run standalone.

No re-exports: import the specific module (they pull matplotlib at import time).
"""
