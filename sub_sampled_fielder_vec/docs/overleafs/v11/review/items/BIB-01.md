### [BIB/01] Three bibliography entries are never cited
- **Anchor:** references.bib:1
- **Source:** latex-paper-en:bibliography
- **Category:** unused-entry
- **Type:** citation
- **Severity:** editorial
- **Context:** The bibliography file, which the paper cites from throughout.
- **Item:** Three entries are defined and cited nowhere: BandeltDress1992, Griffing2012 and
  stone2009fiedler. Griffing belongs to the distance-based route, which v10 removed.
- **Affects:** none
- **Patch:** none
- **Status:** open
- **Verified:** yes -- verify_bib.py reports them and a grep for each key over sections/
  returns no \cite
