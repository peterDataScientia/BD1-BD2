# Observed QC outputs, 26 November 2026

Transcribed from user-provided CHPC terminal output; not a new execution.

- PROPKA 3.5.1; PDB2PQR 3.6.2; Python 3.9.6; urllib3 1.26.20; requests 2.32.5. OpenSSL 1.0.2k-fips necessitated the temporary urllib3 pin. MGLTools changed default Python to 2.5.6; explicit virtual-environment executable restored correct interpreter.
- Histidine side chains originally contained HD1 and HE2. Final states: BD1 A77 HID; BD2 B369 HIE, B388 HID, B396 HID, B437 HID. PDBQT N/NA assignments match neutral tautomers.
- BD2 chain B residues 331–348 and GLU460 unresolved; ASP459 OXT reconstructed for truncated model. BD1 GLU168 OXT preserved.
- Heavy atoms before/after PDB2PQR: BD1 1062/1062; BD2 905/906. No removed heavy atoms. Only added atom BD2 B459 OXT.
- Heavy-atom RMSD before/after preparation: BD1 0.097705 Å; BD2 0.105952 Å. Two atoms displaced >0.1 Å per domain (BD1 ASN116 ND2/OD1: 2.360072/2.137322 Å; BD2 GLN353 NE2/OE1: 2.383794/2.115866 Å). Ligand minimum distances before/after: 21.145/21.109 Å and 24.375/24.375 Å. Water minimum distances: 17.629/17.591 Å and 20.974/20.974 Å.
- PDB→PDBQT: 1062 and 906 matched heavy atoms, zero RMSD and maximum displacement; no added or removed heavy atoms.
- Final wet PDBQT: 1295+5=1300 and 1105+5=1110 protein/water records. Both water QC tests passed: protein records unchanged, water coordinates/charges/types unchanged, five waters each. Waters are rigid oxygen-only representations.
- Earlier wet baseline redocking: BD1 -9.700 kcal/mol, RMSD 1.045 Å; BD2 -9.273 kcal/mol, RMSD 1.238 Å. Revised pH 7.5 repeat redocking not yet confirmed submitted or completed.
