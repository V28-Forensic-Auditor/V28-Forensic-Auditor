# 🩸 V28 Forensic Auditor
**Independent Forensic Division. We don't predict, we falsify. (独立法医审计部。我们不预测，我们只证伪。)**

### 🔬 What is V28?
V28 is a **Pure-StdLib, Zero-Dependency, White-Box Audit Framework** designed to evaluate the physical and kinematic integrity of AI-generated biomolecular structures. 

Unlike black-box AI models, V28 uses hard-coded physical axioms (van der Waals radii, Euclidean distances, chirality tensors) to prove when a structural artifact is physically impossible.

### 🏆 Core Achievements
- **Extensive FPR Benchmarking:** Achieved a **3.26%** false-positive rate against **100** experimental RCSB PDB structures (0.00% on NMR).
- **Protein-Ligand Interception:** Detected a fatal **1.336 Å** steric clash in classic experimental data (e.g., 3PGK), proving physical impossibilities exist even in human-deposited data.
- **Official Corroboration:** An audit issue submitted to the official RFdiffusion repository (RosettaCommons) regarding missing chain breaks was confirmed by core developers.
- **MCP Agent Ready:** Fully encapsulated as a Model Context Protocol (MCP) server, allowing LLMs (Claude/Cursor/Trae) to invoke V28 directly for automated physical QA.

### 📄 Full White Paper & Data
The complete v1.5 audit report, evidence chain, and adversarial test cases are permanently archived on Zenodo:
👉 **DOI:** [10.5281/zenodo.24276346](https://doi.org/10.5281/zenodo.24276346)

*Note: The core probe source code is kept private (as `.pyc` bytecodes) for commercial third-party auditing purposes. For enterprise private deployment or independent auditing services, please contact the author.*
