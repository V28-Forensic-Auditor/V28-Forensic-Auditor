# 🩸 V28 Forensic Auditor
**Independent Forensic Division. We don't predict, we falsify. (独立法医审计部。我们不预测，我们只证伪。)**

### 🔬 What is V28?
V28 Forensic Auditor is a specialized Python-based toolkit built to rigorously evaluate the validity of AI-generated biomolecular designs, particularly in the era of protein language models and diffusion-based generators.

Unlike black-box AI scoring methods, V28 employs strict, first-principles physics. It utilizes hard-coded physical axioms—including van der Waals radii, Euclidean distances, and chirality tensors—to mathematically prove whether a generated structure violates fundamental physical laws. Its core mission is to identify and falsify "physically impossible" artifacts that might otherwise slip through traditional validation methods.

The project is built on **Pure-StdLib** and **Zero-Dependency**, making it completely self-contained, reproducible, and easy to integrate into existing scientific pipelines.

### 🏆 Core Achievements
- **Extensive FPR Benchmarking:** Achieved a **3.26%** false-positive rate against **100** experimental RCSB PDB structures (0.00% on NMR).
- **Protein-Ligand Interception:** Detected a fatal **1.336 Å** steric clash in classic experimental data (e.g., 3PGK), proving physical impossibilities exist even in human-deposited data.
- **Official Corroboration:** [An audit issue](https://github.com/RosettaCommons/RFdiffusion/issues/458) submitted to the official RFdiffusion repository (RosettaCommons) regarding missing chain breaks was confirmed by core developers.
- **MCP Agent Ready:** Fully encapsulated as a Model Context Protocol (MCP) server, allowing LLMs (Claude/Cursor/Trae) to invoke V28 directly for automated physical QA.

### 📄 Full White Paper & Data
The complete v1.5 audit report, evidence chain, and adversarial test cases are permanently archived on Zenodo:
👉 **DOI:** [10.5281/zenodo.24276346](https://doi.org/10.5281/zenodo.24276346)

> 🏛️ **Official Citation & Verification:**
> - **Cite this work (Zenodo):** [10.5281/zenodo.24276346](https://doi.org/10.5281/zenodo.24276346)
> - **Author's Academic Profile (ORCID):** [0009-0002-4588-3006](https://orcid.org/0009-0002-4588-3006)
> 
> *⚠️ Disclaimer: This engine is a forensic validator, not a predictive model. Core probe source codes remain private (`.pyc` bytecode) for commercial third-party auditing. Please contact the author for independent physical compliance audits.*

**Contact for Commercial Deployment & MCP Integration:** v28.engine@outlook.com
